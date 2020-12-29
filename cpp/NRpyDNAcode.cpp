#include "nr3python.h"
#include "heapscheduler.h"
#include "dna_funcs.h"
#include <random>
#include <vector>
#include <array>
#include <string>
#include <cstdint>

//  this version 7 is version 6 with bug fixed in decode_c
//  this version 6 doesn't increment salt, but actually finds allowed output chars
//  this version 5 improves DNA constraints and does "fill" when codetext len is specified
//  this Version 4 adds DNA output constraints for GC balance and homopolymer runs
//  this Version 3 adds primers, check for coderate, and check for revcomp

constexpr char ThisVersion[] = "7.01";

using Ullong = unsigned long long;
using Uchar = unsigned char;
using GF4char = unsigned char; // semantically ACGT
// may want to use string if likely that small-string optimizations or etc. would be helpful
using GF4word = std::vector<Uchar>;  // semantically string of ACGT
using GF4reg = unsigned long long; // semantically a compressed GF4word
using Mbit = unsigned char; // semantically VARIABLE NUMBER of plaintext bits
using VecInt = std::vector<int>;
using VecMbit = std::vector<Mbit>;  // message bits unpacked to variable

// global parameters
namespace GP {

    // not normally user-adjustable
    constexpr unsigned int NPREV = 8; // number of hashed previous bits
    constexpr unsigned int NSEQBITS = 10; // number of hashed sequence number bits
    constexpr unsigned int HSALT = 24; // number of hashed bits of salt

    constexpr std::array<unsigned int,4> default_runtimes{ 24, 2500, 110'000, 1'000'000 };
    // user-adjustable parameters relating to runtime constraints
    unsigned int NSALT = default_runtimes[0]; // change salt after this many message bits (thus protecting them)
    unsigned int MAXSEQ = default_runtimes[1]; // maximum number of vbits in a message (one-time work in setcoderate() )
    unsigned int NSTACK = default_runtimes[2];  // initial size of list of hypotheses
    unsigned int HLIMIT = default_runtimes[3]; // limit on number of hypotheses tried before failure
    std::int32_t NSP = NSALT; // updated by setcoderate()
    std::int32_t lastpattnumber; // updated by setcoderate()

    constexpr std::array<unsigned int,4> default_dnas{ 12, 8, 4, 4 };
    // user-adjustable parameters for DNA-related encoding constraints
    std::int32_t DNAWINDOW = default_dnas[0]; // window in which DNA constraints imposed
    std::int32_t MAXGC = default_dnas[1]; // max GC in window
    std::int32_t MINGC = default_dnas[2]; // min GC in window
    std::int32_t MAXRUN = default_dnas[3]; // max length of homopolymers

    // rewards and penalties applied at each position when decoding
    constexpr std::array<double,5> default_scores{-0.13, 1., 1., 1., 0.};
    double reward = default_scores[0];
    double substitution = default_scores[1];
    double deletion = default_scores[2];
    double insertion = default_scores[3];
    double dither = default_scores[4];

    // related to the hypothesis/heap data structures:
    // std::int32_t nhypo{0};
    std::int32_t errcode{0};
    // size_t nnstak{0};
}

GF4word leftprimer, rightprimer;
std::vector<Ullong> primersalt;
 // contains number of bits in each vbit: 0, 1, or 2 (all initialized to 1)
 // set in setcoderate_C but accessed in various places
std::vector<std::uint8_t> pattarr;


Ran ran; // (11015);

HeapScheduler<double,std::int32_t> heap;
// set by traceback/traceback_fulldata (from decode_C); 
std::int32_t finaloffset, finalseq;
double finalscore;

// GF4reg dnawinmask((Ullong(1) << 2 * DNAWINDOW) - 1);
// GF4reg dnaoldmask((Ullong(1) << 2 * (DNAWINDOW-1)) - 1); // used to set oldest to "A"
// constexpr GF4reg acgtacgt{0x1b1b1b1b1b1b1b1bllu}; // "ACGTACGTACGTACGT" used for initialization

// TODO: turn digest into a simple class templated by NPREV, NSEQBITS, HSALT
// static Ranhash ranhash;


// std::int32_t digest(Ullong bits, std::int32_t seq, Ullong salt, std::int32_t mod) {
//     constexpr static Ullong seqnomask{(1LLU << GP::NSEQBITS) - 1};
//     const auto masked_seq = (static_cast<Ullong>(seq) & seqnomask) << GP::NPREV;
//     const auto res = ranhash.int64(((masked_seq | bits) << GP::HSALT) | salt) % mod;
// 	return res;
// }

// side effects: sets values of leftprimer, rightprimer, and primersalt
void findprimersalt(const char* leftpr, const char* rightpr) { // set salt to match a leftprimer
	const auto np = strlen(leftpr), mp = strlen(rightpr);
	constexpr char ACGT[] = "ACGTacgt";
    constexpr int num_try{100};
	std::array<unsigned,256> ACGTvalue{0};

	leftprimer.resize(np);
	rightprimer.resize(mp);
	primersalt.resize(np);

	for (size_t i = 0; i < 8; i++) ACGTvalue[ACGT[i]] = i % 4;
	for (size_t k = 0; k < np; k++) leftprimer[k] = ACGTvalue[leftpr[k]];
	for (size_t k = 0; k < mp; k++) rightprimer[k] = ACGTvalue[rightpr[k]];
	for (size_t k = 0; k < np; k++) {
		for (int i = 0; i < num_try; i++) { // try up to 100 times (what??)
			auto regout = digest(0ULL, k, i, 4);
			if (regout == leftprimer[k]) {
				primersalt[k] = i;
				break;
			}
		}
	}
}

// TODO: reparameterize (pass in patarr)
// global: accesses pattarr, GP::MAXSEQ
// nmb - number of message bits?
std::int32_t vbitlen(std::int32_t nmb) {  // how long is message in vbits?  (patarr must already be set)
	std::int32_t ksize{0}, nn{0};
    while( ksize++ ) {  // how many Mbits do we need?
		if (nn >= nmb) break;
		if (ksize >= GP::MAXSEQ) NRpyException("vbitlen: MAXSEQ too small");
		nn += pattarr[ksize];
	}
	return ksize;
}

static PyObject* getversion(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyObject(ThisVersion);
}

static PyObject* getparams(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyTuple(
		NRpyObject(GP::NSALT),
		NRpyObject(GP::MAXSEQ),
		NRpyObject(GP::NSTACK),
		NRpyObject(GP::HLIMIT),
		NULL);
}

static PyObject* restoreparams(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	GP::NSALT  = GP::default_runtimes[0];
	GP::MAXSEQ = GP::default_runtimes[1];
	GP::NSTACK = GP::default_runtimes[2];
	GP::HLIMIT = GP::default_runtimes[3];
	return NRpyObject(std::int32_t(0));
}

static PyObject* setparams(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 4) {
		NRpyException("setparams takes exactly 4 arguments");
		return NRpyObject(int(1));
	}
	GP::NSALT = NRpyInt(args[0]);
	GP::MAXSEQ = NRpyInt(args[1]);
	GP::NSTACK = NRpyInt(args[2]);
	GP::HLIMIT = NRpyInt(args[3]);
	return NRpyObject(std::int32_t(0));
}

static PyObject* getdnaconstraints(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyTuple(
		NRpyObject(GP::DNAWINDOW),
		NRpyObject(GP::MAXGC),
		NRpyObject(GP::MINGC),
		NRpyObject(GP::MAXRUN),
		NULL
    );
}

static PyObject* restorednaconstraints(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	GP::DNAWINDOW = GP::default_dnas[0];
	GP::MAXGC = GP::default_dnas[1];
	GP::MINGC = GP::default_dnas[2];
	GP::MAXRUN = GP::default_dnas[3];
	return NRpyObject(std::int32_t(0));
}

static PyObject* setdnaconstraints(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 4) {
		NRpyException("setdnaconstraints takes exactly 4 arguments");
		return NRpyObject(std::int32_t(1));
	}
	GP::DNAWINDOW = NRpyInt(args[0]);
	GP::MAXGC = NRpyInt(args[1]);
	GP::MINGC = NRpyInt(args[2]);
	GP::MAXRUN = NRpyInt(args[3]);
	return NRpyObject(std::int32_t(0));
}

static PyObject* getscores(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	return NRpyTuple(
		NRpyObject(GP::reward),
		NRpyObject(GP::substitution),
		NRpyObject(GP::deletion),
		NRpyObject(GP::insertion),
		NRpyObject(GP::dither),
		NULL
    );
}

static PyObject* restorescores(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	GP::reward = GP::default_scores[0];
	GP::substitution = GP::default_scores[1];
	GP::deletion = GP::default_scores[2];
	GP::insertion = GP::default_scores[3];
	GP::dither = GP::default_scores[4];
	return NRpyObject(std::int32_t(0));
}

static PyObject* setscores(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 5) {
		NRpyException("setscores takes exactly 5 arguments");
		return NRpyObject(std::int32_t(1));
	}
	GP::reward = NRpyDoub(args[0]);
	GP::substitution = NRpyDoub(args[1]);
	GP::deletion = NRpyDoub(args[2]);
	GP::insertion = NRpyDoub(args[3]);
	GP::dither = NRpyDoub(args[4]);
	return NRpyObject(std::int32_t(0));
}

static PyObject* minstrandlen(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	const auto nbytes = NRpyInt(args[0]);
	auto len = vbitlen(8 * nbytes) + rightprimer.size();
	return NRpyObject(len);
}

// given the 'seed' value, create a hashed value
// TODO: use std::hash
// is there a point of exporting this?
// static PyObject* hashint(PyObject *self, PyObject *pyargs) {
// 	NRpyArgs args(pyargs);
// 	if (args.size() != 1) {
// 		NRpyException("hashint takes exactly 1 argument");
// 		return NRpyObject(std::int32_t(1));
// 	}
// 	auto nn = NRpyInt(args[0]);
// 	auto hash = ranhash.int32(Ullong(nn));
// 	return NRpyObject(hash);
// }

// side effect: populates pattarr, sets GP::reward
// globals: accesses leftprimer, GP::NSALT, GP::MAXSEQ
void setcoderate_C(int pattnumber, const char* left_primer, const char* right_primer) { // some standard patterns
	findprimersalt(left_primer,right_primer);
    std::vector<std::uint8_t> pattern;
    static double rewards[] = {-0.035,-0.082,-0.127,-0.229,-0.265,-0.324};
    GP::reward = rewards[pattnumber];
    
	if (pattnumber == 1) { // rate 0.75
        pattern = {2,1};
	}
	else if (pattnumber == 2) { // rate 0.6
        pattern = {2,1,1,1,1};
	}
	else if (pattnumber == 3) { // rate 0.5
        pattern = {1};
	}
	else if (pattnumber == 4) { // rate 0.333
        pattern = {1,1,0};
	}
	else if (pattnumber == 5) { // rate 0.25
        pattern = {1,0};
	}
	else if (pattnumber == 6) { // rate 0.166
        pattern = {1,0,0};
	}
	pattarr.assign(GP::MAXSEQ + 2, 1);
	const auto n_pattern{pattern.size()}, LPRIMER{leftprimer.size()};
    std::fill(std::begin(pattarr),std::begin(pattarr) + LPRIMER,0);
	for (unsigned int i = LPRIMER; i < GP::MAXSEQ; i++) pattarr[i] = pattern[i % n_pattern];
	const auto vsalt = vbitlen(GP::NSALT);
	GP::NSP = vsalt + LPRIMER;
    GP::lastpattnumber = pattnumber;
}

// globals: sets GP::errcode
// accessses heap, MAXSEQ, pattarr
// returns 'nfinal', the # of total hypotheses considered
// TODO: turn into method
int shoveltheheap(const GF4word& text, int nmessbits,size_t hlimit) {
	// given the heap, keep processing it until offset limit, hypothesis limit, or an error is reached
    const auto seqmax = vbitlen(nmessbits);
    const size_t limit{text.size()};
	int qq, qqmax = -1, ofmax = -1;
    constexpr std::array<int,3> skews{0,-1,1};
	GP::errcode = 0;
	while (true) {
        const auto res = heap.pop();
        double currscore = res.first;
        qq = res.second;
		const auto hp = hypostack[qq];
		auto seq = hp.seq;
		
		if (seq > GP::MAXSEQ) NRpyException("shoveltheheap: MAXSEQ too small");
		if (hp.offset > ofmax) { // keep track of farthest gotten to
			ofmax = hp.offset;
			qqmax = qq;
		}
		if (currscore > 1.e10) break; // heap is empty (TODO: change this)
		if (hp.offset >= limit-1) break; // errcode 0 (nominal success)
		if (nmessbits > 0 && seq >= seqmax-1) break; // ditto when no. of message bits specified
        if (hypostack.size() > hlimit) { 
			GP::errcode = 2;
			return qqmax;
		}

		const int nguess = 1 << pattarr[seq + 1]; // i.e., 1, 2, or 4
        for(auto skew : skews) {
            for (Uchar mbit = 0; mbit < nguess; mbit++) {
                Hypothesis h();
                if(h.init_from_predecessor(text, qq, mbit, skew)) {
                    hypostack.push_back(h);
                    heap.push(h.score,hypostack.size());
                }
            }
        }
	}
	return qq; // final position
}

static PyObject* setcoderate(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 3) {
		NRpyException("setcoderate takes exactly 3 arguments");
		return NRpyObject(int(1));
	}
	const auto pattnumber = NRpyInt(args[0]);
	if (pattnumber < 1 || pattnumber > 6) {
		NRpyException("setcoderate arg must be in range 1 to 6");
		return NRpyObject(std::int32_t(1));
	}
	const char *leftpr = NRpyCharP(args[1]);
	const char *rightpr = NRpyCharP(args[2]);
	setcoderate_C(pattnumber, leftpr, rightpr);
	return NRpyObject(std::int32_t(0));
}

// globals: accesses rightprimer, pattarr
// n should be length of message
VecMbit unpackvbits(const char *message, size_t n, std::int32_t len) {
	int nmb{8*n};
	const auto ksize = std::max(vbitlen(nmb), len - rightprimer.size()); // aim for codetext of length len if possible
	VecMbit ans(ksize,0);
	size_t i = j = 0;
	for (size_t k = 0; k < ksize; k++) {
		for (int k1 = 0; k1 < pattarr[k]; k1++) {
			Uchar bit = i < n ? (message[i] >> (7 - j++)) & 1 : 0;
			if (j == 8) { j = 0; ++i; }
			ans[k] = (ans[k] << 1) | bit;
		}
	}
	return ans;
}

// globals: accesses pattarr
VecUchar packvbits(const VecMbit& vbits, int nmessbits) {
	
    const auto ksize = vbits.size();
	if (ksize > GP::MAXSEQ) throw("packvbits: MAXSEQ too small");
    // number of bits
    const auto nn = std::accumulate(std::begin(pattarr),std::end(pattarr),0);
	// for (size_t k = 0; k < ksize; k++) nn += pattarr[k]; 
	nn = std::min(nn, nmessbits); // no more than the specified number of bits
	nn = (nn + 7) / 8; // number of bytes
	VecUchar ans(nn,0);
	int i = j = 0;
	for (size_t k = 0; k < ksize; k++) {
		for (auto k1 = pattarr[k] - 1; k1 >= 0; k1--) {
			Uchar bit = (vbits[k] >> k1) & 1;
			ans[i] = ans[i] | (bit << (7-j++));
			if (j == 8) {
				j = 0;
				if (++i == nn) break; 
			}
		}
		if (i == nn) break;
	}
	return ans;
}

// globals: accesses GP::DNAWINDOW,GP::MAXSEQ,leftprimer,rightprimer,primersalt,
// GP::NSP,pattarr
GF4word encode_C(const char *message, size_t n, int len=0) { // dnac
	
    static constexpr Ullong prevmask{(Ullong(1) << GP::NPREV) - 1};
    static constexpr Ullong saltmask{(Ullong(1) << GP::HSALT) - 1};
    static constexpr Ullong acgtacgt{0x1b1b1b1b1b1b1b1bllu}; // "ACGTACGTACGTACGT"
    const Ullong dnawinmask = (1ULL << 2 * GP::DNAWINDOW) - 1;
	GF4word vbits = unpackvbits(message, n, len);
	const int nm{vbits.size()}; // number of variable bits encoded
	if (nm > GP::MAXSEQ) throw("encode: MAXSEQ too small");

    const size_t LPRIMER{leftprimer.size()};
	GF4word codetext(nm + rightprimer.size());
	
	Ullong prevbits = 0, salt = 0, newsalt = 0; 
	GF4reg prevcode = acgtacgt; // initialize with no runs and balanced cg
	for (size_t k = 0; k < nm; k++) { // on decoding, k is called seq
		Mbit messagebit = vbits[k];
		if (k < LPRIMER) {
			salt = primersalt[k];
		}
		else if (k < GP::NSP) {
			salt = 0;
			newsalt = ((newsalt << 1) & saltmask) ^ messagebit;
		}
		else if (k == GP::NSP) {
			salt = newsalt; // time to update the salt
		}
        if(k < LPRIMER) {
            auto regout = digest(prevbits, k, salt, 4);
            codetext[k] = (regout + messagebit) % 4;
        } else {
            const auto dnac_ok = dnacallowed(
                prevcode,GP::DNAWINDOW,GP::MAXGC,GP::MINGC,GP::MAXRUN
            );
            int mod = dnac_ok.size();
            auto regout = digest(prevbits, k, salt, mod);
		    regout = (regout + messagebit) % mod;
            codetext[k] = dnac_ok[regout];
        }
		auto nbits = pattarr[k];
		prevbits = ((prevbits << nbits) & prevmask) | messagebit; // variable number
		prevcode = ((prevcode << 2) | codetext[k]) & dnawinmask;  // shift code by 1
	}
    std::copy(std::cbegin(rightprimer),std::cend(rightprimer),std::begin(codetext) + nm);
	return codetext;
}

GF4word encode_C(VecUchar& message, int len=0) {
	return encode_C((char*)(&message[0]), message,len);
}
GF4word encode_C(char *message, int len=0) {
	return encode_C(message, strlen(message),len);
}

struct Hypothesis {
	int predi; // index of predecessor in hypostack
	int offset; // next char in message
	int seq; // my position in the decoded message (0,1,...)
	double score; // my -logprob score before update
	Mbit messagebit; // last decoded up to now
	Ullong prevbits, salt, newsalt;
	GF4reg prevcode;
	
    // TODO: cleanup with default/delete
	Hypothesis() {}
	Hypothesis(int) {} // so that can cast from zero in NRvector constructor

    // globals: accesses hypostack,NPREV,HSALT,leftprimer,GP::DNAWINDOW,
    // pattarr,deletion,reward,substitution,insertion,dither
    // does the same calculations as encode_C but just for the previous char rather than a whole message
	int init_from_predecessor(const GF4reg& codetext,int pred, Mbit mbit, int skew) {
		
        static constexpr Ullong prevmask{(1ULL << GP::NPREV) - 1};
        static constexpr Ullong saltmask{(1ULL << GP::HSALT) - 1};
        if(hypostack.size() <= pred) throw("predecessor does not exist on the hypothesis stack");
        const size_t LPRIMER{leftprimer.size()};
        const Ullong dnawinmask = (1ULL << 2 * GP::DNAWINDOW) - 1;
		
		const Hypothesis hp = hypostack[pred]; // predecessor
		predi = pred;
		messagebit = mbit; // variable number
		seq = hp.seq + 1; // always one longer than predecessor
		if (seq > GP::MAXSEQ) throw("init_from_predecessor: MAXSEQ too small");
		prevbits = hp.prevbits;
		salt = hp.salt;
		Ullong mysalt;
		if (seq < LPRIMER) {
			mysalt = primersalt[seq];
		}
		else if (seq < GP::NSP) {
			mysalt = salt;
			newsalt = ((hp.newsalt << 1) & saltmask) ^ messagebit; // variable bits overlap, but that's ok with XOR
		}
		else if (seq == GP::NSP) {
			mysalt = salt = hp.newsalt; // time to update the salt
		}
		else mysalt = salt;
		offset = hp.offset + 1 + skew;
		if (offset >= codetext.size()) return 0; // i.e., false
		// calculate predicted message under this hypothesis
		prevcode = hp.prevcode;
        int regout;
        if(seq < LPRIMER) {
            regout = digest(prevbits, seq, mysalt, 4);    
            regout = (regout + Uchar(messagebit)) % 4;
        } else {
            const auto dnac_ok = dnacallowed(
                prevcode,GP::DNAWINDOW,GP::MAXGC,GP::MINGC,GP::MAXRUN
            );
            int mod = dnac_ok.size();
            regout = digest(prevbits, seq, mysalt, mod);
            regout = (regout + Uchar(messagebit)) % mod;
            regout = dnac_ok[regout];
        }
		int nbits = pattarr[seq];
		prevbits = ((hp.prevbits << nbits) & prevmask) | messagebit; // variable number
		prevcode = ((prevcode << 2) | regout) & dnawinmask;
		double mypenalty;
		// compare to observed message and score
		if (skew < 0) { //deletion
			mypenalty = GP::deletion;			
		} else {
			bool ismatch = regout == codetext[offset];  // the only place where a check is possible!
            float adj = ismatch ? GP::reward : GP::substitution;
            mypenalty = adj + (skew == 0) * GP::insertion;
		}
		if (dither > 0.) mypenalty += GP::dither * (2.*ran.doub() - 1.);
		score = hp.score + mypenalty;
		return 1; // i.e., true
	}
	void init_root() {
		predi = -1;
		offset = -1;
		seq = -1;
		messagebit = 0; // not really a message bit
		prevbits = 0;
		score = 0.;
		salt = 0;
		newsalt = 0;
        constexpr GF4reg acgtacgt{0x1b1b1b1b1b1b1b1bllu}; // "ACGTACGTACGTACGT"
		prevcode = acgtacgt;
	}	
	int myid() const { return static_cast<int>(this - &hypostack); } // my index in hypostack
};

std::vector<Hypothesis> hypostack;

// globals: heap hypostack,GP::NSTACK
// TODO: turn into a destructor
void release() { // give back heap and hypostack memory
	heap.reinit();
	hypostack.resize(GP::NSTACK);
}

//globals: heap, hypostack
// TODO: turn into constructor
void init_heap_and_stack() {
	
	hypostack.clear();
    Hypothesis h;
    h.init_root();
	hypostack.emplace_back(h);
	heap.rewind();
	heap.push(1.e10, 0);
}


static PyObject* encode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	
	if (args.size() > 2) {
		NRpyException("encode takes 1 or 2 arguments");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("encode requires array with dtype=uint8 \n");
	VecUchar message(args[0]);
	int len =  args.size() > 1 ? NRpyInt(args[1]) : 0;
	const auto codetext = encode_C(message,len);
	return NRpyObject(codetext);
}

static PyObject* encodestring(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);

	if (args.size() > 2) {
		NRpyException("encodestring takes 1 or 2 arguments");
		return NRpyObject(0); // formerly NULL
	}
	const char *message = NRpyCharP(args[0]);
	int len =  args.size() > 1 ? NRpyInt(args[1]) : 0;
	const auto codetext = encode_C(message,len);
	return NRpyObject(codetext);
}

// accesses hypostack
// TODO: turn into method
VecMbit traceback(int nfinal) {
	int kk = 0, q = nfinal;
	while ((q = hypostack[q].predi) > 0) ++kk; // get length of chain (this is the seq field!?)
    // update global params (thjs need not be done here?)
	finalscore = hypostack[nfinal].score;
	finaloffset = hypostack[nfinal].offset;
	finalseq = hypostack[nfinal].seq;
	// todo - refactor so we loop through once
	VecMbit ans(kk + 1,0);
    // populate the answer
	q = nfinal;
	ans[kk--] = hypostack[q].messagebit;
	while ((q = hypostack[q].predi) > 0) {
		ans[kk] = hypostack[q].messagebit;
		--kk;
	}
	return ans;
}

// global containers for fulldata
VecInt allseq;
VecInt allnhypo;
VecInt alloffset;
std::vector<double> allscore;
VecInt allpredi;
VecUchar allmessagebit;
VecInt allprevbits;
VecInt allsalt;
VecInt allnewsalt;

// TODO: turn into method
void traceback_fulldata(const std::vector<Hypothesis>& hypo_stack,int nfinal) {
// TODO: questionable! messagebit might be 0, 1 or 2 bits.  how are you supposed to know?
// see packvbits()
	int kk = 0, q = nfinal;
	while ((q = hypo_stack[q].predi) > 0) ++kk; // get length of chain
	finalscore = hypo_stack[nfinal].score;
	finaloffset = hypo_stack[nfinal].offset;
	finalseq = hypo_stack[nfinal].seq;
	allseq.resize(kk + 1);
	alloffset.resize(kk + 1);
	allscore.resize(kk + 1);
	allnhypo.resize(kk + 1);
	allpredi.resize(kk + 1);
	allmessagebit.resize(kk + 1);
	allprevbits.resize(kk + 1);
	allsalt.resize(kk + 1);
	allnewsalt.resize(kk + 1);

	q = nfinal;
	int k = kk;
	allseq[k] = hypo_stack[q].seq;
	alloffset[k] = hypo_stack[q].offset;
	allscore[k] = hypo_stack[q].score;
	allnhypo[k] = q;
	allpredi[k] = hypo_stack[q].predi;
	allmessagebit[k] = hypo_stack[q].messagebit;
	allprevbits[k] = static_cast<int>(hypo_stack[q].prevbits); // only returning 32 (or 31) bits of these
	allsalt[k] = static_cast<int>(hypo_stack[q].salt);
	allnewsalt[k] = static_cast<int>(hypo_stack[q].newsalt);
	--k;
	while ((q = hypo_stack[q].predi) > 0) {
		allseq[k] = hypo_stack[q].seq;
		alloffset[k] = hypo_stack[q].offset;
		allscore[k] = hypo_stack[q].score;
		allnhypo[k] = q;
		allpredi[k] = hypo_stack[q].predi;
		allmessagebit[k] = hypo_stack[q].messagebit;
		allprevbits[k] = static_cast<int>(hypo_stack[q].prevbits); // only returning 32 (or 31) bits of these
		allsalt[k] = static_cast<int>(hypo_stack[q].salt);
		allnewsalt[k] = static_cast<int>(hypo_stack[q].newsalt);
		--k;
	}
}

int decode_prep(GF4word& text,int nmessbits,size_t hlimit = GP::HLIMIT) {
	init_heap_and_stack();
    // returning the former 'nfinal' global var
    return shoveltheheap(text, nmessbits, hlimit);
}

VecUchar decode_C(GF4word& codetext, int nmessbits,size_t hlimit) {
    auto nfinal = decode_prep(codetext,nmessbits,hlimit); // should nmessbits also be zero here?
	VecMbit trba = traceback(nfinal);
	return packvbits(trba, nmessbits); // truncate only at the end
}

// loops throuh the different predefined 'patterns' to get coderates @ the given parameter values
// TODO: refactor so that the funcs using GP::MAXSEQ, lastpattnumber take them as parameters,
// since this function should just forward args
std::vector<int> gethowfar(int hlimit, int maxseq, GF4word &codetext, const char* leftpr, const char* rightpr) {
	std::array<int,6> ans{0}; // pattern 0 is not defined
	int MAXSEQ_save = GP::MAXSEQ, pattno_save = GP::lastpattnumber;
	GP::MAXSEQ = maxseq;
	for (int i = 0; i < 6; i++) {
		setcoderate_C(i + 1, leftpr, rightpr);
		auto dc = decode_C(codetext,0,hlimit);
		ans[i] = finaloffset;
	}
	// restore the globals
	GP::MAXSEQ = MAXSEQ_save;
    // its a bit odd since we 'restore' the old coderate but with these new primers...
	setcoderate_C(pattno_save, leftpr, rightpr);
	return ans;
}

void decode_fulldata_C(GF4word& codetext) {
    auto nfinal = decode_prep(codetext,0);
    // sets a bunch of vectors
	traceback_fulldata(hypostack,nfinal);
}

static PyObject* decode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	int nmessbits;
	if (args.size() == 1) {
		nmessbits = 0;
	}
	else if (args.size() == 2) {
		nmessbits = NRpyInt(args[1]);
	}
	else {
		NRpyException("decode takes 1 or 2 arguments only");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("decode requires array with dtype=uint8 \n");
	GF4word codetext(args[0]);
	VecUchar plaintext = decode_C(codetext,nmessbits,GP::HLIMIT);
	return NRpyTuple(
		NRpyObject(GP::errcode),
		NRpyObject(plaintext),
		NRpyObject(hypostack.size()),
		NRpyObject(finalscore),
		NRpyObject(finaloffset),
		NRpyObject(finalseq),
		NULL
    );
}

static PyObject* decode_fulldata(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("decode takes 1 argument only");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("decode requires array with dtype=uint8 \n");
	GF4word codetext(args[0]);
	decode_fulldata_C(codetext);
	VecUchar t_allmessagebit(allmessagebit);
	VecInt t_allseq(allseq), t_alloffset(alloffset), t_allpredi(allpredi), t_allprevbits(allprevbits),
	  t_allsalt(allsalt), t_allnewsalt(allnewsalt), t_allnhypo(allnhypo);
	VecDoub t_allscore(allscore);
    // must return copies, because Python gets control of its contents!
	return NRpyTuple(
		NRpyObject(GP::errcode),
		NRpyObject(hypostack.size()),
		NRpyObject(t_allmessagebit),
		NRpyObject(t_allseq),
		NRpyObject(t_alloffset),
		NRpyObject(t_allscore),
		NRpyObject(t_allnhypo),
		NRpyObject(t_allpredi),
		NRpyObject(t_allprevbits),
		NRpyObject(t_allsalt),
		NRpyObject(t_allnewsalt),
		NULL
	);
}

// changed to return the reverse-complement rather than modifying it in-place
static PyObject* revcomp(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("revcomp takes exactly 1 argument");
		return NRpyObject(0);
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("revcomp requires array with dtype=uint8 \n");
	GF4word arr(args[0]);
    auto ans = revcomp(arr);
	return NRpyObject(arr);
}

static PyObject* tryallcoderates(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 5) {
		NRpyException("tryallcoderates takes exactly 5 arguments");
		return NRpyObject(0);
	}
	int hlimit = NRpyInt(args[0]);
	int maxseq = NRpyInt(args[1]);
	if (PyArray_TYPE(args[2]) != PyArray_UBYTE) NRpyException("tryallcoderates requires array with dtype=uint8 \n");
	GF4word codetext(args[2]);
	const char *leftpr = NRpyCharP(args[3]);
	const char *rightpr = NRpyCharP(args[4]);
	VecInt maxoffsets = gethowfar(hlimit, maxseq, codetext, leftpr, rightpr);
	return NRpyObject(maxoffsets);
}

static PyObject* createerrors(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 4) {
		NRpyException("createerrors takes exactly 4 arguments");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) NRpyException("createrrors requires array with dtype=uint8 \n");
	const std::string codetext(args[0]);
	double srate = NRpyDoub(args[1]); 
	double drate = NRpyDoub(args[2]); 
	double irate = NRpyDoub(args[3]);
    const auto ans = dna_errors(srate,drate,irate);
	return NRpyObject(ans);
}

static PyObject* makegoodsense(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 2) {
		NRpyException("makegoodsense takes exactly 2 arguments");
		return NRpyObject(0);
	}
	if (PyArray_TYPE(args[1]) != PyArray_UBYTE) NRpyException("makegoodsense requires array with dtype=uint8 \n");

	const std::string left_primer(NRpyCharP(args[0]));
	const GF4word codeword(args[1]);
	GF4word newcodeword = make_sense(left_primer,codeword);
	return NRpyObject(newcodeword);
}

// this should probably happen automatically?
static PyObject* releaseall(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	release();
	allseq.resize(0);
	alloffset.resize(0);
	allscore.resize(0);
	allpredi.resize(0);
	allmessagebit.resize(0);
	allprevbits.resize(0);
	allsalt.resize(0);
	allnewsalt.resize(0);
	return NRpyObject(int(0));
}

// boilerplate 
static PyMethodDef NRpyDNAcode_methods[] = {
	{ "getversion", getversion, METH_VARARGS,
	"version = getversion()\n get version number as a float" },
	{ "minstrandlen", minstrandlen, METH_VARARGS,
	"minstrandlength = minstrandlen(nbytes)\n get min length of DNA to encode nbytes (then use longer!)" },
	{ "getparams", getparams, METH_VARARGS,
	"(NSALT, MAXSEQ, NSTACK, HLIMIT) = getparams()\n get current values of parameters" },
	{ "restoreparams", restoreparams, METH_VARARGS,
	"restoreparams()\n restore parameters to default values\nNB! must be followed by setcoderate" },
	{ "setparams", setparams, METH_VARARGS,
	"errorcode = setparams(NSALT, MAXSEQ, NSTACK, HLIMIT)\n set new parameter values\nNB! must be followed by setcoderate" },
	{ "getdnaconstraints", getdnaconstraints, METH_VARARGS,
	"(DNAWINDOW, MAXGC, MINGC, MAXRUN) = getdnaconstraints()\n get current values of DNA constraints" },
	{ "restorednaconstraints", restorednaconstraints, METH_VARARGS,
	"restorednaconstraints()\n restore DNA constraints to default values" },
	{ "setdnaconstraints", setdnaconstraints, METH_VARARGS,
	"errorcode = setdnaconstraints(DNAWINDOW, MAXGC, MINGC, MAXRUN)\n set new DNA constraint values\nDNAWINDOW=0 for no constraints" },	
	{ "getscores", getscores, METH_VARARGS,
	"(reward,substitution,deletion,insertion,dither) = getscores()\n get current scoring parameters" },
	{ "restorescores", restorescores, METH_VARARGS,
	"restorescores()\n restore scoring parameters to default values" },
	{ "setscores", setscores, METH_VARARGS,
	"errorcode = setscores(reward,substitution,deletion,insertion,dither)\n set new scoring parameters" },
	{ "setcoderate", setcoderate, METH_VARARGS,
	"errorcode = setcoderate(number, leftprimer, rightprimer)\n\
	 set coderate to one of six values for number=1..6 (0.75, 0.6, 0.5, 0.333, 0.25, 0.166)" },
	{ "encode", encode, METH_VARARGS,
	"int8_dna_array = encode(int8_message_array [, strandlen])\n encode a message with runout to strandlen" },
	{ "encodestring", encodestring, METH_VARARGS,
	"int8_dna_array = encodestring(message_as_string)\n encode a message" },
	{ "decode", decode, METH_VARARGS,
	"(errcode, int8_message_array, nhypo, score, offset, seq) = decode(int8_dna_array[, nmessbits])\n\
	decode a message optionally limited to nmessbits message bits" },
	{ "tryallcoderates", tryallcoderates, METH_VARARGS,
	"maxoffsets = tryallcoderates(hlimit, maxseq, int8_dna_array, leftprimer, rightprimer)\n\
	maxoffsets[i] is maximum offset achieved in trying coderate i (in 1..6) limited by hlimit" },
	{ "decode_fulldata", decode_fulldata, METH_VARARGS,
	"(errcode,nhypo,messagebit,seq,offset,score,hypo,predi,prevbits,salt,newsalt) =\n\
	decode_fulldata(codetext[, nmessbits])" },
	{ "createerrors", createerrors, METH_VARARGS,
	"new_int8_dna_array = createerrors(int8_dna_array, subrate, delrate, insrate)\n\
	create Poisson random errors at specified rates"},
	{ "releaseall", releaseall, METH_VARARGS,
	"errcode = releaseall()\n release memory grabbed by decode_fulldata"},
	{ "revcomp", revcomp, METH_VARARGS,
	"revcomp(int8_dna_array)\n reverse-complement a dna array in place"},
	{ "makegoodsense", makegoodsense, METH_VARARGS,
	"new_int8_dna_array = makegoodsense(left_primer, int8_dna_array)\n\
	return array or its reverse-complement, whichever agrees best with left_primer"},
	{ "hashint", hashint, METH_VARARGS,
	"hashedint = hashint(int)\n hash an integer by same algorithm as used throughout"},
	{ NULL, NULL, 0, NULL }
};
PyMODINIT_FUNC initNRpyDNAcode(void) { // N.B. must rename to agree with module name
	import_array();
	Py_InitModule("NRpyDNAcode", NRpyDNAcode_methods);  // N.B. must rename first arg, not second
}

