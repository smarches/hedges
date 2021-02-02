#include <random>
#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <algorithm>
#include <type_traits>

#include "heapscheduler.h"

using Ullong = unsigned long long;
using Uchar = unsigned char;
// using GF4char = unsigned char; // semantically ACGT
using GF4reg = unsigned long long; // 8 bytes' worth of compressed GF4word values
// using Mbit = unsigned char; // semantically VARIABLE NUMBER of plaintext bits
// using VecMbit = std::vector<Mbit>;  // message bits unpacked to variable
// may want to use string if likely that small-string optimizations or etc. would be helpful
using GF4word = std::vector<Uchar>;  // string of ACGT
using VecUchar = std::vector<Uchar>; // aka std::basic_string<unsigned char>

template<class T>
class simple_matrix {
    static_assert(std::is_trivial_v<T>);
    public:
    std::size_t nrow,ncol;
    std::vector<T> data;
    simple_matrix() = delete;
    simple_matrix(size_t rows,size_t cols) : nrow(rows), ncol(cols) {
        data.reserve(nrow * ncol);
    }
    simple_matrix(size_t rows,size_t cols,T val) : nrow(rows), ncol(cols) {
        data(nrow * ncol,val);
    }
    T& operator()(size_t i,size_t j) {
        if(i >= nrow || j >= ncol) throw "simple_matrix index oob";
        return data[i*ncol + j];
    }
    T operator()(size_t i, size_t j) const {
        if(i >= nrow || j >= ncol) throw "simple_matrix index oob";
        return data[i*ncol + j];
    }
};

// the key step in the encoding of sequential (groups of 8) bytes
// baked-in parameters are 10, 24, 8 respectively
template<unsigned SEQBITS,unsigned HSALT,unsigned NPREV>
class digester {
    static constexpr Ullong seqnomask{(1LLU << SEQBITS) - 1};
    public:
    // this is basically a hash function so it could be replaced with any
    // suitable hashing functionality
    Ullong int64(Ullong u) const noexcept {
        Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
        v ^= v >> 21;
        v ^= v << 37; 
        v ^= v >> 4;
        v *= 4768777513237032717LL;
        v ^= v << 20; 
        v ^= v >> 41; 
        v ^= v << 5;
        return v;
    }
    std::int32_t digest(
        Ullong bits, 
        std::int32_t seq, 
        Ullong salt, 
        std::int32_t mod
    ) const noexcept {
        // get first 'seqbits' bits of seq and shift left by nprev
        const auto masked_seq = (static_cast<Ullong>(seq) & seqnomask) << NPREV;
        return int64(((masked_seq | bits) << HSALT) | salt) % mod;
    }
    unsigned seqbits() const noexcept { return SEQBITS; }
    unsigned hsalt()   const noexcept { return HSALT; }
    unsigned nprev()   const noexcept { return NPREV; }
};

unsigned int bytepopcount(Uchar byte) noexcept;

GF4word revcomp(const GF4word& arr);

unsigned countgc(const Ullong prev,const Ullong mask) noexcept;

unsigned count_gc(Ullong prev,Ullong mask) noexcept;

std::string dna_errors(const std::string& s,double sub_rate,double del_rate,double ins_rate);

double primerscore(const std::string& a,const std::string& b,double mispen,double gappen, double skewpen);

GF4word make_sense(const std::string& left_primer,const GF4word& codeword);

// simple helper class which generates values for the different 'patterns'
// rather than storing a (large) fixed-sized array and helplessly failing (and needing extra error cheking)
// if the string is larger, it seems much simpler to just generate them as needed
class pattarr {
    size_t prsize, psum;
    std::vector<std::uint8_t> pvec;
    public:
    pattarr(){}; // formality so it can be declared as a class data member
    pattarr(unsigned pattern,unsigned primer_size) : prsize(primer_size) {

        if(pattern < 1 || pattern > 6) throw("pattarr: pattern must be in [1,..,6]");
        if (pattern == 1) {      // rate 0.75
            pvec = {2,1};
        }
        else if (pattern == 2) { // rate 0.6
            pvec = {2,1,1,1,1};
        }
        else if (pattern == 3) { // rate 0.5
            pvec = {1};
        }
        else if (pattern == 4) { // rate 0.333
            pvec = {1,1,0};
        }
        else if (pattern == 5) { // rate 0.25
            pvec = {1,0};
        }
        else if (pattern == 6) { // rate 0.166
            pvec = {1,0,0};
        }
        psum = std::accumulate(std::cbegin(pvec),std::cend(pvec),0);
    }
    // returns value which would be at index i if we repeated pvec
    // enough times (including 'prsize' left-padding of zeros)
    unsigned operator[](size_t i) const noexcept {
        return i < prsize ? 0 : pvec[(i - prsize) % pvec.size()];
    }
    // cumulative sum of values as if accessed sequentially via the
    // index operator
    unsigned vbitsum(size_t i) const noexcept {
        const auto pvsize = pvec.size();
        if(i < prsize || pvsize == 0) return 0;
        auto n_use = i - prsize;
        size_t n_group = n_use / pvsize;
        unsigned rem = n_use % pvsize;
        return n_group * psum + pvec[rem];
    }
};

using hypovec = std::vector<Hypothesis<8,24>>;
class encoder {
    
    double reward;
    unsigned _nsalt, _pattern;
    static constexpr double rewards[] = {-0.035,-0.082,-0.127,-0.229,-0.265,-0.324};
    static constexpr Ullong acgtacgt{0x1b1b1b1b1b1b1b1bllu}; // "ACGTACGTACGTACGT"
    // 10, 24, 8 were the NSEQBITS, HSALT, NPREV global parameters
    digester<10,24,8> D;
public:
    pattarr Pat;
    GF4word leftprimer, rightprimer;
    std::vector<Ullong> primersalt;
    encoder() = delete;
    encoder(
        std::string& lprimer,std::string& rprimer,
        unsigned pattern = 1,unsigned nsalt = 24
    ) : _nsalt(nsalt), _pattern(pattern) {
        set_primers(lprimer,rprimer);
        set_coderate(pattern);
    }

    void set_coderate(unsigned pattern);
    void set_primers(const std::string& leftp,const std::string& rightp);
    size_t vbitlen(size_t num_msg_bits) const noexcept;
    size_t NSP() const noexcept { return vbitlen(_nsalt) + leftprimer.size();}
    VecUchar packvbits(const VecUchar& vbits, size_t msg_bits) const;
    VecUchar unpackvbits(const std::string& message, std::uint32_t len) const;
    size_t strand_len(size_t bytes) const noexcept;
    GF4word encode(
        std::string& message,
        size_t len,
        const unsigned DNAWINDOW = 12,
        const unsigned MAXGC = 8,
        const unsigned MINGC = 4,
        const unsigned MAXRUN = 4
    ) const;

    int search_heap(hypovec& hstack,std::vector<Uchar>& text,unsigned msg_bits,unsigned hlimit);
    std::string decode(GF4word data) const;
};

// default values are 8 and 24 resp.
template<unsigned NPREV, unsigned HSALT>
class Hypothesis {
    static constexpr Ullong prevmask{(1ULL << NPREV) - 1};
    static constexpr Ullong saltmask{(1ULL << HSALT) - 1};
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<double> runif(-1,1.);
    Ullong prevbits, salt, newsalt;
    public:
    size_t pred_ix, my_ix; // index of predecessor in hypostack
    int offset; // next char in message
    int seq; // my position in the decoded message (0,1,...)
    double score; // my -logprob score before update (what update??)
    Uchar messagebit; // last decoded up to now
    GF4reg prevcode;
    // TODO: cleanup with default/delete
    Hypothesis() {

        pred_ix = 0; // the root element is its own predecessor
        my_ix = 0;
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

    // does the same calculations as encode_C but just for the previous char rather than a whole message
    // TODO: make this a constructor overload?
    int init_from_predecessor(
        const GF4reg& codetext,
        const Hypothesis& hp,
        const Encoder& E,
        Uchar mbit, 
        int skew
    ) {

        const size_t LPRIMER{leftprimer.size()};
        const Ullong dnawinmask = (1ULL << 2 * GP::DNAWINDOW) - 1;

        pred_ix = hp.pred_ix;
        messagebit = mbit; // variable number
        seq = hp.seq + 1;  // always one longer than predecessor
        prevbits = hp.prevbits;
        salt = hp.salt;
        Ullong mysalt;
        const size_t nsp = E.NSP(), LPRIMER = E.leftprimer.size();
        
        if (seq < LPRIMER) {
            mysalt = E.primersalt[seq];
        }
        else if (seq < nsp) {
            mysalt = salt;
            newsalt = ((hp.newsalt << 1) & saltmask) ^ messagebit; // variable bits overlap, but that's ok with XOR
        }
        else if (seq == nsp) {
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
        int nbits = E.Pat[seq];
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
        if (dither > 0.) mypenalty += GP::dither * runif(rng);
        score = hp.score + mypenalty;
    }
};
