#include <type_traits>
#include <iostream>
#include "dna_funcs.h"

// not being used, but might someday!
// see https://stackoverflow.com/questions/9949935/calculate-number-of-bits-set-in-byte#25808559
unsigned int bytepopcount(Uchar byte) noexcept {
    static constexpr Uchar NIBBLE_LOOKUP[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
    return NIBBLE_LOOKUP[byte & 0x0F] + NIBBLE_LOOKUP[byte >> 4];
}

// this simple helper is not exported in the header
std::uint8_t char2i(char c) noexcept {
    switch(c) {
        case 'a':
        case 'A':
            return 0;
        case 'c':
        case 'C':
            return 1;
        case 'g':
        case 'G':
            return 2;
        case 't':
        case 'T':
            return 3;
        default:
            // should we map to zero or do something else?
            return 4;
    }
}

GF4word charvec_from_dna(const std::string& s) {
    GF4word res;
    res.reserve(s.size());
    for(const auto c: s) {
        res.push_back(char2i(c));
    };
    return res;
}

// converter from vector<Uchar> to std::string
std::string string_from_charvec(const GF4word& v) {
    static constexpr char ACGT[] = "ACGT";
    std::string ans;
    ans.reserve(v.size());
    for(const auto c: v ) {
        ans.push_back( c > 3 ? 'N' : ACGT[c]);
    };
    return ans;
}

// reverse complement for GF4word (aka vector of bytes)
GF4word revcomp(const GF4word& arr) {
    GF4word ans = arr;
    constexpr Uchar TGCA[] = { 3,2,1,0 };
    std::reverse(std::begin(ans),std::end(ans));
    std::for_each(std::begin(ans),std::end(ans),[TGCA](auto& e){
        e = e > 3 ? e : TGCA[e];
    });
    return ans;
}

// this could be constexpr but we'd never use it as such
// note, this only produces the correct answer if ALL input bits are
// in {0,1,2,3}!
unsigned countgc(const Ullong prev,const Ullong mask) noexcept {
        Ullong reg = prev & mask;
        // makes ones for GC, zeros for AT
        reg = (reg ^ (reg >> 1)) & 0x5555555555555555ull;
        reg -= ((reg >> 1) & 0x5555555555555555ull);
        reg = (reg & 0x3333333333333333ull) + (reg >> 2 & 0x3333333333333333ull);
        return ((reg + (reg >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

/* count # of G/C when the input 8 bytes are viewed as 1-byte chars,
 * remembering that A,C,G,T are represented as 0,1,2,3 here
 * the order that the bytes are iterated over could depend on the CPU being
 * little endian or big endian, but fortunately order of traversal doesn't matter
 * for this count.
 */
unsigned count_gc(Ullong prev,Ullong mask) noexcept {
    Ullong prev_mask = prev & mask;
    unsigned res = 0;
    for(unsigned i=0; i<8;++i) {
        Ullong imask = (prev_mask >> 8*i) & 0xff;
        res += imask == 1 || imask == 2;
    }
    return res;
}

// length of return value is used as the 'mod' value in digest()
// values in return val are the allowable nucleotides to insert in the sequence
// satisfying GC and run-length constraints
std::vector<std::uint8_t> dnacallowed(
    GF4reg prev,
    const DNA_params p
    ) {
    using VU = std::vector<std::uint8_t>;
    
    if (p.DNA_window == 0) {
        return VU{0,1,2,3}; // no constraint
    }
    // used to set oldest to "A", assuming DNAWINDOW is less
    // than 32 (otherwise it'll get shifted away entirely)
    const Ullong dnaoldmask((1ULL << 2 * (p.DNA_window - 1)) - 1);

    // check how many times the initial 2-bit pattern in x is repeated
    // by default checks up to all 32 2-bit slots but can stop early
    // direction of checking is from the low (least significant) bits up.
    auto run_len = [](Ullong x,unsigned stop = 32) -> unsigned {
        const auto val0 = x & 3;
        unsigned n_run = 1;
        while(n_run < stop) {
            x >>= 2;
            if( (x & 3) == val0 ) n_run++; 
            else break;
        }
        return n_run;
    };
    // is there a run and, if so, of what
    bool isrun = run_len(prev,p.max_run) >= p.max_run;
    const auto last = prev & 3;
    const auto gccount = countgc(prev,dnaoldmask);
    // the horrible logic tree:
    VU ans;
    if (gccount >= p.max_GC) {
        ans = {0,3}; // A, T OK
        if (isrun) {
            if (last == 0) {
                ans.erase(ans.begin()); // only T ok
            } else if (last == 3) {
                ans.pop_back(); // only A ok
            }
        }
    }
    else if (gccount <= p.min_GC) {
        ans = {1,2}; // C, G OK
        if (isrun) {
            if (last == 1) {
                ans.erase(ans.begin()); // only G ok
            } else if (last == 2) {
                ans.pop_back(); // only C ok
            }
        }
    }
    else {  // no GC constraints
        ans = {0,1,2,3};
        if (isrun) { // remove the 'runner':
            ans.erase(ans.begin() + last);
        }
    }
    return ans;
}

// misc. utility function that creates errors in a DNA string
std::string dna_errors(
    const std::string& s,
    double sub_rate,
    double del_rate,
    double ins_rate
) {
    static std::mt19937 rng{std::random_device{}()};
    static std::uniform_real_distribution<> runif(0.,1.);
    static std::uniform_int_distribution<> rint(0,3);

    size_t n = 0, nn = s.size();
    std::vector<unsigned char> ans;
    // reserve a sufficient capacity for most cases:
    const size_t size_guess = nn * (1 - del_rate) * (1 + ins_rate) * 1.33;
    ans.reserve(size_guess);
    while (n < nn) {
        if (runif(rng) < ins_rate) { // insertion
            ans.push_back(rint(rng));
            continue;
        }
        if (runif(rng) < del_rate) { // deletion
            ++n;
            continue;
        }
        if (runif(rng) < sub_rate) { //substitution or errorfree
            ans.push_back(rint(rng));
        } else {
            ans.push_back(s[n++]);
        }
    }
    return string_from_charvec(ans);
}

// returns penalty of match (large is bad)
// rather than converting to ACGT strings, consider using the 0,1,2,3 representation
double primerscore(
    const std::string& a, const std::string& b,
    double mispen = 1., double gappen = 1., double skewpen = 1.
    ) {
    const size_t na{a.size()}, nb{b.size()};
    simple_matrix<double> cost(na + 1, nb + 1);
    cost(0,0) = 0.;
    for (size_t i = 1; i <= na; i++) cost(i,0) = cost(i - 1,0) + skewpen;
    for (size_t i = 1; i <= nb; i++) cost(0,i) = cost(0,i - 1) + skewpen;
    for (size_t i = 1; i <= na; i++) {
        for (size_t j = 1; j <= nb; j++) {
            double dn = cost(i - 1,j) + ((j == nb) ? skewpen : gappen);
            double rt = cost(i,j - 1) + ((i == na) ? skewpen : gappen);
            double dg = cost(i - 1,j - 1) + ((a[i - 1] == b[j - 1]) ? -1. : mispen);
            cost(i,j) = std::min(std::min(dn, rt), dg);
        }
    }
    return cost(na,nb);
}


// reverse complement codeword if that makes left_primer agree better
// note the inconsistent types here: a string (ACGT) and a 0123 encoding...fix
GF4word make_sense(const std::string& left_primer,const GF4word& codeword) {
    auto rcodeword = revcomp(codeword);
    std::string scode{string_from_charvec(codeword)}, srcode{string_from_charvec(rcodeword)};
    const double lscore = primerscore(left_primer, scode);
    const double rscore = primerscore(left_primer, srcode);
    return (rscore <= lscore) ? rcodeword : codeword;
}


// method definitions for encoder class (defined in header)

void encoder::set_coderate(unsigned pattern) {
    if(pattern < 1 || pattern > 6) {
        throw("encoder: invalid pattern number outside 1 - 6 range.");
    }
    reward = rewards[pattern - 1];
    Pat = pattarr(pattern,leftprimer.size());
}

void encoder::set_primers(const std::string& leftp,const std::string& rightp) {
    const size_t nL{leftp.size()};
    leftprimer = charvec_from_dna(leftp);
    rightprimer = charvec_from_dna(rightp);
    primersalt.reserve(nL);
    // give primersalt a default value
    std::fill(std::begin(primersalt),std::begin(primersalt)+nL,0);
    constexpr unsigned num_try{100};
    for (size_t k = 0; k < nL; k++) {
        auto match = leftprimer[k];
        for (unsigned i = 0; i < num_try; i++) { // try up to 100 times (what??)
            auto regout = D.digest(0ULL, k, i, 4);
            if (regout == match) {
                primersalt[k] = i;
                break;
            }
        }
    }
    // since the pattern 'array' has a dependence on leftprimer, need
    // to reset that, too:
    set_coderate(_pattern);
}

size_t encoder::vbitlen(size_t num_msg_bits) const noexcept {
    return Pat.vbitsum(num_msg_bits);
}

VecUchar encoder::packvbits(const VecUchar& vbits, size_t msg_bits = 0) const {
    
    const auto ksize = vbits.size();
    // no more than the specified number of bits if it is set
    size_t nn = msg_bits > 0 ? msg_bits : Pat[ksize];
    const auto n_bytes = (nn + 7) / 8; // number of bytes
    VecUchar ans(n_bytes,0);
    unsigned int i = 0, j = 0;
    for (size_t k = 0; k < ksize; k++) {
        for (auto k1 = Pat[k] - 1; k1 >= 0; k1--) {
            Uchar bit = (vbits[k] >> k1) & 1;
            ans[i] = ans[i] | (bit << (7-j++));
            if (j == 8) {
                j = 0;
                if (++i == n_bytes) break; 
            }
        }
        if (i == n_bytes) break;
    }
    return ans;
}

VecUchar encoder::unpackvbits(const std::string& message, std::uint32_t len=0) const {
    const size_t n = message.size();
    size_t nmb = 8*n; // # of bits in the message
    auto ksize = vbitlen(nmb);
    // guard against unsigned underflow
    size_t size2 = len > rightprimer.size() ? len - rightprimer.size() : 0;
    // aim for codetext of length len if possible
    ksize = std::max(ksize,size2);
    VecUchar ans(ksize,0);
    size_t i = 0, j = 0;
    for (size_t k = 0; k < ksize; k++) {
        for (size_t k1 = 0; k1 < Pat[k]; k1++) {
            Uchar bit = i < n ? (message[i] >> (7 - j++)) & 1 : 0;
            if (j == 8) { j = 0; ++i; }
            ans[k] = (ans[k] << 1) | bit;
        }
    }
    return ans;
}

size_t encoder::strand_len(size_t bytes) const noexcept {
    return vbitlen(8 * bytes) + rightprimer.size();
}

// encode a message!
GF4word encoder::encode(
    std::string& message,
    size_t len,
    DNA_params p
) const {
    static const Ullong prevmask{ (1ULL << D.nprev()) - 1 };
    static const Ullong saltmask{ (1ULL << D.hsalt()) - 1 };
    const Ullong dnawinmask { (1ULL << 2 * p.DNA_window) - 1 };
    
    GF4word vbits = unpackvbits(message, len);
    const auto nm{vbits.size()}; // number of variable bits encoded
    const size_t LPRIMER{leftprimer.size()};
    GF4word codetext(nm + rightprimer.size());
    Ullong prevbits = 0, salt = 0, newsalt = 0; 
    GF4reg prevcode = acgtacgt; // initialize with no runs and balanced cg
    const size_t nsp = NSP();
    for (size_t k = 0; k < nm; k++) { // on decoding, k is called seq
        auto messagebit = vbits[k];
        if (k < LPRIMER) {
            salt = primersalt[k];
        }
        else if (k < nsp) {
            salt = 0;
            newsalt = ((newsalt << 1) & saltmask) ^ messagebit;
        }
        else if (k == nsp) {
            salt = newsalt; // time to update the salt
        }
        if(k < LPRIMER) {
            auto regout = D.digest(prevbits, k, salt, 4);
            codetext[k] = (regout + messagebit) % 4;
        } else {
            const auto dnac_ok = dnacallowed(prevcode,p);
            int mod = dnac_ok.size();
            auto regout = D.digest(prevbits, k, salt, mod);
            regout = (regout + messagebit) % mod;
            codetext[k] = dnac_ok[regout];
        }
        auto nbits = Pat[k];
        prevbits = ((prevbits << nbits) & prevmask) | messagebit; // variable number
        prevcode = ((prevcode << 2) | codetext[k]) & dnawinmask;  // shift code by 1
    }
    std::copy(std::cbegin(rightprimer),std::cend(rightprimer),std::begin(codetext) + nm);
    return codetext;
}

int encoder::search_heap(
    hypovec& hypostack,
    std::vector<Uchar>& text,
    unsigned msg_bits,
    unsigned hlimit) const {
    
    // given the heap, keep processing it until offset limit, hypothesis limit, or an error is reached
    const auto seqmax = vbitlen(msg_bits);
    const size_t limit{text.size()};
    if(limit < 1) return 0;
    HeapScheduler<double,size_t> heap(limit);
    int qq, qqmax = -1, ofmax = -1;
    constexpr std::array<int,3> skews{0,-1,1};
    while (true) {
        const auto [currscore, val2] = heap.pop();
        qq = val2;
        const auto hp = hypostack[qq];
        auto seq = hp.seq;
        
        if (hp.offset > ofmax) { // keep track of farthest gotten to
            ofmax = hp.offset;
            qqmax = qq;
        }
        if (currscore == heap.default_value()) break; // heap is empty (TODO: change this)
        if (static_cast<unsigned>(hp.offset) + 1 >= limit) break; // errcode 0 (nominal success)
        if (msg_bits > 0 && static_cast<size_t>(seq) + 1 >= seqmax) break; // ditto when no. of message bits specified
        // TODO: fix error codes
        if (hypostack.size() > hlimit) { 
            return qqmax;
        }

        const int nguess = 1 << Pat[seq + 1]; // 1, 2, or 4
        for(auto skew : skews) {
            for (Uchar mbit = 0; mbit < nguess; mbit++) {
                Hypothesis<8,24> h;
                if(h.init_from_predecessor(text, hp, qq, mbit, skew)) {
                    hypostack.push_back(h);
                    heap.push(h.score);
                }
            }
        }
    }
    return qq; // final position
}

std::string encoder::decode(GF4word data,unsigned msg_bits,unsigned hlimit) const {
    hypovec hstack;
    int final_pos = search_heap(hstack,data, msg_bits, hlimit);
    // update global params (this need not be done here?)
    // TODO - return a struct with the message and these extra fields
    // auto finalscore = hstack[final_pos].score;
    // auto finaloffset = hstack[final_pos].offset;
    // auto finalseq = hstack[final_pos].seq;
    std::string ans;
    // populate the answer
    size_t q = final_pos;
    ans.push_back(hstack[q].messagebit);
    while (hstack[q].pred_ix > 0) {
        ans.push_back(hstack[q].messagebit);
        q = hstack[q].pred_ix;
    }
    std::reverse(std::begin(ans),std::end(ans));
    return ans;
}
