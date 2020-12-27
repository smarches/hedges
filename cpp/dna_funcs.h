#include <random>
#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <algorithm>
#include <type_traits>

using Ullong = unsigned long long;
using Uchar = unsigned char;
using GF4char = unsigned char; // semantically ACGT
// may want to use string if likely that small-string optimizations or etc. would be helpful
using GF4word = std::vector<Uchar>;  // semantically string of ACGT
using GF4reg = unsigned long long; // semantically a compressed GF4word
using Mbit = unsigned char; // semantically VARIABLE NUMBER of plaintext bits
using VecInt = std::vector<int>;
using VecMbit = std::vector<Mbit>;  // message bits unpacked to variable
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
// this is basically a hash function so it could be replaced with any
// suitable hashing functionality
template<unsigned SEQBITS,unsigned HSALT,unsigned NPREV>
class digester {
    static constexpr Ullong seqnomask{(1LLU << SEQBITS) - 1};
    Ullong int64(Ullong u) const noexcept {
		Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
		v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
		v *= 4768777513237032717LL;
		v ^= v << 20; v ^= v >> 41; v ^= v << 5;
		return v;
	}
    public:
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

int bytepopcount(Uchar byte) noexcept;

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

class encoder {
    
    GF4word leftprimer, rightprimer;
    std::vector<Ullong> primersalt;
    size_t NSP;
    double reward;
    static constexpr double rewards[] = {-0.035,-0.082,-0.127,-0.229,-0.265,-0.324};
    static constexpr Ullong acgtacgt{0x1b1b1b1b1b1b1b1bllu}; // "ACGTACGTACGTACGT"
    digester<10,24,8> D;
    pattarr Pat;
public:
    encoder() = delete;
    encoder(
        std::string& lprimer,std::string& rprimer,
        unsigned pattern = 1,size_t nsalt = 24
    ) {
        set_primers(lprimer,rprimer);
        set_coderate(pattern);
        NSP = vbitlen(nsalt) + leftprimer.size();
    }

    void set_coderate(unsigned pattern);
    void set_primers(const std::string& leftp,const std::string& rightp);
    size_t vbitlen(size_t num_msg_bits) const noexcept;
    VecUchar packvbits(const VecMbit& vbits, size_t msg_bits) const;
    VecMbit unpackvbits(const std::string& message, std::uint32_t len) const;
    GF4word encode(
        std::string& message,
        size_t len,
        const unsigned DNAWINDOW = 12,
        const unsigned MAXGC = 8,
        const unsigned MINGC = 4,
        const unsigned MAXRUN = 4
    ) const;

};

