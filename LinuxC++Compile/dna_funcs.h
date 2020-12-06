#include <random>
#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <algorithm>

using Ullong = unsigned long long;
using Uchar = unsigned char;
using GF4char = unsigned char; // semantically ACGT
// may want to use string if likely that small-string optimizations or etc. would be helpful
using GF4word = std::vector<Uchar>;  // semantically string of ACGT
using GF4reg = unsigned long long; // semantically a compressed GF4word
using Mbit = unsigned char; // semantically VARIABLE NUMBER of plaintext bits
using VecInt = std::vector<int>;
using VecMbit = std::vector<Mbit>;  // message bits unpacked to variable

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
    Ullong int64(Ullong u) {
		Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
		v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
		v *= 4768777513237032717LL;
		v ^= v << 20; v ^= v >> 41; v ^= v << 5;
		return v;
	}
    public:    
    std::int32_t digest(Ullong bits, std::int32_t seq, Ullong salt, std::int32_t mod) {
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

std::string dna_errors(const std::string& s,double sub_rate,double del_rate,double ins_rate);

double primerscore(const std::string& a,const std::string& b,double mispen,double gappen, double skewpen);

GF4word make_sense(const std::string& left_primer,const GF4word& codeword);

class encoder;
