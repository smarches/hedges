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

// not being used, but might someday!
int bytepopcount(Uchar byte) {
  static constexpr Uchar NIBBLE_LOOKUP[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
  return NIBBLE_LOOKUP[byte & 0x0F] + NIBBLE_LOOKUP[byte >> 4];
}

// in-place reverse complement for GF4word (aka vector of bytes)
void revcomp_C(GF4word& arr) {
	const auto len = arr.size();
	static constexpr Uchar TGCA[] = { 3,2,1,0 };
    std::reverse(std::begin(arr),std::end(arr));
    std::for_each(std::begin(arr),std::end(arr),[](auto& e){
        e = e > 3 ? e : TGCA[e];
    });
}

// length of return value is used as the 'mod' value in digest()
// values in return val are the allowable nucleotides to insert in the sequence
std::vector<std::uint8_t> dnacallowed(
    GF4reg prev,
    const unsigned DNAWINDOW,
    const unsigned MAXGC,
    const unsigned MINGC,
    const unsigned MAXRUN 
    ) {
    using VU = std::vector<std::uint8_t>;
	// returns the number of allowed ACGTs and puts them in dnac_ok
	if (DNAWINDOW == 0) {
        VU res{0,1,2,3};
        return res;
    }
    const Ullong dnaoldmask((1ULL << 2 * (DNAWINDOW-1)) - 1); // used to set oldest to "A"
    
    auto countgc = [dnaoldmask](const Ullong prev,const Ullong mask) -> Ullong {
        // get GCcount
        Ullong reg = prev & dnaoldmask;
        reg = (reg ^ (reg >> 1)) & 0x5555555555555555ull; // makes ones for GC, zeros for AT
        // popcount inline:
        reg -= ((reg >> 1) & 0x5555555555555555ull);
        reg = (reg & 0x3333333333333333ull) + (reg >> 2 & 0x3333333333333333ull);
        return ((reg + (reg >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
    };

	const int gccount = countgc(prev,dnaoldmask);

	// is there a run and, if so, of what
	Ullong reg = prev >> 2;
	bool isrun = false;
    int nrun = 1, last = prev & 3;
	while ((reg & 3) == last) {
		++nrun;
		if (nrun >= MAXRUN) { isrun = true; break; }
		reg >>= 2;
	}
	// the horrible logic tree:
    VU ans;
    ans.reserve(4);
	if (gccount >= MAXGC) {
		ans.push_back(0); // A is ok
		ans.push_back(3); // T is ok
		if (isrun) {
			if (last == 0) {
                ans.erase(ans.begin()); // only T ok
			} else if (last == 3) {
                ans.pop_back(); // only A ok
			}
		}
	}
	else if (gccount <= MINGC) {
		ans.push_back(1); // C is ok
		ans.push_back(2); // G is ok
		if (isrun) {
			if (last == 1) {
				ans.erase(ans.begin()); // only G ok
			} else if (last == 2) {
				ans.pop_back(); // only C ok
			}
		}
	}
	else {  // no GC constraints
		ans = VU({0,1,2,3});
		if (isrun) { // remove the 'runner':
            ans.erase(ans.begin() + last);
		}
	}
	return ans;
}
