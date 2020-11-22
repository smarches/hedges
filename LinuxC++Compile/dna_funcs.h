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

// reverse complement for GF4word (aka vector of bytes)
GF4word revcomp(const GF4word& arr) {
	GF4word ans = arr;
	static constexpr Uchar TGCA[] = { 3,2,1,0 };
    std::reverse(std::begin(ans),std::end(ans));
    std::for_each(std::begin(ans),std::end(ans),[](auto& e){
        e = e > 3 ? e : TGCA[e];
    });
    return ans;
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

std::string dna_errors(const std::string& s,double sub_rate,double del_rate,double ins_rate) {
    static std::mt19937 rng{std::random_device{}};
    static std::uniform_real_distrubition runif(0.,1.);
    static std::uniform_int_distribution rint(0,3);

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
			// ans[k++] = (codetext[n++] + (rint(rng) % 3) + 1) % 4;
            ans.push_back(rint(rng));
		} else {
			ans.push_back(s[n++]);
		}
	}
    return ans;
}

double primerscore(const std::string& a,const std::string& b) {
	// returns penalty of match (large is bad)
    // TODO: sync with global vars??
	static constexpr char ACGT[] = "ACGT";
	const double mispen = 1., gappen = 1., skewpen = 1.;
	const size_t ia{ain.size()}, ib{b.size()};
    // TODO: use Eigen or Armadillo?
	MatDoub cost(ia + 1, ib + 1);
	cost[0][0] = 0.;
	for (size_t i = 1; i <= ia; i++) cost[i][0] = cost[i - 1][0] + skewpen;
	for (size_t i = 1; i <= ib; i++) cost[0][i] = cost[0][i - 1] + skewpen;
	for (size_t i = 1; i <= ia; i++) {
        for (size_t j = 1; j <= ib; j++) {
		    double dn = cost[i - 1][j] + ((j == ib) ? skewpen : gappen);
		    double rt = cost[i][j - 1] + ((i == ia) ? skewpen : gappen);
		    double dg = cost[i - 1][j - 1] + ((ain[i - 1] == ACGT[bin[j - 1]]) ? -1. : mispen);
		    cost[i][j] = std::min(std::min(dn, rt), dg);
        }
	}
	return cost[ia][ib];
}

GF4word make_sense(const std::string& left_primer,const GF4word& codeword) {
	// reverse complement codeword if that makes left_primer agree better
	GF4word rcodeword = revcomp(codeword);
	const double lscore = primerscore(left_primer, codeword);
	const double rscore = primerscore(left_primer, rcodeword);
	return (rscore <= lscore) ? rcodeword : codeword;
}

std::uint8_t char2i(char c) {
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

// class responsible for performing encoding
class encoder {
    GF4word leftprimer, rightprimer;
    std::vector<Ullong> primersalt;
    std::vector<std::uint8_t> pattarr;
    size_t maxL;
    double reward;
    constexpr double rewards[] = {-0.035,-0.082,-0.127,-0.229,-0.265,-0.324};
    
public:
    encoder() = delete;
    encoder(
        std::string& lprimer,std::string& rprimer,
        size_t maxLen = 2500,int pattern = 0
    ) : maxL(maxLen) {
        set_primers(lprimer,rprimer);
        set_coderate(pattern);
    }

    void set_coderate(int pattern) {
        reward = rewards[pattern];
        std::vector<std::uint8_t> pvec;
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
        pattarr.assign(maxL + 2, 1);
        const auto n_pattern{pattern.size()}, LPRIMER{leftprimer.size()};
        std::fill(std::begin(pattarr),std::begin(pattarr) + LPRIMER,0);
        for (unsigned int i = LPRIMER; i < maxL; i++) {
            pattarr[i] = pattern[i % n_pattern];s
        }
    }
    
    void set_primers(const std::string& leftp,const std::string& rightp) {
        const auto nL{leftp.size()}, nR{rightp.size()};
        leftprimer.reserve(nL);
        rightprimer.reserve(nR);
        primersalt.assign(nL,0);
        for(const auto e: leftp)  leftprimer.push_back(char2i(e));
        for(const auto e: rightp) rightprimer.push_back(char2i(e));
        constexpr int num_try{100};
        for (size_t k = 0; k < nL; k++) {
            auto match = leftprimer[k];
            for (size_t i = 0; i < num_try; i++) { // try up to 100 times (what??)
                auto regout = digest(0ULL, k, i, 4);
                if (regout == match) {
                    primersalt[k] = i;
                    break;
                }
            }
	    }
    }
};