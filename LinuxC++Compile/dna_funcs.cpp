#include <type_traits>
#include "dna_funcs.h"

// not being used, but might someday!
int bytepopcount(Uchar byte) noexcept {
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
auto countgc(const Ullong prev,const Ullong mask) noexcept {
        Ullong reg = prev & mask;
        // makes ones for GC, zeros for AT
        reg = (reg ^ (reg >> 1)) & 0x5555555555555555ull;
        reg -= ((reg >> 1) & 0x5555555555555555ull);
        reg = (reg & 0x3333333333333333ull) + (reg >> 2 & 0x3333333333333333ull);
        return ((reg + (reg >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

// length of return value is used as the 'mod' value in digest()
// values in return val are the allowable nucleotides to insert in the sequence
// satisfying GC and run-length constraints
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
    // used to set oldest to "A"
    const Ullong dnaoldmask((1ULL << 2 * (DNAWINDOW-1)) - 1);
	const auto gccount = countgc(prev,dnaoldmask);

	// is there a run and, if so, of what
	Ullong reg = prev >> 2;
	bool isrun = false;
    unsigned nrun = 1, last = prev & 3;
	while ((reg & 3) == last) {
		++nrun;
		if (nrun >= MAXRUN) { isrun = true; break; }
		reg >>= 2;
	}
	// the horrible logic tree:
    VU ans;
	if (gccount >= MAXGC) {
        ans = {0,3}; // A, T OK
		if (isrun) {
			if (last == 0) {
                ans.erase(ans.begin()); // only T ok
			} else if (last == 3) {
                ans.pop_back(); // only A ok
			}
		}
	}
	else if (gccount <= MINGC) {
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
    static std::uniform_real_distribution runif(0.,1.);
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
// note the confusing types here: a string (ACGT) and a 0123 encoding...fix
GF4word make_sense(const std::string& left_primer,const GF4word& codeword) {
	auto rcodeword = revcomp(codeword);
    std::string scode{string_from_charvec(codeword)}, srcode{string_from_charvec(rcodeword)};
	const double lscore = primerscore(left_primer, scode);
	const double rscore = primerscore(left_primer, srcode);
	return (rscore <= lscore) ? rcodeword : codeword;
}

// class responsible for performing encoding
class encoder {
    GF4word leftprimer, rightprimer;
    std::vector<Ullong> primersalt;
    std::vector<std::uint8_t> pattarr;
    size_t maxL, NSP;
    double reward;
    static constexpr double rewards[] = {-0.035,-0.082,-0.127,-0.229,-0.265,-0.324};
    static constexpr Ullong acgtacgt{0x1b1b1b1b1b1b1b1bllu}; // "ACGTACGTACGTACGT"
    digester<10,24,8> D;
    
public:
    encoder() = delete;
    encoder(
        std::string& lprimer,std::string& rprimer,
        size_t maxLen = 2500,size_t nsalt = 24,int pattern = 0
    ) : maxL(maxLen) {
        set_primers(lprimer,rprimer);
        set_coderate(pattern);
        NSP = vbitlen(nsalt) + leftprimer.size();
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
        const auto n_pattern{pvec.size()}, LPRIMER{leftprimer.size()};
        std::fill(std::begin(pattarr),std::begin(pattarr) + LPRIMER,0);
        for (unsigned int i = LPRIMER; i < maxL; i++) {
            pattarr[i] = pvec[i % n_pattern];
        }
    }
    
    void set_primers(const std::string& leftp,const std::string& rightp) {
        const auto nL{leftp.size()};
        leftprimer = charvec_from_dna(leftp);
        rightprimer = charvec_from_dna(rightp);
        constexpr int num_try{100};
        for (size_t k = 0; k < nL; k++) {
            auto match = leftprimer[k];
            for (size_t i = 0; i < num_try; i++) { // try up to 100 times (what??)
                auto regout = D.digest(0ULL, k, i, 4);
                if (regout == match) {
                    primersalt[k] = i;
                    break;
                }
            }
	    }
    }

    size_t vbitlen(size_t num_msg_bits) const {
        size_t ksize{0}, nn{0};
        while( ksize++ ) {  // how many Mbits do we need?
            if (nn >= num_msg_bits) break;
            if (ksize >= pattarr.size()) throw "vbitlen: MAXSEQ too small";
            nn += pattarr[ksize];
	    }
	    return ksize;
    }

    VecMbit unpackvbits(const std::string& message, std::uint32_t len=0) {
        const size_t n{message.size()};
        size_t nmb{8*n}; // # of bits in the message
        auto ksize = vbitlen(nmb);
        // guard against unsigned underflow
        size_t size2 = len > rightprimer.size() ? len - rightprimer.size() : 0;
        // aim for codetext of length len if possible
        ksize = std::max(ksize,size2);
        VecMbit ans(ksize,0);
        size_t i = 0, j = 0;
        for (size_t k = 0; k < ksize; k++) {
            for (int k1 = 0; k1 < pattarr[k]; k1++) {
                Uchar bit = i < n ? (message[i] >> (7 - j++)) & 1 : 0;
                if (j == 8) { j = 0; ++i; }
                ans[k] = (ans[k] << 1) | bit;
            }
        }
        return ans;
    }
    
    // encode a message!
    GF4word encode(
        std::string& message,
        size_t len,
        const unsigned DNAWINDOW = 12,
        const unsigned MAXGC = 8,
        const unsigned MINGC = 4,
        const unsigned MAXRUN = 4
    ) {
        static const Ullong prevmask{ (1ULL << D.nprev()) - 1 };
        static const Ullong saltmask{ (1ULL << D.hsalt()) - 1 };
        const Ullong dnawinmask { (1ULL << 2 * DNAWINDOW) - 1 };
        GF4word vbits = unpackvbits(message, len);
	    const auto nm{vbits.size()}; // number of variable bits encoded
	    if (nm > pattarr.size()) throw("encode: MAXSEQ too small");

        const size_t LPRIMER{leftprimer.size()};
	    GF4word codetext(nm + rightprimer.size());
	
	    Ullong prevbits = 0, salt = 0, newsalt = 0; 
	    GF4reg prevcode = acgtacgt; // initialize with no runs and balanced cg
        for (size_t k = 0; k < nm; k++) { // on decoding, k is called seq
            Mbit messagebit = vbits[k];
            if (k < LPRIMER) {
                salt = primersalt[k];
            }
            else if (k < NSP) {
                salt = 0;
                newsalt = ((newsalt << 1) & saltmask) ^ messagebit;
            }
            else if (k == NSP) {
                salt = newsalt; // time to update the salt
            }
            if(k < LPRIMER) {
                auto regout = D.digest(prevbits, k, salt, 4);
                codetext[k] = (regout + messagebit) % 4;
            } else {
                const auto dnac_ok = dnacallowed(prevcode,DNAWINDOW,MAXGC,MINGC,MAXRUN);
                int mod = dnac_ok.size();
                auto regout = D.digest(prevbits, k, salt, mod);
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
};
