#include <random>
#include <cstdint>
// Copyright (C) 1992-2007 Numerical Recipes Software
// Permission is granted to use this file in association with the HEDGES DNA Error Correcting Code only
// For additional permissions or licenses, see http://numerical.recipes

static constexpr double magic_float{2.32830643653869629E-10};
static constexpr double magic_float2{5.42101086242752217E-20};
static constexpr std::uint64_t magic_ull{4101842887655102017LL};

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(magic_ull), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}

	Ran() : v(magic_ull), w(1) {
		Ullong j = std::random_device{}();
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}

	Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; 
        v ^= v << 31; 
        v ^= v >> 8;
		w = 4294957665U * (w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); 
        x ^= x >> 35; 
        x ^= x << 4;
		return (x + v) ^ w;
	}
	double doub() { return magic_float2 * int64(); }
	std::uint32_t int32() { return static_cast<std::uint32_t>(int64()); }

	// Ullong timeseed() { // should change about every millisecond (only)
	// 	time_t x;
	// 	return 1000 * Ullong(time(&x)) + Ullong(clock());
	// }

};

struct Ranq1 {
	Ullong v;
	Ranq1(Ullong j) : v(magic_ull) {
		v ^= j;
		v = int64();
	}
	Ullong int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	double doub() { return magic_float2 * int64(); }
	std::uint32_t int32() { return static_cast<std::uint32_t>(int64()); }
};

struct Ranq2 {
	Ullong v,w;
	Ranq2(Ullong j) : v(magic_ull), w(1) {
		v ^= j;
		w = int64();
		v = int64();
	}
	Ullong int64() {
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		return v ^ w;
	}
	double doub() { return magic_float2 * int64(); }
	std::uint32_t int32() { return static_cast<std::uint32_t>(int64()); }
};

struct Ranhash {
	Ullong int64(Ullong u) {
		Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
		v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
		v *= 4768777513237032717LL;
		v ^= v << 20; v ^= v >> 41; v ^= v << 5;
		return  v;
	}
	std::uint32_t int32(Ullong u)
		{ return (std::uint32_t)(int64(u) & 0xffffffff) ; }
	double doub(Ullong u)
		{ return magic_float2 * int64(u); }
};

struct Ranbyte {
	int s[256],i,j,ss;
	Ranbyte(int u) {
		std::uint32_t v = 2244614371U ^ u;
		for (i=0; i<256; i++) {s[i] = i;}
		for (j=0, i=0; i<256; i++) {
			ss = s[i];
			j = (j + ss + (v >> 24)) & 0xff;
			s[i] = s[j]; s[j] = ss;
			v = (v << 24) | (v >> 8);
		}
		i = j = 0;
        // calling just for the side effect?
		for (int k=0; k<256; k++) int8();
	}
	std::uint8_t int8() {
		i = (i+1) & 0xff;
		ss = s[i];
		j = (j+ss) & 0xff;
		s[i] = s[j]; s[j] = ss;
		return static_cast<std::uint8_t>(s[(s[i]+s[j]) & 0xff]);
	}
	std::uint32_t int32() {
		std::uint32_t v = 0;
		for (int k=0; k<4; k++) {
			i = (i+1) & 0xff;
			ss = s[i];
			j = (j+ss) & 0xff;
			s[i] = s[j]; s[j] = ss;
			v = (v << 8) | s[(s[i]+s[j]) & 0xff];
		}
		return v;
	}
	double doub() {
		return magic_float * ( int32() +
			   magic_float * int32() );
	}
};

struct Ranfib {
	double dtab[55], dd;
	int inext, inextp;
	Ranfib(Ullong j) : inext(0), inextp(31) {
		Ranq1 init(j);
		for (int k=0; k<55; k++) dtab[k] = init.doub();
	}
	double doub() {
		if (++inext == 55) inext = 0;
		if (++inextp == 55) inextp = 0;
		dd = dtab[inext] - dtab[inextp];
		if (dd < 0) dd += 1.0;
		return (dtab[inext] = dd);
	}
	unsigned long int32() { 
        return static_cast<unsigned long>(doub() * 4294967295.0);
    }
};

struct Ranlim32 {
	std::uint32_t u,v,w1,w2;
	Ranlim32(std::uint32_t j) : v(2244614371U), w1(521288629U), w2(362436069U) {
		u = j ^ v; int32();
		v = u; int32();
	}
	std::uint32_t int32() {
		u = u * 2891336453U + 1640531513U;
		v ^= v >> 13; 
        v ^= v << 17; 
        v ^= v >> 5;
		w1 = 33378 * (w1 & 0xffff) + (w1 >> 16);
		w2 = 57225 * (w2 & 0xffff) + (w2 >> 16);
		std::uint32_t x = u ^ (u << 9); 
        x ^= x >> 17; 
        x ^= x << 6;
		std::uint32_t y = w1 ^ (w1 << 17); 
        y ^= y >> 15; 
        y ^= y << 5;
		return (x + v) ^ (y + w2);
	}
	double doub() { return magic_float * int32(); }
    // so is 'doub' a 'false doub'??
	double truedoub() {
		return magic_float * ( int32() + magic_float * int32() );
	}
};

