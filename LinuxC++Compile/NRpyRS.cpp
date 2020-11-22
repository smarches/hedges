#include "nr3python.h"
#include <vector>
#include <random>
#include <algorithm>
#include "reed_solomon_schifra.h"

SchifraCode<255, 32> rs;

static std::mt19937 rng{std::random_device{}()};
static std::uniform_real_distribution runif(0.,1.);

using VecUchar = std::vector<unsigned char>;

static PyObject* rsencode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	if (args.size() != 1) {
		NRpyException("rsencode takes 1 argument only");
		return NRpyObject(0); // formerly NULL
	}
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) {
		NRpyException("rsencode requires array with dtype=uint8 \n");
		return NRpyObject(0);
	}
	VecUchar message(args[0]);
	if (message.size() != 255) {
		NRpyException("rsencode requires input array of size exactly 255");
		return NRpyObject(0);
	}
	VecUchar codetext = rs.encode(message);
	return NRpyObject(codetext);
}

static PyObject* rsdecode(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	VecUchar received;
	std::vector<int> locations;
	if (PyArray_TYPE(args[0]) != PyArray_UBYTE) {
		NRpyException("rsdecode requires array with dtype=uint8 \n");
		return NRpyObject(0);
	}
	if (args.size() == 1) {
		received = VecUchar(args[0]);
		locations.resize(0);
	} else if (args.size() == 2) {
		received = VecUchar(args[0]);
		locations = std::vector<int>(args[1]);
	} else {
		NRpyException("rsencode takes 1 or 2 arguments");
		return NRpyObject(0);
	}
	int errs_detected, errs_corrected, err_code;
	bool recoverable;
	VecUchar decoded = rs.decode(received, locations, errs_detected, errs_corrected, err_code, recoverable);
	return NRpyTuple(
		NRpyObject(decoded),
		NRpyObject(errs_detected),
		NRpyObject(errs_corrected),
		NRpyObject(err_code),
		NRpyObject(Int(recoverable)),
		NULL
	);
}

static PyObject* makeerasures(PyObject *self, PyObject *pyargs) {
    
	NRpyArgs args(pyargs);
	const VecUchar incodeword(args[0]);
	size_t nerase{NRpyInt(args[1])};
    const auto L{incodeword.size()};
	VecUchar codeword;
    codeword.reserve(L - nerase);
	std::vector<int> dlocs;
    dlocs.reserve(nerase);
	size_t j{0};
    for (size_t i=0;i<L;i++) {
        double p = static_cast<double>(nerror - j) / (L - i);
        if (runif(rng) < p) {
            dlocs.push_back(i); // a deletion site
			j++; 
		} else {
			codeword.push_back(incodeword[i]); // wraparound; anything to change it is ok
        }
    }
	return NRpyTuple(
		NRpyObject(codeword),
		NRpyObject(dlocs),
		NULL
	);
}

static PyObject* makeerrors(PyObject *self, PyObject *pyargs) {
	NRpyArgs args(pyargs);
	VecUchar incodeword(args[0]);
	VecUchar codeword(incodeword); // so won't alter incodeword
	int nerror(NRpyInt(args[1])), merror{0};
	const auto nn{codeword.size()};
	for (size_t i = 0; i < nn; i++) {
		double p = static_cast<double>(nerror - merror) / (nn - i);
		if (runif(rng) < p) {
			codeword[i] += 1; // wraparound; anything to change it is ok
			++merror;
		}
	}
	return NRpyObject(codeword);
}


// standard boilerplate
static PyMethodDef NRpyRS_methods[] = {
	{ "rsencode", rsencode, METH_VARARGS,
	"codetext = rsencode(uint8_array_length_255)" },
	{ "rsdecode", rsdecode, METH_VARARGS,
	"(decoded, errs_detected, errs_corrected, err_code, recoverable) = rsdecode(uint8_array_length_255[, std::vector<int> erasure_locations])" },
	{ "makeerasures", makeerasures, METH_VARARGS,
	"(newcodetext,locations) = makeerasures(codetext,nerasures)" },
	{ "makeerrors", makeerrors, METH_VARARGS,
	"newcodetext = makeerrors(codetext,nerrors)" },
	{ NULL, NULL, 0, NULL }
};
PyMODINIT_FUNC initNRpyRS(void) {
	import_array();
	Py_InitModule("NRpyRS", NRpyRS_methods);
}

