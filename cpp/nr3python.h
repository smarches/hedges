/* nr3python.h */
// version 0.5 checks wordlength of passed integer arrays
// version 0.4
// This file is a version of nr3.h with hooks that
// make it easy to interface to Python
// See http://www.nr.com/nr3_python_tutorial.html
#ifndef _NR3_H_
#define _NR3_H_
#ifndef Py_PYTHON_H
#include "Python.h"
#endif
#ifndef Py_ARRAYOBJECT_H
#include "numpy/arrayobject.h"
#endif

constexpr static bool _CHECKBOUNDS_{true};
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
#define _USEPYERRORCLASS_ 1
#define _TURNONFPES_ 1

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <ctime>
// #include <fcntl.h>
#include <string>
#include <cctype>

// NaN: see https://en.cppreference.com/w/cpp/types/numeric_limits/quiet_NaN
static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
// somewhat silly method where we explicitly set the bits of the exponent:
// unsigned int proto_nan[2]={0xffffffff, 0x7fffffff};
// double NaN = *( double* )proto_nan;


// Python glue Part I starts here (Part II at end of file)

PyObject* NRpyException(const char *str, int die=1, int val=0) {
	PySys_WriteStderr("ERROR: ");
	PySys_WriteStderr(str);
	PySys_WriteStderr("\n");
	//PyErr_BadArgument();
	//PyErr_Format(PyExc_RuntimeError,str,val);
	//PyErr_Print();
	//PyErr_SetInterrupt();
	PyErr_CheckSignals(); // causes a KeyboardInterrupt, only way I know to get back to the interpreter
	//PyErr_SetInterrupt();
	//PyErr_CheckSignals(); // maybe twice works better!
	return Py_None;
}

constexpr char NRpyMainName[] = "__main__";

PyObject* NRpyGetByName(char *name, char *dict = nullptr) {
// get a PyObject from Python namespace (__main__ by default)
	if (dict == nullptr) dict = NRpyMainName;
	PyObject* pymodule = PyImport_AddModule(dict);
	PyObject* dictobj = PyModule_GetDict(pymodule);
	PyObject* bb = PyDict_GetItemString(dictobj,name);
	if (! bb) NRpyException("Failed to fetch a Python object by name.");
	return bb;
}

struct NRpyArgs {
// make arguments from Python individually addressable by [] subscript
	PyObject* pyargs;
	NRpyArgs(PyObject* pyaargs) : pyargs(pyaargs) {}
	int size() const {
        return static_cast<int>(PyTuple_Size(pyargs));
    }
	PyObject* operator[](int i) {
		if (i < PyTuple_Size(pyargs)) return PyTuple_GetItem(pyargs,i);
		// Returns a borrowed (unprotected) ref.  Args refs owned by calling routine.
		else NRpyException("Index out of range in NRpyArgs.");
		return Py_None; // unreachable?
	}
};

// explicitly construct scalars and strings from PyObjects or Python namespace

bool NRpyIsNumber(PyObject* ob) {
    return (PyInt_Check(ob) || PyFloat_Check(ob));
}

int NRpyInt(PyObject* ob) {
	if (ob == Py_None) return 0; // NaN is in fact 0 for integral types
	if (NRpyIsNumber(ob)) return static_cast<int>(PyInt_AsLong(ob));
	else NRpyException("NRpyInt argument is not a number.");
	return 0;
}

int NRpyInt(char *name, char *dict = nullptr) {
	return NRpyInt(NRpyGetByName(name,dict));
}

double NRpyDoub(PyObject* ob) {
	if (ob == Py_None) return NaN;
	else if (NRpyIsNumber(ob)) return static_cast<double>(PyFloat_AsDouble(ob)); // casts ob to double
	else NRpyException("NRpyDoub argument is not a number.");
	return 0.;
}
double NRpyDoub(char *name, char *dict = nullptr) {
	return NRpyDoub(NRpyGetByName(name,dict));
}

// type char* (string)
char* NRpyCharP(PyObject *ob) {
	if (PyString_Check(ob)) return PyString_AsString(ob);
	else NRpyException("NRpyCharP argument is not a string.");
	return NULL;
}
char* NRpyCharP(char *name, char *dict = nullptr) {
	return NRpyCharP(NRpyGetByName(name,dict));
}

// type encapsulated function pointer (note different syntax so that templating can work)
template<class T>
void NRpyCFunction(T* fptr, PyObject* ob) {
	if (! PyCapsule_CheckExact(ob)) NRpyException("NRpyCFunction arg is not a C++ function capsule.");
	fptr = static_cast<T*>(PyCapsule_GetPointer(ob,nullptr));
	return;
}

// wrapper class for Python List, implementing only simple operations
struct NRpyList {
	PyObject* p;
	int n;
	bool isnew;
	NRpyList(int nn) : p(PyList_New(nn)), n(nn), isnew(true) {
		for (int i=0;i<nn;i++) {
			Py_INCREF(Py_None); // needed?
			PyList_SetItem(p,i,Py_None);
		}
		if (! PyList_Check(p)) NRpyException("NRpyList not successfully created.");
	}
    // we don't enforce that pp is actually a *list* as opposed to any other PyObject?
	NRpyList(PyObject *pp) : p(pp), isnew(false) {
		if (p == nullptr) p = Py_None;
		n = PyList_Check(p) ? PyList_Size(p) : 0;
	}
	int size() const noexcept { return n; }
	template <class T> 
    int set(int i, T val) {
		int flag = PyList_SetItem(p, i, NRpyObject(val));
		return flag;
	}
    // Python list can be accessed via negative index
	NRpyList operator[](int i) {
		if (i >= n || i < -n) NRpyException("NRpyList subscript out of range.");
        if (i < 0) i = n - i;
		return NRpyList(PyList_GetItem(p,i));
		// Returns a borrowed (unprotected) ref, but assumes List is bound by calling routine.
	}
};

// cast list to integer value of its 1st element
// what's the purpose of this?
int NRpyInt(NRpyList& list) {
	return NRpyInt(PyList_GetItem(list.p,0)); 
}
// ToDo: also make NRpyList constructors for NRvector and NRmatrix

// wrapper class for Python Dict
struct NRpyDict {
	PyObject* p;
	bool isnew;
	NRpyDict() : p(PyDict_New()), isnew(true) {}
	NRpyDict(PyObject *pp) : p(pp), isnew(false) {
		if (! PyDict_Check(pp)) NRpyException("Argument not a dict in NRpyDict constructor.");
	}
	template <class T, class U> 
    int set(T key, U val) {
		PyObject *thekey = NRpyObject(key), *theval = NRpyObject(val);
		int flag = PyDict_SetItem(p, thekey, theval);
		Py_DECREF(thekey); // because PyDict_SetItem INCREFs both key and val
		Py_DECREF(theval);
		return flag; // returns 0 for success, -1 for failure
	}
	template <class T> 
    PyObject* get(const T key) {
		PyObject *thekey = NRpyObject(key);
		PyObject *theval = PyDict_GetItem(p, thekey); // borrowed ref
		Py_DECREF(thekey);
		if (theval) return theval; // intended use is immediate conversion so borrowed ref is ok
		else return Py_None; // ditto, won't be decremented because won't be returned to Python
	}
};

// overloaded functions to turn anything into a PyObject (vector and matrix are in Part II below)

template<class T> 
PyObject* NRpyObject(T &a) {
// default applies to all function objects or other structs
	PyObject *thing = PyCapsule_New((void*)a,NULL,NULL);
	return thing;
}

PyObject* NRpyObject(const double a) {return PyFloat_FromDouble(a);}
PyObject* NRpyObject(const int a) {return PyInt_FromLong(a);}
PyObject* NRpyObject(const unsigned long long a) {return PyInt_FromSize_t(a);}
PyObject* NRpyObject(const bool a) {return PyBool_FromLong(a);}
PyObject* NRpyObject(const char *a) {return PyString_FromString(a);} // string is copied
PyObject* NRpyObject() {return Py_BuildValue("");} // Python None
PyObject* NRpyObject(NRpyList& a) {
	if (! a.isnew) Py_INCREF(a.p); // make a new reference to return (should happen rarely)
	return a.p;
}
PyObject* NRpyObject(NRpyDict &a) {
	if (! a.isnew) Py_INCREF(a.p); // make a new reference to return (should happen rarely)
	return a.p;
}
// NRpyObjects are generally return values.  Hence they pass clear title to ownership.
// If you create an NRpyObject and don't return it, you need to Py_DECREF it later.

// send an object into Python namespace (except for scalars, will become shared ref)
template <class T>
void NRpySend(T& a, char *name, char *dict=nullptr) {
	if (dict == nullptr) dict = NRpyMainName;
	PyObject* pymodule = PyImport_AddModule(dict);
	PyObject* dictobj = PyModule_GetDict(pymodule);
	PyObject* aa = NRpyObject(a);
	int ret = PyDict_SetItemString(dictobj, name, aa);
	if (ret) NRpyException("Failed to share an NR object with Python.");
	Py_XDECREF(aa); // because dictobj now has the responsibility
}

// templated check of a PyObject's type (used in initpyvec and initpymat below)
// TODO: use some type traits or concepts to clean this up?
template <class T> inline int NRpyTypeOK(PyObject *a) {return 0;}
template <> inline int NRpyTypeOK<double>(PyObject *a) {return PyArray_ISFLOAT(a);}
template <> inline int NRpyTypeOK<int>(PyObject *a) {
	return PyArray_ISINTEGER(a) && (PyArray_TYPE(a) == NPY_INT32);
}
template <> inline int NRpyTypeOK<char>(PyObject *a) {return PyArray_ISINTEGER(a);}
template <> inline int NRpyTypeOK<unsigned char>(PyObject *a) {return PyArray_ISINTEGER(a);}


// templated return a PyObject's type (used in NRpyObject on vector and matrix args)
template <class T> inline int NRpyDataType() {return PyArray_INT;}
template <> inline int NRpyDataType<double>() {return PyArray_DOUBLE;}
template <> inline int NRpyDataType<int>() {return PyArray_INT;}
template <> inline int NRpyDataType<char>() {return PyArray_BYTE;}
template <> inline int NRpyDataType<unsigned char>() {return PyArray_UBYTE;}

// templated cast a PyObject's type (used in NRpyPyFunction for return type)
template <class T>  T NRpyCast(PyObject *a) {return (T*)nullptr; }
template <> double NRpyCast<double>(PyObject *a) {return NRpyDoub(a);}
template <> int NRpyCast<int>(PyObject *a) {return NRpyInt(a);}
template <> char* NRpyCast<char*>(PyObject *a) {return NRpyCharP(a);}

// end Python glue Part I  (see Part II at end of file)

// macro-like inline functions

// template<class T>
// T SQR(const T a) {return a*a;}

// this is NOT the usual sign function...
// needs cleanup so that comparing e.g. <int,float> does the right thing
// template<class T,class U>
// T SIGN(const T a, const U b) {
//     return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
// }

// exception handling

#ifdef _USENRERRORCLASS_
struct NRerror {
	char *message;
	char *file;
	int line;
	NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err) {
	printf("ERROR: %s\n     in file %s at line %d\n",
		err.message, err.file, err.line);
	exit(1);
}
#elif defined _USEPYERRORCLASS_
#define throw(message) NRpyException(message,1,0);
#else
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#endif

// usage example:
//
//	try {
//		somebadroutine();
//	}
//	catch(NRerror s) {NRcatch(s);}
//
// (You can of course substitute any other catch body for NRcatch(s).)


// Vector and Matrix Classes

// TODO: either jettison this or add .begin(), .end()
// and a few utils like push_back and move constructor

#ifdef _USESTDVECTOR_
#define NRvector vector
#else

template <class T>
class NRvector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	bool ownsdata; // true for normal NRvector, false if Python owns the data
	PyObject *pyident; // if I don't own my data, who does?
	NRvector();
	explicit NRvector(int n);		// Zero-based array
	NRvector(PyObject *a);			// construct from Python array
	NRvector(char *name, char *dict = nullptr); // construct from name in Python scope
	void initpyvec(PyObject *a); // helper function used by above
	NRvector(int n, const T &a);	//initialize to constant value
	NRvector(int n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor
	NRvector& operator=(const NRvector &rhs);	//assignment
	typedef T value_type; // make T available externally
	T& operator[](const int i);	//i'th element
	const T& operator[](const int i) const;
	int size() const;
	void resize(int newn, bool preserve=false); // resize 
	//void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T& a); // resize and assign a constant value
	void assign(char *name, char *dict = nullptr); // assign to a name in Python scope
	~NRvector();
};

// NRvector definitions

template <class T>
NRvector<T>::NRvector() : nn(0), v(nullptr), ownsdata(true){}

template <class T>
NRvector<T>::NRvector(int n) : nn(n), ownsdata(true), v(n>0 ? static_cast<T*>(PyMem_Malloc(n*sizeof(T))) : nullptr) {}

template <class T>
void NRvector<T>::initpyvec(PyObject *a) {
	ownsdata = false;
	pyident = a;
	if (! PyArray_CheckExact(a)) NRpyException("PyObject is not an Array in NRvector constructor.");
	if (! PyArray_ISCARRAY_RO(a)) NRpyException("Python Array must be contiguous (e.g., not strided).");
	if (! NRpyTypeOK<T>(a)) NRpyException("Python Array type does not agree with NRvector type.");
	int i, ndim = PyArray_NDIM(a);
	nn = 1;
	for (i=0;i<ndim;i++) nn *= int(PyArray_DIMS(a)[i]);
	v = (nn>0 ? (T*)PyArray_DATA(a) : nullptr);
}
template <class T> NRvector<T>::NRvector(PyObject *a) {
	initpyvec(a);
}
template <class T> NRvector<T>::NRvector(char *name, char *dict ) {
	initpyvec(NRpyGetByName(name,dict));
}

template <class T>
NRvector<T>::NRvector(int n, const T& a) : 
    nn(n), ownsdata(true), 
    v(n>0 ? static_cast<T*>(PyMem_Malloc(n*sizeof(T))) : nullptr)
{
	for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), ownsdata(true), v(n>0 ? (T*)PyMem_Malloc(n*sizeof(T)) : nullptr)
{
	for(int i=0; i<n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), ownsdata(true),
	v(nn>0 ? (T*)PyMem_Malloc(nn*sizeof(T)) : nullptr) {
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
NRvector<T>& NRvector<T>::operator=(const NRvector<T> &rhs) {
	if (this != &rhs) {
		if (nn != rhs.nn) {
			resize(rhs.nn);
			nn=rhs.nn;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

//subscripting
template <class T>
T& NRvector<T>::operator[](const int i) {
    if constexpr (_CHECKBOUNDS_) {
        if (i<0 || i>=nn) {
	        throw("NRvector subscript out of bounds");
        }
    }
	return v[i];
}

template<class T>
const T& NRvector<T>::operator[](const int i) const { //subscripting
    if constexpr(_CHECKBOUNDS_) {
        if (i<0 || i>=nn) {
	        throw("NRvector subscript out of bounds");
        }
    }
	return v[i];
}

template <class T>
int NRvector<T>::size() const noexcept {
	return nn;
}

template <class T>
void NRvector<T>::resize(int newn, bool preserve) {
	if (newn != nn) {
		if (ownsdata) {
			if (preserve) {
				const auto nmin = std::min(nn,newn);
				T *vsave = v;
				//v = newn > 0 ? new T[newn] : nullptr;
				v = newn > 0 ? (T*)PyMem_Malloc(newn*sizeof(T)) : nullptr;
				for (int i=0;i<nmin;i++) v[i] = vsave[i];
				for (int i=nmin;i<newn;i++) v[i] = T{};
				//if (vsave != nullptr) delete[] (vsave);
				if (vsave != nullptr) PyMem_Free(vsave);
				nn = newn;
			} else {
				nn = newn;
				if (v != nullptr) PyMem_Free(v);
				v = nn > 0 ? (T*) PyMem_Malloc(nn*sizeof(T)) : nullptr;
			}
		} else { // Python
			if (preserve) NRpyException("resize Python array with preserve contents not implemented");
			nn = newn;
			int dm[1];
			dm[0] = newn;
			PyArray_Dims mydims;
			mydims.ptr = (npy_intp*)dm;
			mydims.len = 1;
			PyArray_Resize((PyArrayObject *)pyident, &mydims, 0, NPY_CORDER);
			// the return value is garbage, or maybe PyNone, contrary to Numpy docs
			// I think it's a Numpy bug, but following is correct
			v = nn>0 ? (T*)PyArray_DATA(pyident) : nullptr;
		}
	}
}

template <class T>
void NRvector<T>::assign(int newn, const T& a) {
	resize(newn);
	for (int i=0;i<nn;i++) v[i] = a;
}

template <class T>
void NRvector<T>::assign(char *name, char *dict ) {
	if (! ownsdata) NRpyException("Attempt to assign Python array to another Python array.");
	if (v != nullptr) PyMem_Free(v);
	initpyvec(NRpyGetByName(name,dict));
}

template <class T>
NRvector<T>::~NRvector()
{
	if (v != nullptr && ownsdata) {
		PyMem_Free(v);
	}
}

// end of NRvector definitions

#endif //ifdef _USESTDVECTOR_

template <class T>
class NRmatrix {
private:
	int nn;
	int mm;
	T **v;
public:
	bool ownsdata; // 1 for normal NRmatrix, 0 if Python owns the data
	PyObject *pyident;
	NRmatrix();
	NRmatrix(int n, int m);			// Zero-based array
	NRmatrix(PyObject *a); // construct from Python array
	NRmatrix(char *name, char *dict = nullptr); // construct from name in Python scope
	void initpymat(PyObject *a); // helper function used by above
	NRmatrix(int n, int m, const T &a);	//Initialize to constant
	NRmatrix(int n, int m, const T *a);	// Initialize to array
	NRmatrix(const NRmatrix &rhs);		// Copy constructor
	NRmatrix & operator=(const NRmatrix &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	void assign(char *name, char *dict = nullptr); // assign to a Python name and scope
	~NRmatrix();
};

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), ownsdata(true), v(nullptr) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), ownsdata(true), v(n>0 ? new T*[n] : nullptr)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? (T*)PyMem_Malloc(nel*sizeof(T)) : nullptr;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
void NRmatrix<T>::initpymat(PyObject *a) {
	pyident = a;
	ownsdata = false;
	if (! PyArray_CheckExact(a)) NRpyException("PyObject is not an Array in NRmatrix constructor.");
	if (! PyArray_ISCARRAY_RO(a)) NRpyException("Python Array must be contiguous (e.g., not strided).");
	if (PyArray_NDIM(a) != 2) NRpyException("Python Array must be 2-dim in NRmatrix constructor.");
	if (! NRpyTypeOK<T>(a)) NRpyException("Python Array type does not agree with NRmatrix type.");
	int i, nel;
	nn = int(PyArray_DIMS(a)[0]);
	mm = int(PyArray_DIMS(a)[1]);
	nel = mm*nn;
	v = (nn>0 ? new T*[nn] : nullptr);
	if (v) v[0] = nel>0 ? (T*)PyArray_DATA(a) : nullptr;
	for (i=1;i<nn;i++) v[i] = v[i-1] + mm;
}
template <class T> NRmatrix<T>::NRmatrix(PyObject *a) {
	initpymat(a);
}
template <class T> NRmatrix<T>::NRmatrix(char *name, char *dict ) {
	initpymat(NRpyGetByName(name,dict));
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), ownsdata(true), v(n>0 ? new T*[n] : nullptr)
{
	const auto nel = m * n;
	if (v) v[0] = nel>0 ? (T*)PyMem_Malloc(nel*sizeof(T)) : nullptr;
	for (int i=1; i< n; i++) v[i] = v[i-1] + m;
	for (int i=0; i< n; i++) 
        for (int j=0; j<m; j++) 
            v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), ownsdata(true), v(n>0 ? new T*[n] : nullptr)
{
	const auto nel = m * n;
	if (v) v[0] = nel>0 ? (T*)PyMem_Malloc(nel*sizeof(T)) : nullptr;
	for (int i=1; i< n; i++) v[i] = v[i-1] + m;
	for (int i=0; i< n; i++) 
        for (int j=0; j<m; j++) 
            v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), ownsdata(true), v(nn>0 ? new T*[nn] : nullptr)
{
	const auto nel = mm * nn;
	if (v) v[0] = nel>0 ? (T*)PyMem_Malloc(nel*sizeof(T)) : nullptr;
	for (int i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (int i=0; i< nn; i++) 
        for (int j=0; j<mm; j++)
            v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs) {
	if (this != &rhs) {
		if (nn != rhs.nn || mm != rhs.mm) {
			resize(rhs.nn,rhs.mm);
			nn = rhs.nn;
			mm = rhs.mm;
		}
		for (int i=0; i<nn; i++) 
            for (int j=0; j<mm; j++) 
                v[i][j] = rhs[i][j];
	}
	return *this;
}

//subscripting: pointer to row i
template <class T>
inline T* NRmatrix<T>::operator[](const int i) {
    if constexpr(_CHECKBOUNDS_) {
        if (i<0 || i>=nn) {
            throw("NRmatrix subscript out of bounds");
        }
    }
	return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](const int i) const {
    if constexpr(_CHECKBOUNDS_) {
        if (i<0 || i>=nn) {
            throw("NRmatrix subscript out of bounds");
        }
    }
	return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const noexcept
{
	return nn;
}

template <class T>
inline int NRmatrix<T>::ncols() const noexcept
{
	return mm;
}

template <class T>
void NRmatrix<T>::resize(int newn, int newm) {
	
	if (newn != nn || newm != mm) {
		nn = newn;
		mm = newm;
		const int nel = mm*nn;
		if (ownsdata) {
			if (v != nullptr) {
				PyMem_Free(v[0]);
				delete[] (v);
			}
			v = nn > 0 ? new T*[nn] : nullptr;
			if (v) v[0] = nel>0 ? (T*)PyMem_Malloc(nel*sizeof(T)) : nullptr;
		} else {
			if (v != nullptr) delete[] (v);
			int dm[2];
			dm[0] = newn; dm[1] = newm;
			PyArray_Dims mydims;
			mydims.ptr = (npy_intp*)dm;
			mydims.len = 2;
			PyArray_Resize((PyArrayObject *)pyident, &mydims, 0, NPY_CORDER);
			// the return value is garbage, or maybe PyNone, contrary to Numpy docs
			// I think it's a Numpy bug, but following is correct
			v = (nn>0 ? new T*[nn] : nullptr);
			if (v) v[0] = nel>0 ? (T*)PyArray_DATA(pyident) : nullptr;
		}
		for (int i=1; i<nn; i++) v[i] = v[i-1] + mm;
	}
}

template <class T>
void NRmatrix<T>::assign(int newn, int newm, const T& a) {
	resize(newn,newm);
	for (size_t i=0; i< nn; i++) 
        for (size_t j=0; j<mm; j++) 
            v[i][j] = a;
}

template <class T>
void NRmatrix<T>::assign(char *name, char *dict) {
	if (! ownsdata) NRpyException("Attempt to assign Python matrix to another Python matrix");
	if (v != nullptr) {
		PyMem_Free(v[0]);
		delete[] (v);
	}
	initpymat(NRpyGetByName(name,dict));
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
	if (v != nullptr) {
		if (ownsdata) PyMem_Free(v[0]); // pointer to the data
		delete[] (v); // pointer to the pointers
	}
}

template <class T>
class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(nullptr) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n]) {

	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(int j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(int i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(int j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
	return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
	if (v != nullptr) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}


// basic type names (redefine if your bit lengths don't match)

// typedef int Int; // 32 bit integer
// typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

// typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

// typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

// typedef bool Bool;

// vector types

typedef const NRvector<int> VecInt_I;
typedef NRvector<int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<unsigned int> VecUint_I;
typedef NRvector<unsigned int> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<char> VecChar_I;
typedef NRvector<char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<char*> VecCharp_I;
typedef NRvector<char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<double> VecDoub_I;
typedef NRvector<double> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<double*> VecDoubp_I;
typedef NRvector<double*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<bool> VecBool_I;
typedef NRvector<bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<int> MatInt_I;
typedef NRmatrix<int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<unsigned int> MatUint_I;
typedef NRmatrix<unsigned int> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<char> MatChar_I;
typedef NRmatrix<char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<double> MatDoub_I;
typedef NRmatrix<double> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<bool> MatBool_I;
typedef NRmatrix<bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<double> Mat3DDoub_I;
typedef NRMat3d<double> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

// Python glue Part II begins here

// NRpyObject for vector and matrix
template<class T>
PyObject* NRpyObject(NRvector<T> &a) {
	if (! a.ownsdata) {Py_INCREF(a.pyident); return a.pyident;}
	npy_int nd = 1;
	npy_intp dims[1];
	dims[0] = a.size();
	PyObject *thing;
	if (dims[0] > 0) {
		thing = PyArray_SimpleNewFromData(nd, dims, NRpyDataType<T>(), &a[0]);
	} else {
		thing = PyArray_SimpleNew(nd, dims, NRpyDataType<T>()); // zero size
	}
	PyArray_FLAGS(thing) |= NPY_OWNDATA;
	a.ownsdata = false;
	a.pyident = thing;
	return thing;
}

template<class T>
PyObject* NRpyObject(NRmatrix<T> &a) {
	if (! a.ownsdata ) {Py_INCREF(a.pyident); return a.pyident;}
	npy_int nd = 2;
	npy_intp dims[2];
	dims[0] = a.nrows(); dims[1] = a.ncols();
	PyObject *thing;
	if (dims[0]*dims[1] > 0) {
		thing = PyArray_SimpleNewFromData(nd, dims, NRpyDataType<T>(), &a[0][0]);
	} else {
		thing = PyArray_SimpleNew(nd, dims, NRpyDataType<T>()); // zero size
	}
	PyArray_FLAGS(thing) |= NPY_OWNDATA;
	a.ownsdata = false;
	a.pyident = thing;
	return thing;
}

// PyObject(tuple) must go down here because it uses an NRvector
PyObject* NRpyTuple(PyObject *first, ...) {
	int nargs = 1;
	std::vector<PyObject*> argslist;
	argslist.push_back(first);
	va_list vl;
	va_start(vl,first);
	for (;;) {
		argslist.push_back(va_arg(vl,PyObject*));
		if (argslist[nargs] == nullptr) break;
		else ++nargs;
	}
	va_end(vl);
	PyObject* tuple = PyTuple_New(nargs);
	// assumes that tuple will be returned to Python, so we are not responsible for it.
	for(int i=0;i<nargs;i++) PyTuple_SetItem(tuple, i, argslist[i]);
	return tuple;
}
// If you create an NRpyTuple and don't return it, you need to Py_DECREF it later.

// macros used to make and use persistent objects

#define NRpyCONNECT(CLASS,METHOD) \
static PyObject *METHOD(PyObject *self, PyObject *pyargs) { \
	NRpyArgs args(pyargs); \
	CLASS *p = (CLASS *)PyCapsule_GetPointer(args[0],nullptr); \
	return p->METHOD(args); }
// To-Do: check that args[0] exists and is a PyCapsule before using it

// destructor to register with PyCapsule, calls actual destructor
// then, constructor calls actual constructor to create instance, returns it
#define NRpyCONSTRUCTOR(CLASS,METHOD) \
void NRpyDestroy(PyObject *myself) { \
	((CLASS*)PyCapsule_GetPointer(myself,nullptr))->~CLASS(); } \
static PyObject *METHOD(PyObject *self, PyObject *pyargs) { \
	NRpyArgs args(pyargs); \
	CLASS *instance = new CLASS(args); \
	return PyCapsule_New(instance,nullptr,NRpyDestroy); }

// functor class to help with calling Python functions from within C++ modules
template <class R> 
struct NRpyPyFunction {
	PyObject *ob;
	int argcount;
	NRpyPyFunction() : ob(nullptr), argcount(0) {}
	NRpyPyFunction(PyObject *obb) : ob(obb) {
		if(! PyCallable_Check(ob)) NRpyException("NRpyPyFunction: non-callable object.");
		PyCodeObject *code = (PyCodeObject *)PyFunction_GetCode(ob);
		argcount = code->co_argcount; // caution, uses not-officially-exposed value
	}
	inline void argcheck(int argn) {
		if (argn != argcount)
			NRpyException("NRpyPyFunction: should be %d args in PyFunction call.",1,argcount);
	}
	// constructors for 0 to 4 args.  You can add more if you want.
	R operator()() {
		argcheck(0);
		return NRpyCast<R>(PyObject_CallObject(ob,nullptr));
	}
	template <class T> R operator()(T x1) {
		PyObject *tuple;
		argcheck(1);
		tuple = NRpyTuple(NRpyObject(x1),nullptr);
		PyObject *tmp = PyObject_CallObject(ob,tuple);
		if (tmp == nullptr) NRpyException("Error in evaluating a Python function called from C++");
		return NRpyCast<R>(tmp);
	}
	template <class T, class U> R operator()(T x1, U x2) {
		PyObject *tuple;
		argcheck(2);
		tuple = NRpyTuple(NRpyObject(x1),NRpyObject(x2),nullptr);
		PyObject *tmp = PyObject_CallObject(ob,tuple);
		if (tmp == nullptr) NRpyException("Error in evaluating a Python function called from C++");
		return NRpyCast<R>(tmp);
	}
	template <class T, class U, class V> R operator()(T x1, U x2, V x3) {
		PyObject *tuple;
		argcheck(3);
		tuple = NRpyTuple(NRpyObject(x1),NRpyObject(x2),NRpyObject(x3),nullptr);
		PyObject *tmp = PyObject_CallObject(ob,tuple);
		if (tmp == nullptr) NRpyException("Error in evaluating j Python function called from C++");
		return NRpyCast<R>(tmp);
	}
	template <class T, class U, class V, class W> R operator()(T x1, U x2, V x3, W x4) {
		PyObject *tuple;
		argcheck(4);
		tuple = NRpyTuple(NRpyObject(x1),NRpyObject(x2),NRpyObject(x3),NRpyObject(x4),nullptr);
		PyObject *tmp = PyObject_CallObject(ob,tuple);
		if (tmp == nullptr) NRpyException("Error in evaluating a Python function called from C++");
		return NRpyCast<R>(tmp);
	}
};

// functor class for calling a function of one arg that may be either C or Python
// T is return type, U is argument type
template <class T, class U>
struct NRpyAnyFunction {
	T (*cfunc)(U);
	NRpyPyFunction<T> ftor;
	int ispy;
	
	NRpyAnyFunction(PyObject *ob) {
		if(PyCallable_Check(ob)) {
			ispy = 1;
			ftor = NRpyPyFunction<T>(ob);
		} else if (PyCapsule_CheckExact(ob)) {
			ispy = 0;
			NRpyCFunction(cfunc,ob);
		} else NRpyException("Not a function object of either type in NRpyAnyFunction.");
	}
	inline T operator()(U x) {
		if (ispy) return ftor(x);
		else return cfunc(x);
	}
};

#endif /* _NR3_H_ */

