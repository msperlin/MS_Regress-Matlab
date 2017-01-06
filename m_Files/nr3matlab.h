/* nr3matlab.h */
// version 0.8
// This file is a version of nr3.h with hooks that
// make it easy to write Matlab mex files, in particular
// ones that use NR3 routines.
// See http://www.nr.com/nr3_matlab.html

#include "mex.h"

#define _CHECKBOUNDS_ 1

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <typeinfo.h>

using namespace std;

// macro-like inline functions

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

// exception handling when executing underneath Matlab

#ifdef _MSC_VER
#define throw(message) \
{char msg[1024]; sprintf_s(msg,1024,"%s in file %s at line %d\n", \
message,__FILE__,__LINE__); mexErrMsgTxt(msg);}
#else
#define throw(message) \
{char msg[1024]; sprintf(msg,"%s in file %s at line %d\n", \
message,__FILE__,__LINE__); mexErrMsgTxt(msg);}
#endif

// basic type names (put here so can use for Matlab stuff)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: you should test by verifying that (NaN != NaN) is true (see nr3.h)
static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

// get mxClassID of any type T

template <class T> inline mxClassID mxT() {return mxUNKNOWN_CLASS;}
template <> inline mxClassID mxT<Doub>() {return mxDOUBLE_CLASS;}
template <> inline mxClassID mxT<float>() {return mxSINGLE_CLASS;}
template <> inline mxClassID mxT<Int>() {return mxINT32_CLASS;}
template <> inline mxClassID mxT<Uint>() {return mxUINT32_CLASS;}
template <> inline mxClassID mxT<Char>() {return mxCHAR_CLASS;}
template <> inline mxClassID mxT<Uchar>() {return mxUINT8_CLASS;}
template <> inline mxClassID mxT<Llong>() {return mxINT64_CLASS;}
template <> inline mxClassID mxT<Ullong>() {return mxUINT64_CLASS;}
template <> inline mxClassID mxT<Bool>() {
	if (sizeof(Bool)==1) return mxLOGICAL_CLASS;
	else throw("bool and mxLOGICAL_CLASS have incompatible sizes");
}
inline mxClassID mxT(const mxArray *p) {return mxGetClassID(p);} 

// functions to map Matlab scalars

template <class T>
const T& mxScalar(const mxArray *prhs) {
	if (mxGetClassID(prhs) != mxT<T>())
		throw("attempt to assign scalar ref to wrong type");
	return *(T*)mxGetData(prhs);
}

template <class T>
T& mxScalar(mxArray* &plhs) {
	plhs = mxCreateNumericMatrix(1,1,mxT<T>(),mxREAL);
	return *(T*)mxGetData(plhs);
}

template <class T>
const T& mxScalar(const char *varname) {
	const mxArray* mxptr = mexGetVariablePtr("base",varname);
	if (mxptr == NULL) throw("attempt to get nonexistent variable");
	if (mxGetClassID(mxptr) != mxT<T>())
		throw("attempt to assign scalar ref to wrong type");
	return *(T*)mxGetData(mxptr);
}

template <class T>
void mxScalar(T val, const char *varname) {
	mxClassID mxclass = mxT<T>();
	if (mxclass == mxUNKNOWN_CLASS) throw("no corresponding Matlab type");
	mxArray* mxdum = mxCreateNumericMatrix(1,1,mxclass,mxREAL);
	*(T*)mxGetData(mxdum) = val;
	if (mexPutVariable("base",varname,mxdum))
		throw("failed to send data to Matlab variable by name");
	mxDestroyArray(mxdum);
}



// Vector and Matrix Classes

template <class T>
class NRvector {
private:
	int nn;	// size of array. upper index is nn-1
	int own; // 1 if own data, 0 if Matlab's
	T *v;
public:
	NRvector();
	explicit NRvector(int n);		// Zero-based array
	NRvector(int n, const T &a);	//initialize to constant value
	NRvector(int n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor
	NRvector & operator=(const NRvector &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	~NRvector();

	NRvector(const mxArray *prhs); // map Matlab rhs to vector (read-only)
	NRvector(int n, mxArray* &plhs); // create Matlab lhs and map to vector
	NRvector(const char *varname); // import Matlab variable by name (read-only)	
	void put(const char *varname); // export vector to a named Matlab variable

};

template <class T>
NRvector<T>::NRvector(const mxArray* prhs) : own(0) {
	if (mxGetClassID(prhs) != mxT<T>())
		throw("constructing VecDoub from a different Matlab type prhs");
	nn = mxGetNumberOfElements(prhs);
	v = (T*)mxGetData(prhs);
}

template <class T>
NRvector<T>::NRvector(int n, mxArray* &plhs) : nn(n), own(0) {
	mxClassID mxclass = mxT<T>();
	if (mxclass == mxUNKNOWN_CLASS) throw("no corresponding Matlab type for plhs");
	plhs = mxCreateNumericMatrix(1,nn,mxclass,mxREAL);
	v = (T*)mxGetData(plhs);
}

template <class T>
NRvector<T>::NRvector(const char *varname) : own(0) {
	const mxArray* mxptr = mexGetVariablePtr("base",varname);
	if (mxptr == NULL) throw("attempt to get nonexistent variable");
	if (mxGetClassID(mxptr) != mxT<T>())
		throw("constructing a VecDoub from a different Matlab type");
	nn = mxGetNumberOfElements(mxptr);
	v = (T*)mxGetData(mxptr);
}

template <class T>
void NRvector<T>::put(const char *varname) {
	mxClassID mxclass = mxT<T>();
	if (mxclass == mxUNKNOWN_CLASS) throw("no corresponding Matlab type");
	mxArray* mxdum = mxCreateNumericMatrix(1,1,mxclass,mxREAL);
	void* sav = mxGetData(mxdum);
	mxSetN(mxdum,nn);
	mxSetData(mxdum,v);
	if (mexPutVariable("base",varname,mxdum))
		throw("failed to send data to Matlab variable by name");
	mxSetData(mxdum,sav);
	mxSetN(mxdum,1);
	mxDestroyArray(mxdum);
}

template <class T>
NRvector<T>::NRvector() : nn(0), own(1), v(NULL) {}

template <class T>
NRvector<T>::NRvector(int n) : nn(n), own(1), v(n>0 ? new T[n] : NULL) {}

template <class T>
NRvector<T>::NRvector(int n, const T& a) : nn(n), own(1), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), own(1), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = *a++;
}

template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), own(1), v(nn>0 ? new T[nn] : NULL)
{
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (!own) throw("resize of mxArray by assignment not allowed");
			if (v != NULL) delete [] (v);
			nn=rhs.nn;
			v= nn>0 ? new T[nn] : NULL;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
inline T & NRvector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T & NRvector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRvector<T>::size() const
{
	return nn;
}

template <class T>
void NRvector<T>::resize(int newn)
{
	if (newn != nn) {
		if (!own) throw("resize of mxArray not allowed");
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T>
void NRvector<T>::assign(int newn, const T& a)
{
	if (newn != nn) {
		if (!own) throw("resize of mxArray by assign method not allowed");
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i=0;i<nn;i++) v[i] = a;
}

template <class T>
NRvector<T>::~NRvector()
{
	if (own && v != NULL) delete[] (v);
}

// end of NRvector definitions

template <class T>
class NRmatrix {
private:
	int nn;
	int mm;
	int own; // owned by self 1, vs Matlab 0
	T **v;
public:
	NRmatrix();
	NRmatrix(int n, int m);			// Zero-based array
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
	~NRmatrix();

	NRmatrix(const mxArray *prhs); // map Matlab rhs to matrix (read-only)
	NRmatrix(int n, int m, mxArray* &plhs); // create Matlab lhs and map to matrix
	NRmatrix(const char *varname); // import Matlab variable by name (read-only)	
	void put(const char *varname); // export matrix to a named Matlab variable

};


template <class T>
void NRmatrix<T>::put(const char *varname) {
	mxClassID mxclass = mxT<T>();
	if (mxclass == mxUNKNOWN_CLASS) throw("no corresponding Matlab type");
	mxArray* mxdum = mxCreateNumericMatrix(1,1,mxclass,mxREAL);
	void* sav = mxGetData(mxdum);
	mxSetN(mxdum,nn);
	mxSetM(mxdum,mm);
	mxSetData(mxdum,v[0]);
	if (mexPutVariable("base",varname,mxdum))
		throw("failed to send data to Matlab variable by name");
	mxSetData(mxdum,sav);
	mxSetN(mxdum,1);
	mxSetM(mxdum,1);
	mxDestroyArray(mxdum);
}

template <class T>
NRmatrix<T>::NRmatrix(const char *varname) : own(0) {
	const mxArray* mxptr = mexGetVariablePtr("base",varname);
	if (mxptr == NULL) throw("attempt to get nonexistent variable");
	if (mxGetClassID(mxptr) != mxT<T>())
		throw("constructing an NRmatrix from a different type Matlab variable");
	nn = mxGetN(mxptr);
	mm = mxGetM(mxptr);
	v = nn>0 ? new T*[nn] : NULL;
	if (v) v[0] = (T*)mxGetData(mxptr);
	for (int i=1;i<nn;i++) v[i] = v[i-1] + mm;
}

template <class T>
NRmatrix<T>::NRmatrix(const mxArray *prhs) : own(0) {
	if (mxGetClassID(prhs) != mxT<T>())
		throw("constructing NRmatrix from a different type Matlab prhs");
	nn = mxGetN(prhs);
	mm = mxGetM(prhs);
	v = nn>0 ? new T*[nn] : NULL;
	if (v) v[0] = (T*)mxGetData(prhs);
	for (int i=1;i<nn;i++) v[i] = v[i-1] + mm;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, mxArray* &plhs) : nn(n), mm(m), own(0) {
	mxClassID mxclass = mxT<T>();
	if (mxclass == mxUNKNOWN_CLASS) throw("no corresponding Matlab type for plhs");
	plhs = mxCreateNumericMatrix(mm,nn,mxclass,mxREAL);
	v = nn>0 ? new T*[nn] : NULL;
	if (v) v[0] = (T*)mxGetData(plhs);
	for (int i=1;i<nn;i++) v[i] = v[i-1] + mm;
}

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), own(1), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), own(1), v(n>0 ? new T*[n] : NULL)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), own(1), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), own(1), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), own(1), v(nn>0 ? new T*[nn] : NULL)
{
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (!own) throw("resize of mxArray by assignment not allowed");
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* NRmatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRmatrix<T>::ncols() const
{
	return mm;
}

template <class T>
void NRmatrix<T>::resize(int newn, int newm)
{
	int i,nel;
	if (newn != nn || newm != mm) {
		if (!own) throw("resize of mxArray not allowed");
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
}

template <class T>
void NRmatrix<T>::assign(int newn, int newm, const T& a)
{
	int i,j,nel;
	if (newn != nn || newm != mm) {
		if (!own) throw("resize of mxArray by assign function not allowed");
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
	if (v != NULL) {
		if (own) delete[] (v[0]);
		delete[] (v);
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
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
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
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}

// vector types

typedef const NRvector<Int> VecInt_I;
typedef NRvector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const NRvector<Uint> VecUint_I;
typedef NRvector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const NRvector<Llong> VecLlong_I;
typedef NRvector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const NRvector<Ullong> VecUllong_I;
typedef NRvector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const NRvector<Char> VecChar_I;
typedef NRvector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const NRvector<Char*> VecCharp_I;
typedef NRvector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const NRvector<Uchar> VecUchar_I;
typedef NRvector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const NRvector<Doub> VecDoub_I;
typedef NRvector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const NRvector<Doub*> VecDoubp_I;
typedef NRvector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const NRvector<Complex> VecComplex_I;
typedef NRvector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const NRvector<Bool> VecBool_I;
typedef NRvector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;
