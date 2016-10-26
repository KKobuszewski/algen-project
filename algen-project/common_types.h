#pragma once

#ifndef _COMMON_TYPES_H_
#define _COMMON_TYPES_H_
///////////////////////////////////////////////////////////////////////////////////////////////
// Include:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// External Libs:
#include <fftw3.h>
#include "common_defined.h"

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// DEFINE:
///////////////////////////////////////////////////////////////////////////////////////////////

/* * * * * * * * * * * * * * * * *
* Real & Imag
*/
#define Real 0
#define Imag 1

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// TYPES:
///////////////////////////////////////////////////////////////////////////////////////////////

/* * * * * * * * * * * * * * * * *
* Default major data type.
* float / double / long double
*/
typedef double Data;

/* * * * * * * * * * * * * * * * *
* Default complex type.
*/
class Cplx_core {
private:
	fftw_complex fftw_cplx;
public:
	INLINE fftw_complex& get();
public:
	INLINE double&    operator [] (int);
	INLINE Cplx_core& operator =  (Cplx_core);
	INLINE Cplx_core& operator += (Cplx_core);
	INLINE Cplx_core& operator -= (Cplx_core);
public:
	INLINE Cplx_core& operator + (Cplx_core);
	INLINE Cplx_core& operator - (Cplx_core);
	INLINE Cplx_core& operator * (Cplx_core);
	INLINE Cplx_core& operator / (Cplx_core);
public:
	INLINE Cplx_core& operator * (float);
	INLINE Cplx_core& operator / (float);
};
typedef Cplx_core Cplx;

/* * * * * * * * * * * * * * * * *
* Basis functions.
* They take arguments: double* point, double* arguments, ndims, nparams
* They return value of type double __complex__
*/
typedef void (*basis_ptr) (Cplx*, double*, double*, unsigned, unsigned);


///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// CPLX TYPE SUPPORT:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// basic macros CORE version:

// USE WITH CAUTION!!!
#define CPLX_MULTIPLY_CORE(_out, _z1, _z2)													  \
	_out[Real] = _z1[Real] * _z2[Real] - _z1[Imag] * _z2[Imag];								  \
	_out[Imag] = _z1[Real] * _z2[Imag] + _z1[Imag] * _z2[Real];

// USE WITH CAUTION!!!
#define CPLX_DIVIDE_CORE(_out, _z1, _z2)													  \
	_out[Real] = (_z1[Real] * _z2[Real] + _z1[Imag] * _z2[Imag]) / (_z2[Real] * _z2[Real] + _z2[Imag] * _z2[Imag]); \
	_out[Imag] = (_z1[Imag] * _z2[Real] - _z1[Real] * _z2[Imag]) / (_z2[Real] * _z2[Real] + _z2[Imag] * _z2[Imag]);

///////////////////////////////////
// basic macros standard version:

#define CPLX_ASSIGN(_out, _z) {																  \
	_out[Real] = _z[Real];																	  \
	_out[Imag] = _z[Imag];																	  \
}

#define CPLX_ADD_ASSIGN(_out, _z) {															  \
	_out[Real] += _z[Real];																	  \
	_out[Imag] += _z[Imag];																	  \
}

#define CPLX_ADD(_out, _z1, _z2) {															  \
	_out[Real] = _z1[Real] + _z2[Real];														  \
	_out[Imag] = _z1[Imag] + _z2[Imag];														  \
}

#define CPLX_SUBTRACT(_out, _z1, _z2) {														  \
	_out[Real] = _z1[Real] - _z2[Real];														  \
	_out[Imag] = _z1[Imag] - _z2[Imag];														  \
}

#define CPLX_MULTIPLY(_out, _z1, _z2) {														  \
	Cplx _tmp_out;																			  \
	CPLX_MULTIPLY_CORE(_tmp_out, _z1, _z2);													  \
	_out[Real] = _tmp_out[Real];															  \
	_out[Imag] = _tmp_out[Imag];															  \
}

#define CPLX_DIVIDE(_out, _z1, _z2) {														  \
	Cplx _tmp_out;																			  \
	CPLX_DIVIDE_CORE(_tmp_out, _z1, _z2);													  \
	_out[Real] = _tmp_out[Real];															  \
	_out[Imag] = _tmp_out[Imag];															  \
}

#define CPLX_SCALE(_out, _z, _scale) {														  \
	_out[Real] = _z[Real] * _scale;															  \
	_out[Imag] = _z[Imag] * _scale;															  \
}

#define CPLX_CONJUGATE(_z) {																  \
	_z[Imag] = (-1) * _z[Imag];																  \
}

#define CPLX_MODULE_2(_z)																	  \
	_z[Real] * _z[Real] + _z[Imag] * _z[Imag]

#define CPLX_MODULE(_z)																		  \
	sqrt( CPLX_MODULE_2(_z) )

///////////////////////////////////
// overloaded class operators (Cplx - Cplx only):

INLINE double& Cplx_core::operator [] (int i) {
	return fftw_cplx[i];
}

INLINE Cplx_core& Cplx_core::operator = (Cplx_core z) {
	CPLX_ASSIGN(this->fftw_cplx, z.fftw_cplx);
	return *this;
}

INLINE Cplx_core& Cplx_core::operator += (Cplx_core z) {
	CPLX_ADD_ASSIGN(this->fftw_cplx, z.fftw_cplx);
	return *this;
}

INLINE Cplx_core& Cplx_core::operator -= (Cplx_core z) {
	this->fftw_cplx[Real] -= z.fftw_cplx[Real];
	this->fftw_cplx[Imag] -= z.fftw_cplx[Imag];
	return *this;
}

///////////////////////////////////
// overloaded math operators (Cplx - Cplx only):

INLINE Cplx_core& Cplx_core::operator + (Cplx_core z) {
	CPLX_ADD(this->fftw_cplx, this->fftw_cplx, z);
	return *this;
}

INLINE Cplx_core& Cplx_core::operator - (Cplx_core z) {
	CPLX_ADD(this->fftw_cplx, this->fftw_cplx, z);
	return *this;
}

INLINE Cplx_core& Cplx_core::operator * (Cplx_core z) {
	CPLX_MULTIPLY(this->fftw_cplx, this->fftw_cplx, z);
	return *this;
}

INLINE Cplx_core& Cplx_core::operator / (Cplx_core z) {
	CPLX_DIVIDE(this->fftw_cplx, this->fftw_cplx, z);
	return *this;
}

///////////////////////////////////
// overloaded math operators: (primitive_type - Cplx):

INLINE Cplx_core& Cplx_core::operator * (float scale) {
	CPLX_SCALE(this->fftw_cplx, this->fftw_cplx, scale);
	return *this;
}

INLINE Cplx_core& operator * (float scale, Cplx_core& z) {
	CPLX_SCALE(z, z, scale);
	return z;
}

INLINE Cplx_core& Cplx_core::operator / (float scale) {
	scale = 1.0f / scale;
	CPLX_SCALE(this->fftw_cplx, this->fftw_cplx, scale);
	return *this;
}

INLINE Cplx_core& operator / (float scale, Cplx_core& z) {
	Cplx_core reverse;
	reverse[Real] = scale;
	reverse[Imag] = 0;
	CPLX_DIVIDE(z, reverse, z);
	return z;
}

///////////

#endif

