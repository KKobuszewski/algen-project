#pragma once

#ifndef _COMMON_TYPES_H_
#define _COMMON_TYPES_H_
///////////////////////////////////////////////////////////////////////////////////////////////
// Include:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// External Libs:
#include <fftw3.h>

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// Types:
///////////////////////////////////////////////////////////////////////////////////////////////

/* * * * * * * * * * * * * * * * *
* Default major data type.
* float / double / long double
*/
typedef double Data;

/* * * * * * * * * * * * * * * * *
* Default complex type.
*/
typedef fftw_complex Cplx;

/* * * * * * * * * * * * * * * * *
* Basis functions.
* They take arguments: double* point, double* arguments, ndims, nparams
* They return value of type double __complex__
*/
typedef void (*basis_ptr) (Cplx*, double*, double*, unsigned, unsigned);

///////////

#endif
