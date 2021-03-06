#ifndef __FFTW_INTERFACE_H__
#define __FFTW_INTERFACE_H__

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include <fftw3.h>

#include "common_header.h"

// creating plan
void create_plan_1d(size_t xdim);

// applying transforms
void fft_Z2Z_forward_1d(Cplx *in_arr, Cplx *out_arr, const size_t xdim);
void fft_Z2Z_inverse_1d(Cplx *in_arr, Cplx *out_arr, const size_t xdim);

// applying differentiation
void fftdiff_1d(double *in_arr, double *dfdx, double *d2fd2x, const size_t xdim);

// destroying plans & cleaning memory
void destroy_plan(void);

#endif
