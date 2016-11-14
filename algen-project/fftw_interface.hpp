#ifndef __FFTW_INTERFACE_HPP__
#define __FFTW_INTERFACE_HPP__

#include <complex>
#include <fftw3.h>
#include <cmath>

// creating plan
void create_plan_1d(size_t xdim);

// applying transforms
void fft_Z2Z_forward_1d(std::complex<double> *, std::complex<double> *, const size_t);
void fft_Z2Z_inverse_1d(std::complex<double> *, std::complex<double> *, const size_t);

// applying differentiation
void fftdiff_1d(double *in_arr, double *dfdx, double *d2fd2x, const size_t xdim);

// destroying plans & cleaning memory
void destroy_plan(void);

#endif