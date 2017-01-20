#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "fftw_interface.hpp"

#ifndef M_PI
#define M_PI 3.141592653589793
#endif



fftw_complex *fftw_input = NULL;
fftw_complex *fftw_input2 = NULL;
fftw_complex *fftw_output = NULL;
fftw_complex *fftw_output2 = NULL;
fftw_plan pf, pi, pi2;


void create_plan_1d(size_t xdim)
{
    fftw_input   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim);
    fftw_input2  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim);
    fftw_output  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim);
    fftw_output2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim);

    pf  = fftw_plan_dft_1d(xdim, fftw_input, fftw_output, FFTW_FORWARD,  FFTW_PATIENT);
    pi  = fftw_plan_dft_1d(xdim, fftw_output, fftw_input, FFTW_BACKWARD, FFTW_PATIENT);
    pi2 = fftw_plan_dft_1d(xdim, fftw_output2, fftw_input2, FFTW_BACKWARD, FFTW_PATIENT);
}

void create_plan_2d(size_t xdim, size_t ydim)
{
    fftw_input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim*ydim);
    fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim*ydim);

    pf = fftw_plan_dft_2d(xdim, ydim, fftw_input, fftw_output, FFTW_FORWARD,  FFTW_PATIENT);
    pi = fftw_plan_dft_2d(xdim, ydim, fftw_output, fftw_input, FFTW_BACKWARD, FFTW_PATIENT);
}

void create_plan_3d(size_t xdim, size_t ydim, size_t zdim)
{
    fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim*ydim*zdim);
    fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xdim*ydim*zdim);

    pf = fftw_plan_dft_3d(xdim, ydim, zdim, fftw_input, fftw_output, FFTW_FORWARD,  FFTW_PATIENT);
    pi = fftw_plan_dft_3d(xdim, ydim, zdim, fftw_output, fftw_input, FFTW_BACKWARD, FFTW_PATIENT);
}

void fft_Z2Z_forward_1d(std::complex<double> *in_arr_cpp, std::complex<double> *out_arr_cpp, const size_t xdim)
{
    fftw_complex* in_arr  = reinterpret_cast<fftw_complex*>(in_arr_cpp );
    fftw_complex* out_arr = reinterpret_cast<fftw_complex*>(out_arr_cpp);
    
    fftw_plan plan = fftw_plan_dft_1d(xdim, in_arr, out_arr, FFTW_FORWARD,  FFTW_ESTIMATE);
    fftw_execute(plan);
}

void fft_Z2Z_inverse_1d(std::complex<double> *in_arr_cpp, std::complex<double> *out_arr_cpp, const size_t xdim)
{
    fftw_complex* in_arr  = reinterpret_cast<fftw_complex*>(in_arr_cpp );
    fftw_complex* out_arr = reinterpret_cast<fftw_complex*>(out_arr_cpp);
    
    fftw_plan plan = fftw_plan_dft_1d(xdim, in_arr, out_arr, FFTW_BACKWARD,  FFTW_ESTIMATE);
    fftw_execute(plan);
}

void fft_D2Z_forward_1d(double *in_arr, std::complex<double> *out_arr_cpp, const size_t xdim)
{
    std::complex<double>* cpp_input  = reinterpret_cast<std::complex<double>*>(fftw_input );
    std::complex<double>* cpp_output = reinterpret_cast<std::complex<double>*>(fftw_output);
    
    #pragma omp simd
    for (size_t ix = 0; ix < xdim; ix++)  { std::complex<double> tmp(in_arr[ix], 0.); cpp_input[ix] = tmp; }
    
    fftw_execute(pf);
    
    #pragma omp simd
    for (size_t ix = 0; ix < xdim; ix++)    out_arr_cpp[ix] = cpp_output[ix];
}

void fft_Z2D_inverse_1d(std::complex<double> *in_arr_cpp, double *out_arr, const size_t xdim)
{
    fftw_complex* in_arr  = reinterpret_cast<fftw_complex*>(in_arr_cpp);
    
    fftw_plan plan = fftw_plan_dft_1d(xdim, in_arr, in_arr, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    
    #pragma omp simd
    for (size_t ix = 0; ix < xdim; ix++)    out_arr[ix] = in_arr_cpp[ix].real();
}

void fftdiff_1d(double *in_arr, double *dfdx, double *d2fd2x, const size_t xdim)
{
    size_t ix;
    int jx;
    
    std::complex<double>* cpp_input   = reinterpret_cast<std::complex<double>*>(fftw_input  );
    std::complex<double>* cpp_input2  = reinterpret_cast<std::complex<double>*>(fftw_input2 );
    std::complex<double>* cpp_output  = reinterpret_cast<std::complex<double>*>(fftw_output );
    std::complex<double>* cpp_output2 = reinterpret_cast<std::complex<double>*>(fftw_output2);
    
    // copy input array to fftw input (double to complex)
    #pragma omp simd
    for (ix=0; ix < xdim; ix++) { std::complex<double> tmp(in_arr[ix], 0.); cpp_input[ix] = tmp; }
    
    
    // execute fftw forward
    fftw_execute(pf);
    
    // copy fourier transform
    memcpy(fftw_output2, fftw_output, xdim * sizeof(std::complex<double>));
    
    #pragma omp parallel sections default(shared) private(ix,jx) num_threads(2)
    {
    #pragma omp section
    {
        /* *** first derivative *** */
        // multiply by i*k_x / xdim (division by number of samples in case to normalize)
        #pragma omp simd
        for ( ix = 0; ix < xdim / 2; ix++ )
        {
            std::complex<double> tmp( 0., (2. * ( double ) M_PI / ((double) xdim) * ( double ) ix) / ((double) xdim) );
            cpp_output[ix] *= tmp;
        }
        jx = - ix;
        #pragma omp simd
        for ( ix = xdim / 2; ix < xdim; ix++ ) 
        {
            std::complex<double> tmp( 0., (2. * ( double ) M_PI / ((double) xdim) * ( double ) jx) / ((double) xdim) );
            cpp_output[ix] *= tmp;
            jx++;
        }
        
        // execute fftw inverse
        fftw_execute(pi);
        
        // copy result to double* array
        #pragma omp simd
        for (ix=0; ix < xdim; ix++) dfdx[ix] = cpp_input[ix].real();
    }
    #pragma omp section
    {
        /* *** second derivative *** */
        // multiply by -k_x*k_x / xdim (division by number of samples in case to normalize)
        #pragma omp simd
        for ( ix = 0; ix < xdim / 2; ix++ )
        {
            cpp_output2[ix] *= -1. * (2. * ( double ) M_PI / ((double) xdim) * ( double ) ix)
                                   * (2. * ( double ) M_PI / ((double) xdim) * ( double ) ix) / ((double) xdim);
        }
        jx = - ix;
        #pragma omp simd
        for ( ix = xdim / 2; ix < xdim; ix++ )
        {
            cpp_output2[ix] *= -1. * (2. * ( double ) M_PI / ((double) xdim) * ( double ) jx)
                                   * (2. * ( double ) M_PI / ((double) xdim) * ( double ) jx) / ((double) xdim);
            jx++;
        }
        
        // execute fftw inverse
        fftw_execute(pi2);
        
        // copy result to double* array
        #pragma omp simd
        for (ix=0; ix < xdim; ix++) d2fd2x[ix] = cpp_input2[ix].real();
    }
    }
}

void destroy_plan(void)
{
    #pragma omp parallel num_threads(4) shared(pf,pi,pi2,fftw_input,fftw_output,fftw_output2)
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                fftw_destroy_plan(pf);
                if (fftw_input)   free(fftw_input);
            }
            #pragma omp section
            {
                fftw_destroy_plan(pi);
                if (fftw_output)  free(fftw_output);
            }
            #pragma omp section
            {
                if (fftw_output2) {  fftw_destroy_plan(pi2); free(fftw_output2);  }
            }
            #pragma omp section
            {
                if (fftw_input2)  free(fftw_input2);
            }
        }
    }
}
