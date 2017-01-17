#ifndef __HYPERGOEM_BASIS_H__
#define __HYPERGOEM_BASIS_H__

#include <stdlib.h>
#include <math.h>

#include <complex>
#include <gsl_sf_hyperg.h>

// #include "basis.h"
#include "math_extensions.h"

/*
 * params[0] - c
 * params[1] - normalization
 */
std::complex<double> basis_1D_0_F_1(double* x, double* params, unsigned tot_dims, unsigned nparams)
{
    std::complex<double> res( double gsl_sf_hyperg_0F1(params[0], x[0]), 0. );
    return params[1]*res;
}

std::complex<double> basis_1D_1_F_1
{
    std::complex<double> res( gsl_sf_hyperg_1F1 (params[0], params[1], double x), 0. );
}




#endif