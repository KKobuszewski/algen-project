#ifndef __PLANEWAVE_BASIS_H__
#define __PLANEWAVE_BASIS_H__

#include <stdlib.h>
#include <error.h>
#include <math.h>

#include <complex>

// #include "basis.h"
#include "math_extensions.h"

// inline unsigned int factorial(const unsigned int n)
// {
//     unsigned int retval = 1;
//     for (int i = n; i > 1; --i) retval *= i;
//     return retval;
// }

/*
 * This function could be use only for one-particle case! (due to optimizations)
 * 
 * 
 */
template <unsigned dims, unsigned n> 
std::complex<double> planewave_basis_func(double* x, double* params, unsigned tot_dims, unsigned nparams)
{
    // tot_dims ignored!
    switch (dims)
    {
    case 1:
    {
        const std::complex<double> exponent(0.,2.*M_PI*params[0]*x[0]/params[1]);
        return   exp( exponent ) / sqrt( params[1] );
    }
    case 2:
    {
        const std::complex<double> exponent(0.,2.*M_PI*(params[0]*x[0] + params[1]*x[1])/params[2]);
        return   exp( exponent ) / sqrt( params[2] );
    }
    case 3:
    {
        const std::complex<double> exponent(0.,2.*M_PI*(params[0]*x[0] + params[1]*x[1] + params[2]*x[2])/params[3]);
        return   exp( exponent ) / sqrt( params[3] );
    }
    default:
    {
        if (dims != nparams-1) { fprintf(stderr,"\nERROR! (planewave_basis_func): Wrong number of parameters!\n\n"); exit(EXIT_FAILURE); }
        std::complex<double> res = 1.;
        for (unsigned ii = dims-1; ii == 0; ii--)
        {
            const std::complex<double> exponent(0.,2.*M_PI*params[ii]*x[ii]/params[nparams-1]);
            res = res*exp(exponent);
        }
        return res / sqrt( params[nparams-1] );
    }
    }
}

#endif