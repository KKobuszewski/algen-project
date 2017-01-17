#ifndef __HO_BASIS_H__
#define __HO_BASIS_H__

#include <stdlib.h>
#include <math.h>

#include <complex>
#include <boost/math/special_functions/hermite.hpp>

// #include "basis.h"
#include "math_extensions.h"


/*
 * This function could be use only for one-particle case! (due to optimizations)
 * 
 * TODO: IF set dims = 0/-1 -> uses tot_dims and not templated dims !!!
 * 
 * @param double* x         - array containing point in phase space
 * @param double* params    - array containing parameters of wavefunction
 * @param unsigned tot_dims - number of elements in array x 
 * @param unsigned nparams  - number of elements in array params
 */
template <unsigned n, unsigned dims> std::complex<double> ho_basis_func(double* x, double* params, unsigned tot_dims, unsigned nparams)
{
    // tot_dims ignored!
    switch (dims)
    {
    case 1:
    {
        return   boost::math::hermite(n, x[0]) * exp(-.5*x[0]*x[0]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) );
    }
    /*case 2:
    { 
        return   boost::math::hermite(n, x[0]) * exp(-.5*x[0]*x[0]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) )
               * boost::math::hermite(n, x[1]) * exp(-.5*x[1]*x[1]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) );
    }
    case 3:
    {
        return   boost::math::hermite(n, x[0]) * exp(-.5*x[0]*x[0]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) )
               * boost::math::hermite(n, x[1]) * exp(-.5*x[1]*x[1]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) )
               * boost::math::hermite(n, x[2]) * exp(-.5*x[2]*x[2]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) );
    }*/
    default:
    {
        // TODO: Implement for cases of many particles!!! -> permanent/determinant
        double ret = 1.;
        for (unsigned ii = 0; ii < tot_dims; ii++)
        {
            ret *= boost::math::hermite(n, x[ii]) * exp(-.5*x[ii]*x[ii]) / sqrt( pow(2,n)*factorial(n)*sqrt(M_PI) );
        }
        return ret;
    }
    }
}


#define MAX_HO_BASIS 100


/*
template<unsigned n>
struct ho_basis
{
    ptr_func_t ho_basis_func<n>
}*/

#endif