#ifndef __BASIS_H__
#define __BASIS_H__

#include <stdlib.h>
#include <math.h>

#include <complex>

#include "ho_basis.h"
#include "planewave_basis.h"


/*
 * 
 */
std::complex<double> lorentzian_basis_func(double* x, double* params, unsigned tot_dims, unsigned nparams)
{
    const double gamma = params[0];
    const double x0    = 0.;//params[1];
    return sqrt(.25*gamma*gamma*gamma/M_PI) / ( (x[0]-x0)*(x[0]-x0) + .25*gamma*gamma );
}



typedef std::complex<double> (*ptr_func_t)(double*, double*, unsigned, unsigned);


ptr_func_t basis_1D[MAX_HO_BASIS+12] =
{   
    lorentzian_basis_func,ho_basis_func<0,1>,
    ho_basis_func<1,1>,
    ho_basis_func<2,1>,
    ho_basis_func<3,1>,
    ho_basis_func<4,1>,
    ho_basis_func<5,1>,
    ho_basis_func<6,1>,
    ho_basis_func<7,1>,
    ho_basis_func<8,1>,
    ho_basis_func<9,1>,
    ho_basis_func<10,1>,
    ho_basis_func<11,1>,
    ho_basis_func<12,1>,
    ho_basis_func<13,1>,
    ho_basis_func<14,1>,
    ho_basis_func<15,1>,
    ho_basis_func<16,1>,
    ho_basis_func<17,1>,
    ho_basis_func<18,1>,
    ho_basis_func<19,1>,
    ho_basis_func<20,1>,
    ho_basis_func<21,1>,
    ho_basis_func<22,1>,
    ho_basis_func<23,1>,
    ho_basis_func<24,1>,
    ho_basis_func<25,1>,
    ho_basis_func<26,1>,
    ho_basis_func<27,1>,
    ho_basis_func<28,1>,
    ho_basis_func<29,1>,
    ho_basis_func<30,1>,
    ho_basis_func<31,1>,
    ho_basis_func<32,1>,
    ho_basis_func<33,1>,
    ho_basis_func<34,1>,
    ho_basis_func<35,1>,
    ho_basis_func<36,1>,
    ho_basis_func<37,1>,
    ho_basis_func<38,1>,
    ho_basis_func<39,1>,
    ho_basis_func<40,1>,
    ho_basis_func<41,1>,
    ho_basis_func<42,1>,
    ho_basis_func<43,1>,
    ho_basis_func<44,1>,
    ho_basis_func<45,1>,
    ho_basis_func<46,1>,
    ho_basis_func<47,1>,
    ho_basis_func<48,1>,
    ho_basis_func<49,1>,
    ho_basis_func<50,1>,
    ho_basis_func<51,1>,
    ho_basis_func<52,1>,
    ho_basis_func<53,1>,
    ho_basis_func<54,1>,
    ho_basis_func<55,1>,
    ho_basis_func<56,1>,
    ho_basis_func<57,1>,
    ho_basis_func<58,1>,
    ho_basis_func<59,1>,
    ho_basis_func<60,1>,
    ho_basis_func<61,1>,
    ho_basis_func<62,1>,
    ho_basis_func<63,1>,
    ho_basis_func<64,1>,
    ho_basis_func<65,1>,
    ho_basis_func<66,1>,
    ho_basis_func<67,1>,
    ho_basis_func<68,1>,
    ho_basis_func<69,1>,
    ho_basis_func<70,1>,
    ho_basis_func<71,1>,
    ho_basis_func<72,1>,
    ho_basis_func<73,1>,
    ho_basis_func<74,1>,
    ho_basis_func<75,1>,
    ho_basis_func<76,1>,
    ho_basis_func<77,1>,
    ho_basis_func<78,1>,
    ho_basis_func<79,1>,
    ho_basis_func<80,1>,
    ho_basis_func<81,1>,
    ho_basis_func<82,1>,
    ho_basis_func<83,1>,
    ho_basis_func<84,1>,
    ho_basis_func<85,1>,
    ho_basis_func<86,1>,
    ho_basis_func<87,1>,
    ho_basis_func<88,1>,
    ho_basis_func<89,1>,
    ho_basis_func<90,1>,
    ho_basis_func<91,1>,
    ho_basis_func<92,1>,
    ho_basis_func<93,1>,
    ho_basis_func<94,1>,
    ho_basis_func<95,1>,
    ho_basis_func<96,1>,
    ho_basis_func<97,1>,
    ho_basis_func<98,1>,
    ho_basis_func<99,1>,
    planewave_basis_func<1,0>,
    planewave_basis_func<1,1>,
    planewave_basis_func<1,2>,
    planewave_basis_func<1,3>,
    planewave_basis_func<1,4>,
    planewave_basis_func<1,5>,
    planewave_basis_func<1,6>,
    planewave_basis_func<1,7>,
    planewave_basis_func<1,8>,
    planewave_basis_func<1,9>,
    planewave_basis_func<1,10>};


#endif