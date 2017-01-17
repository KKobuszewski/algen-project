#ifndef __HAMILTONIAN_H__
#define __HAMILTONIAN_H__

#include <stdlib.h>
#include <stdio.h>
#include <error.h>
#include <math.h>
#include <omp.h>

#include <complex>

#include "basis.h"

#define DX ((double) 9e-6)
#define MAX_PARTICLE 16
#define MAX_DIM      4



#ifdef HARMONIC_POTENTIAL
double omega[MAX_DIM] = {1.,1.,1.};

inline double potential(double* x, const unsigned ndims, unsigned nparticles)
{
    //if (pot_nparams == 0) fprintf(stderr,"WARNING: (basis.h/potential)\tWrong number of params!\n");
    
    double pot = 0.;
    #pragma omp simd reduction(+:pot)
    for (unsigned ii=0; ii < nparticles*ndims; ii++)
    {
        pot += omega[ii%ndims]*x[ii]*x[ii] + .05*omega[ii%ndims]*x[ii]*x[ii]*omega[ii%ndims]*x[ii]*x[ii];
    }
    return .5*pot;  // 1/2 m w^2 \sum_i \sum_j x_ij^2
}
#endif
#ifdef PLAIN
inline double potential(double* x, const unsigned ndims, unsigned nparticles) {return 0.;}
#endif


#ifndef SPIN // spinless hamiltonian
inline std::complex<double> kinetic(ptr_func_t psi, double* x, double* params, unsigned ndims, unsigned nparams, unsigned nparticles)
{
    std::complex<double> T;
    double tmp_x[MAX_DIM*MAX_PARTICLE];
    
    #pragma omp simd
    for (int ii=0; ii<ndims*nparticles; ii++) tmp_x[ii] = x[ii];
    
    // evalute \psi(...,\vec{x}_i - dx*\hat{j},...) - 2\psi(...,\vec{x}_i,...) + \psi(...,\vec{x}_i - dx*\hat{j},...) <- no spin
    //#pragma omp simd reduction(+:T)
    for (unsigned ii=0; ii < ndims*nparticles; ii++)
    {
        tmp_x[ii] = x[ii] - DX;
        T += psi(tmp_x, params, ndims*nparticles, nparams) - 2.*psi(x, params, ndims*nparticles, nparams);
        tmp_x[ii] = x[ii] + DX;
        T += psi(tmp_x, params, ndims*nparticles, nparams);
    }
    
    //printf("T: %lf + %lfj\n",T.real(),T.imag());
    
    return -.5*T/DX/DX; // - \hbar^2/2m * 1/dx^2 \sum_i \nabla_i \psi(\vec{x}_1,\vec{x}_2, ..., \vec{x}_i, ..., vec{x}_Npart)
}
#endif



inline std::complex<double> hamiltonian_elem(ptr_func_t psi_i,ptr_func_t psi_j, double* x, double* params, unsigned ndims, unsigned nparams, unsigned nparticles)
{
    std::complex<double> h(0.,0.);
    h = kinetic(psi_j, x, params, ndims, nparams, nparticles) + potential(x, ndims, nparticles)*psi_j(x, params, ndims*nparticles, nparams);
    
    return std::conj(psi_i(x, params, ndims*nparticles, nparams)) * h;
}

inline double real_hamiltonian_elem(ptr_func_t psi_i,ptr_func_t psi_j, double* x, double* params, unsigned ndims, unsigned nparams, unsigned nparticles)
{
    return hamiltonian_elem(psi_i,psi_j,x,params,ndims,nparams,nparticles).real();
}


inline std::complex<double> identity(ptr_func_t psi_i,ptr_func_t psi_j, double* x, double* params, unsigned ndims, unsigned nparams, unsigned nparticles)
{
    return std::conj(psi_i(x, params, ndims*nparticles, nparams))*psi_j(x, params, ndims*nparticles, nparams);
}

#endif