#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "basis.h"
#include "SimulationManager.hpp"
#include <eigensolver/Eigensolver.hpp>
#include "genetic.hpp"

/*
 * g++ -std=c++11 test_simmang.cpp -o test_simmang.exe -fPIC -fopenmp -mtune=native -march=native -O3 -lopenblas -lgomp -lm -DHARMONIC_POTENTIAL
 */

inline void mean_energy(SimulationManager* sim, double* mean, double* varaince)
{
    const unsigned m = sim->mpopulation;
    double _mean = 0.;
    double _variance = 0.;
    #pragma omp simd
    for (unsigned ii=0; ii<m; ii++)
    {
        _mean += sim->energies[ii][0];
    }
    _mean /= m;
    
    #pragma omp simd
    for (unsigned ii=0; ii<m; ii++)
    {
        _variance += pow(_mean - sim->energies[ii][0],2);
    }
    _variance /= m;
    
    *mean = _mean;
    *varaince = _variance;
}


int main(int argc, char* argv[])
{
    // genetic algorithm parameters
    double pc = 0.7;
    double pm = 0.1;
    
    // population parameters
    unsigned d = 1;                // dimensions of space
    unsigned npart = 1;            // number of particles
    unsigned nbasis = 4;           // number of basis functions for each individual
    unsigned mpopulation = 2;     // number of individuals
    unsigned nparams = 2;          // number of parameters of basis function (maximal)
    
    // system parameters
    unsigned nx = 8196;
    double L = 32.;
    double xl[1];
    double xu[2];
    //double* xl = (double*) malloc( d*sizeof(double) );
    //double* xu = (double*) malloc( d*sizeof(double) );
    double* params = (double*) malloc( nparams*sizeof(double) );
    xl[0] = -L/2.;
    xu[0] =  L/2.;
    params[0] = 1.;
    params[1] = 1.;
    
    // create simulation
    //*SimulationManager* sim = new SimulationManager(d,npart,nbasis,mpopulation,nparams,basis_1D,nbasis);
    SimulationManager sim(d,npart,nbasis,mpopulation,nparams,basis_1D,nbasis*3);
    Eigensolver eig(nbasis,mpopulation,mpopulation);
    
    // evaluate matrices
    sim.evaluate_matrices(nx,params,xl,xu);
    
    printf("\n");
    printf("\n");
    
    // count eigenvalues with Ritz method
    eig.find_generalized_eigvals_batched(sim.hamiltonians,sim.overlaps,sim.eigvecs,sim.energies);
    
    
    
    // find groundstate with genetic algorithm;
    size_t tot_steps = 20LU;
    double e_max, var;
    double e_max_0, var_0;
    double dE;
    mean_energy(&sim, &e_max_0, &var_0);
    printf("%u. \nMean E: %.5lf +/- %.5lf \n",0,e_max_0,var_0);
    e_max = e_max_0;
    sim.print_basis_function_histogram();
    printf("\n");
    printf("\n");
    
    
    
    sim.print_functions();
    sim.swap_basis_func(0,1,0,2);
    
    sim.print_functions();
    sim.swap_basis_func(1,0,1,2);
    
    
    sim.print_functions();
    
    
    return EXIT_SUCCESS;
}