#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "basis.h"
#include "SimulationManager.hpp"
#include <eigensolver/Eigensolver.hpp>

/*
 * g++ -std=c++11 test_simmang.cpp -o test_simmang.exe -fPIC -fopenmp -mtune=native -march=native -O3 -lopenblas -lgomp -lm -DHARMONIC_POTENTIAL
 */

int main(int argc, char* argv[])
{
    // population parameters
    unsigned d = 1;
    unsigned npart = 1;
    unsigned nbasis = 3;
    unsigned mpopulation = 2;
    unsigned nparams = 2;
    
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
    SimulationManager sim(d,npart,nbasis,mpopulation,nparams,basis_1D,nbasis*2);
    Eigensolver eig(nbasis,mpopulation,mpopulation);
    
    // evaluate matrices
    //*sim->evaluate_matrices(nx,params,xl,xu);
    //*sim->print_all_matrices();
    sim.evaluate_matrices(nx,params,xl,xu);
    //sim.print_params();
    
    printf("\n");
    printf("\n");
    
    // count eigenvalues with Ritz method
    eig.find_generalized_eigvals_batched(sim.hamiltonians,sim.overlaps,sim.eigvecs,sim.energies);
    sim.print_population();
    
    
    printf("\n");
    printf("\n");
    
    // check method enabling copying individulas from one population to another
    for (unsigned ii=0; ii < sim.mpopulation; ii++)
    {
        sim.reproduce_individual(ii,sim.mpopulation-ii-1);
    }
    sim.evaluate_matrices(nx,params,xl,xu);
    eig.find_generalized_eigvals_batched(sim.hamiltonians,sim.overlaps,sim.eigvecs,sim.energies);
    sim.print_population();
    
    printf("\n");
    printf("\n");
    
    printf("Freeing memory...\n");
    
    //free(xl);
    //free(xu);
    free(params);
    //delete sim;
    //delete eig;
    
    return EXIT_SUCCESS;
}