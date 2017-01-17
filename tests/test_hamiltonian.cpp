#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>


#include "hamiltonian.h"
#include "basis.h"


/*
 * Compile: g++ -std=c++11 test_hamiltonian.cpp -o test_hamiltonian.exe -fPIC -fopenmp -mtune=native -march=native -O3 -lgomp -lm -DHARMONIC_POTENTIAL
 * 
 * 
 */


//template<typename T1, typename T2>
//struct mc_integrand_types  {  typedef T2 (*op_integrand_t)(T1*, unsigned);  };



inline double sum_array(double* arr, size_t size)
{
    double sum = 0.;
    #pragma omp simd reduction (+:sum)
    for (size_t ii = 0; ii < size; ii++)
        sum += arr[ii];
    return sum;
}


inline std::complex<double> sum_array(std::complex<double>* arr, size_t size)
{
    std::complex<double> sum = 0.;
    //#pragma omp simd reduction (+:sum)
    for (size_t ii = 0; ii < size; ii++)
        sum += arr[ii];
    return sum;
}


#define MAX_BASIS 5

//#pragma gcc pop_options
//#pragma gcc optimize("unroll-loops")
int main(int argc, char* argv[])
{
    // wavefunction parameters
    const unsigned ndims = 1;
    const unsigned npart = 1;
    const unsigned nbasis = MAX_BASIS;
    double *params = NULL;
    unsigned nparams = 0;
    
    // grid
    const size_t nx = 4*4096;
    const double xmin = -16.;
    const double xmax = 16.;
    const double dx = (xmax - xmin)/((double) nx);
    //double *x = (double*) malloc( nx*sizeof(double) );
    
    double* Tpsi = (double*) malloc( nx*sizeof(double) );
    double* Vpsi = (double*) malloc( nx*sizeof(double) );
    std::complex<double>* Hpsi = (std::complex<double>*) malloc( nx*sizeof(std::complex<double>) );
    
    ptr_func_t psi[MAX_BASIS] = {ho_basis_func<0,1>,ho_basis_func<1,1>,ho_basis_func<2,1>,ho_basis_func<3,1>,ho_basis_func<4,1>};
    
    for (unsigned ii=0; ii<nbasis; ii++)
    {
        double ekin, epot, etot;
        
        #pragma omp parallel sections num_threads(3)
        {
        #pragma omp section
        {
            for(unsigned ix=0; ix < nx; ix++)
            {
                double x[1] = {xmin + ix*dx};
                Tpsi[ix] = ( std::conj(psi[ii](x,params,ndims*npart,nparams)) * 
                             kinetic(psi[ii], x, params, ndims, nparams, npart) ).real();
            }
            ekin = dx*sum_array(Tpsi,nx);
        }
        #pragma omp section
        {
            //ptr_func_t psi = ho_basis_func<ii>;
            //ptr_func_t psi = template std::complex<double> ho_basis_func<ii>(double* x, double* params, unsigned ndims, unsigned nparams);
            for(unsigned ix=0; ix < nx; ix++)
            {
                double x[1] = {xmin + ix*dx};
                Vpsi[ix] = ( std::conj(psi[ii](x,params,ndims*npart,nparams)) * 
                             psi[ii](x,params,ndims*npart,nparams) * potential(x, ndims, npart) ).real();
            }
            epot = dx*sum_array(Vpsi,nx);
        }
        #pragma omp section
        {
//             ptr_func_t psi = ho_basis_func<ii>;
//             ptr_func_t psi = template std::complex<double> ho_basis_func<ii>(double* x, double* params, unsigned ndims, unsigned nparams);
            for(unsigned ix=0; ix < nx; ix++)
            {
                double x[1] = {xmin + ix*dx};
                Hpsi[ix] = hamiltonian_elem(psi[ii],psi[ii], x, params, ndims, nparams, npart);
            }
            etot = dx*sum_array(Hpsi,nx).real();
        }
        }
        
        printf("======= N=%u =======================================\n",ii);
        printf("\n");
        printf("ekin:   %.10lf\n",ekin);
        printf("epot:   %.10lf\n",epot);
        printf("etot:   %.10lf\n",etot);
        printf("virial: %.10lf\n",2*epot - 2*ekin);
    }
    
    
    //free(x);
    free(Tpsi);
    free(Vpsi);
    free(Hpsi);
    
    return EXIT_SUCCESS;
}
//#pragma gcc push_options