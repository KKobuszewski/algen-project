/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h> // == <ccomplex>, std::complex<typenaame T> <- <complex>
#include <math.h>

#include <omp.h>

#include <random>


/*
 * Basis functions.
 * They take arguments: double* point, double* arguments, ndims, nparams
 * They return value of type double __complex__
 */
typedef (double __complex__) (*basis_ptr)(double*, double*, unsigned, unsigned);


int main()
{
  const unsigned N = 8;                  // basis size
  const unsigned m = 100;                // population size
  const unsigned d = 1;                  // number of dimensions
  const unsigned npart = 1;              // number of particles

  printf("\n");

  // allocate array of pointers to a function for basis representation
  basis_ptr* ptr_func_arr = (basis_ptr*) malloc( sizeof(basis_ptr) * N );

  // allocate array of pointers to hamiltonians
  double** hamiltonians    = (double**) malloc( sizeof(double*) * m );
  double*  hamiltonian_mem = (double*) malloc( sizeof(double) * m*N*N );
  for (unsigned ii = 0; ii < m; ii++) hamiltonians[ii] = &(hamiltonian_mem[ii*N]);

  // allocate array of pointers to hamiltonians
  double** overlaps     = (double**) malloc( sizeof(double*) * m );
  double*  overlaps_mem = (double*) malloc( sizeof(double) * m*N*N );
  for (unsigned ii = 0; ii < m; ii++) overlaps[ii] = &(overlaps_mem[ii*N]);

  // allocate array of pointers to basis functions coefficients
  double** coeffs     = (double**) malloc( sizeof(double*) * m );
  double*  coeffs_mem = (double*) malloc( sizeof(double) * m*N );
  for (unsigned ii = 0; ii < m; ii++) coeffs[ii] = &(coeffs_mem[ii*N]);
  printf("Total bytes allocated: %lu\n", 3*sizeof(double*)*m + sizeof(double)*m*N*(2*N+1) );
  printf("\n");

  // generate random population
  printf("Initializing population ...\n");
  std::random_device r; // Seed with a real random value, if available
  std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
  std::mt19937_64 mersenne(seed);
  std::uniform_real_distribution<double>uniform_dist(0., 1.);
  register double norm_const = 0;
  //#pragma omp parallel for
  for(unsigned ii = 0; ii < m; ii++)
  {
      #pragma omp simd
      for (unsigned jj = ii; jj < ii*N; jj++)
      {
          register double random = uniform_dist(mersenne);
          coeffs_mem[jj] =  random;
          norm_const += random;
      }
      #pragma omp simd
      for (unsigned jj = ii; jj < ii*N; jj++) coeffs_mem[jj] /= norm_const;
      norm_const = 0;
  }
  printf("\n");


  for (unsigned ii = 0; ii < m; ii++)
  {
      printf("%u individual:\t",ii);
      for (unsigned jj = 0; jj < m; ii++)
      {
          printf("%.5lf ", coeffs_mem[ii*N + jj]);
      }
      printf("\n");
  }



  free(ptr_func_arr);
  free(hamiltonian_mem);
  free(overlaps_mem);
  free(coeffs_mem);

  free(hamiltonians);
  free(overlaps);
  free(coeffs);
  printf("\n");

  return EXIT_SUCCESS;
}
