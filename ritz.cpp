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
#include <complex.h>
#include <math.h>


/*
 * Basis functions.
 * They take arguments: double* point, double* arguments, ndims, nparams
 * They return value of type double __complex__
 */
typedef double __complex__ (*basis_ptr)(double*, double*, unsigned, unsigned);


int main()
{
  const unsigned N = 8;                  // basis size
  const unsigned m = 100;                // population size
  const unsigned d = 1;                  // number of dimensions
  const unsigned npart = 1;              // number of particles

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





  free(ptr_func_arr);
  free(hamiltonian_mem);
  free(overlaps_mem);
  free(coeffs_mem);

  free(hamiltonians);
  free(overlaps);
  free(coeffs);

  return EXIT_SUCCESS;
}
