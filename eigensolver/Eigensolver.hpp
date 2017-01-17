#ifndef __EIGENSOLVER_HPP__
#define __EIGENSOLVER_HPP__

#include <stdlib.h>
#include <assert.h>
#include <error.h>

#include <complex>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuComplex.h>


class Eigensolver
{
private:
    cusolverDnHandle_t* cusolvers;
//     cusparseHandle_t cusparseH;
//     cusparseStatus_t cusparse_status;
    const int m, lda;                                             // eigenvectors/eigenvalues
    cuDoubleComplex *d_A, *d_B, *d_work;
    double *d_W;
    int *devInfo;
    int lwork;
    unsigned _number_matrices;
    unsigned _number_streams;
    cudaStream_t* streams;
public:
    // consturctors
    Eigensolver(unsigned matrix_size, 
                unsigned number_matrices, 
                unsigned number_streams);
    
    // destructors
    ~Eigensolver();
    
    // methods
    void find_generalized_eigvals_batched(std::complex<double> **A, 
                                          std::complex<double> **B,
                                          std::complex<double> **V,
                                          double **E);
    
};


#endif