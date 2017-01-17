#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuComplex.h>

#include <omp.h>

#include "Eigensolver.hpp"
#include "cuerrors.h"



Eigensolver::Eigensolver(unsigned matrix_size, unsigned number_matrices, unsigned number_streams):
m(matrix_size),lda(matrix_size), _number_matrices(number_matrices), _number_streams(number_streams)
{
    //cudaError_t cudaStat = cudaSuccess;
    d_A = NULL;
    d_B = NULL;
    d_W = NULL;
    d_work = NULL;
    devInfo = NULL;
    
    printf("Initializing cusolver for matrices %dx%d and %u streams\n",m,lda,number_streams);
    
    
    // allocate memory for matrices
    CUDA_CALL(  cudaMalloc( (void**) &d_A, sizeof(cuDoubleComplex) *lda*m * _number_streams )  );
    CUDA_CALL(  cudaMalloc( (void**) &d_B, sizeof(cuDoubleComplex) *lda*m * _number_streams )  );
    CUDA_CALL(  cudaMalloc( (void**) &d_W, sizeof(double) * m * _number_streams )              );
    CUDA_CALL(  cudaMalloc( (void**) &devInfo, sizeof(int) * _number_matrices )                );
    printf("Memory allocated.\n");
    
    // create multiple cuda streams and create handles
    streams = (cudaStream_t*) malloc( _number_streams * sizeof(cudaStream_t) );
    cusolvers = (cusolverDnHandle_t*) malloc( _number_streams * sizeof(cusolverDnHandle_t) );
    for (unsigned ii=0; ii < _number_streams; ii++)
    {
        cudaStreamCreate(&streams[ii]);
        cusolvers[ii] = NULL;
        
        assert( CUSOLVER_STATUS_SUCCESS == cusolverDnCreate(&cusolvers[ii]) );
        assert( CUSOLVER_STATUS_SUCCESS == cusolverDnSetStream(cusolvers[ii],streams[ii]) );
    }
    printf("Initializing streams..\n");
    
    // step 3: query working space of sygvd
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1; // A*x = (lambda)*B*x
    const cusolverEigMode_t jobz  = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
    const cublasFillMode_t  uplo  = CUBLAS_FILL_MODE_LOWER;
    lwork = 0;
    CUSOLVER_CALL(  cusolverDnZhegvd_bufferSize(cusolvers[0], itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, &lwork)  );
    CUDA_CALL(  cudaMalloc( (void**) &d_work, sizeof(cuDoubleComplex) * lwork * _number_streams )  );
    
    printf("Cusolver prepared.\n");
    printf("lwork: %d\n",lwork);
    printf("\n");
}

Eigensolver::~Eigensolver()
{
    // destroy cusolverfor (unsigned ii=0; ii < _number_streams; ii++)
    for (unsigned ii=0; ii < _number_streams; ii++)
    {
        assert( CUSOLVER_STATUS_SUCCESS == cusolverDnDestroy(cusolvers[ii]) );
    }
    
    // deallocate memory
    if (d_A    ) cudaFree(d_A);
    if (d_B    ) cudaFree(d_B);
    if (d_W    ) cudaFree(d_W);
    if (devInfo) cudaFree(devInfo);
    if (d_work ) cudaFree(d_work);
    if (cusolvers) free(cusolvers);
    if (streams) free(streams);
}

/*
 * Solves generalized eigenvalue problem for many matrices A_i v_i = \lambda_i B_i v_i.
 * 
 * A - hamiltonians
 * B - overlaps
 * V - 
 * W - 
 * 
 */
void Eigensolver::find_generalized_eigvals_batched(
    std::complex<double> **A, std::complex<double> **B,
    std::complex<double> **V, double **E)
{
    const unsigned batch_size = lda*m;
    const cusolverEigType_t itype = CUSOLVER_EIG_TYPE_1;       // A*x = (lambda)*B*x
    const cusolverEigMode_t jobz  = CUSOLVER_EIG_MODE_VECTOR;  // compute eigenvalues and eigenvectors.
    const cublasFillMode_t  uplo  = CUBLAS_FILL_MODE_LOWER;
    
    if (_number_matrices == 1)
    {
        fprintf(stderr,"Case of single dense matrix not implemented!\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        #pragma omp parallel for num_threads(_number_matrices)
        for (unsigned ii=0; ii < _number_matrices; ii++)
        {
            int info_gpu = 0;
            //cudaError_t cudaStat = cudaSuccess;
            //cusolverStatus_t cusolver_status;
            
            CUDA_CALL(  cudaMemcpy(d_A + batch_size*ii, (cuDoubleComplex*) A[ii],
                                  sizeof(cuDoubleComplex) * batch_size, cudaMemcpyHostToDevice)  );
            CUDA_CALL(  cudaMemcpy(d_B + batch_size*ii, (cuDoubleComplex*) B[ii],
                                  sizeof(cuDoubleComplex) * batch_size, cudaMemcpyHostToDevice)  );
            
            // step 4: compute spectrum of (A,B)
            CUSOLVER_CALL(  cusolverDnZhegvd(cusolvers[ii], itype, jobz, uplo,
                            m, d_A + batch_size*ii, lda, d_B + batch_size*ii, lda, d_W + m*ii, d_work + lwork*ii, lwork, &devInfo[ii])  );
            CUDA_CALL(  cudaDeviceSynchronize()  );
            
            //check the result
            CUDA_CALL(  cudaMemcpy(E[ii], d_W + m*ii, sizeof(double)*m, cudaMemcpyDeviceToHost)  );
            CUDA_CALL(  cudaMemcpy(V[ii], d_A + batch_size*ii, sizeof(cuDoubleComplex)*lda*m, cudaMemcpyDeviceToHost)  );
            CUDA_CALL(  cudaMemcpy(&info_gpu, &devInfo[ii], sizeof(int), cudaMemcpyDeviceToHost)  );
            
            if (info_gpu != 0)
                printf("%u. after hegvd: info_gpu = %d\n", ii, info_gpu);
            assert(0 == info_gpu);
        }
        
//         cusolverSpZcsreigs[Host](cusolverH, m, nnz, cusparseMatDescr_t descrA,
//         const cuDoubleComplex *csrValA,
//         const int *csrRowPtrA,
//         const int *csrColIndA,
//         cuDoubleComplex left_bottom_corner,
//         cuDoubleComplex right_upper_corner,
//         int *num_eigs);
    }
}