/***************************************************************************
 *   Copyright (C) 2015 by                                                 *
 *   WARSAW UNIVERSITY OF TECHNOLOGY                                       *
 *   FACULTY OF PHYSICS                                                    *
 *   NUCLEAR THEORY GROUP                                                  *
 *   AUTHORS: Konrad Kobuszewski, Gabriel Wlazłowski                       *
 *                                                                         *
 *   This file is a part of GPE for GPU project.                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef CUERRORS_H_
#define CUERRORS_H_



#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include <cublas_v2.h>
#include <curand_kernel.h>





/* *********************************************************************************************************** *
 *                                                                                                             *
 *                                            CUDA ERRORS HANDLING                                             *
 *                                                                                                             *
 * *********************************************************************************************************** */


// =================================== CUDA ERROR HANDLING ===================================================

/*
 * This macro enables simple handling of cudaError_t, and passes error as gpe_result_t to gpe_exec macro
 * TODO: if could be in gpe_engine.cuh?
 */                                                                            \

inline void _cuErrCheck(const cudaError_t err, const char *file, const int line)
{
  if(err != cudaSuccess)
      {
          fprintf( stderr, "\n");
          fprintf( stderr, "CUDA ERROR in file=`%s`, line=%d\n",file,line);
          fprintf( stderr, "CUDA ERROR %d: %s\n", err, cudaGetErrorString((cudaError_t)(err)) );
          fprintf( stderr, "\n");
          cudaDeviceReset();
          exit(EXIT_FAILURE);
      }
}

#define CUDA_CALL(err)         _cuErrCheck(err, __FILE__, __LINE__)


// =================================== CUFFT ERROR HANDLING ==================================================

static inline const char *cufftGetErrorString(cufftResult error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            return "CUFFT_INCOMPLETE_PARAMETER_LIST";
        case CUFFT_INVALID_DEVICE:
            return "CUFFT_INVALID_DEVICE";
        case CUFFT_NO_WORKSPACE:
            return "CUFFT_NO_WORKSPACE";
        case CUFFT_PARSE_ERROR:
            return "CUFFT_PARSE_ERROR";
        case CUFFT_NOT_IMPLEMENTED:
            return "CUFFT_NOT_IMPLEMENTED";
        case CUFFT_LICENSE_ERROR:
          return "CUFFT_LICENSE_ERROR";
        default:
          return "CUFFT UNKNOWN ERROR!";
    }
}

inline void _cufftErrChk(cufftResult err, const char *file, const int line)
{
    if( CUFFT_SUCCESS != err) {
        fprintf( stderr, "\n");
        fprintf( stderr, "CUFFT ERROR in file=`%s`, line=%d\n",file,line);
        fprintf( stderr, "CUFFT ERROR %d: %s\n", err, cufftGetErrorString((cufftResult)(err)) );
        fprintf( stderr, "\n");
        cudaDeviceReset();
        exit(EXIT_FAILURE);
    }
}

#define CUFFT_CALL(err)      _cufftErrChk(err, __FILE__, __LINE__)


// =================================== CUBLAS ERROR HANDLING =================================================

static inline const char* cublasGetErrorString(cublasStatus_t error) {
  switch (error)
  {
  case CUBLAS_STATUS_SUCCESS:
    return "CUBLAS_STATUS_SUCCESS";
  case CUBLAS_STATUS_NOT_INITIALIZED:
    return "CUBLAS_STATUS_NOT_INITIALIZED";
  case CUBLAS_STATUS_ALLOC_FAILED:
    return "CUBLAS_STATUS_ALLOC_FAILED";
  case CUBLAS_STATUS_INVALID_VALUE:
    return "CUBLAS_STATUS_INVALID_VALUE";
  case CUBLAS_STATUS_ARCH_MISMATCH:
    return "CUBLAS_STATUS_ARCH_MISMATCH";
  case CUBLAS_STATUS_MAPPING_ERROR:
    return "CUBLAS_STATUS_MAPPING_ERROR";
  case CUBLAS_STATUS_EXECUTION_FAILED:
    return "CUBLAS_STATUS_EXECUTION_FAILED";
  case CUBLAS_STATUS_INTERNAL_ERROR:
    return "CUBLAS_STATUS_INTERNAL_ERROR";
  case CUBLAS_STATUS_NOT_SUPPORTED:
    return "CUBLAS_STATUS_NOT_SUPPORTED";
  case CUBLAS_STATUS_LICENSE_ERROR:
    return "CUBLAS_STATUS_LICENSE_ERROR";
  }
  return "Unknown cublas status";
}

inline void _cublasErrChk(cublasStatus_t err, const char *file, const int line)
{
    if( CUBLAS_STATUS_SUCCESS != err) {
        fprintf( stderr, "\n");
        fprintf( stderr, "CURAND ERROR in file=`%s`, line=%d\n",file,line);
        fprintf( stderr, "CURAND ERROR %d: %s\n", err, cublasGetErrorString((cublasStatus_t)(err)) );
        fprintf( stderr, "\n");
        cudaDeviceReset();
        exit(EXIT_FAILURE);
    }
}

#define CUBLAS_CALL(err)    _cublasErrChk(err, __FILE__, __LINE__)


// =================================== CURAND ERROR HANDLING =================================================

static inline const char* curandGetErrorString(curandStatus_t error) {
  switch (error)
  {
  case CURAND_STATUS_SUCCESS:
    return "CURAND_STATUS_SUCCESS";
  case CURAND_STATUS_VERSION_MISMATCH:
    return "CURAND_STATUS_VERSION_MISMATCH";
  case CURAND_STATUS_NOT_INITIALIZED:
    return "CURAND_STATUS_NOT_INITIALIZED";
  case CURAND_STATUS_ALLOCATION_FAILED:
    return "CURAND_STATUS_ALLOCATION_FAILED";
  case CURAND_STATUS_TYPE_ERROR:
    return "CURAND_STATUS_TYPE_ERROR";
  case CURAND_STATUS_OUT_OF_RANGE:
    return "CURAND_STATUS_OUT_OF_RANGE";
  case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
    return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
  case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
    return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
  case CURAND_STATUS_LAUNCH_FAILURE:
    return "CURAND_STATUS_LAUNCH_FAILURE";
  case CURAND_STATUS_PREEXISTING_FAILURE:
    return "CURAND_STATUS_PREEXISTING_FAILURE";
  case CURAND_STATUS_INITIALIZATION_FAILED:
    return "CURAND_STATUS_INITIALIZATION_FAILED";
  case CURAND_STATUS_ARCH_MISMATCH:
    return "CURAND_STATUS_ARCH_MISMATCH";
  case CURAND_STATUS_INTERNAL_ERROR:
    return "CURAND_STATUS_INTERNAL_ERROR";
  }
  return "Unknown curand status";
}

inline void _curandErrChk(curandStatus_t err, const char *file, const int line)
{
    if( CURAND_STATUS_SUCCESS != err) {
        fprintf( stderr, "\n");
        fprintf( stderr, "CURAND ERROR in file=`%s`, line=%d\n",file,line);
        fprintf( stderr, "CURAND ERROR %d: %s\n", err, curandGetErrorString((curandStatus_t)(err)) );
        fprintf( stderr, "\n");
        cudaDeviceReset();
        exit(EXIT_FAILURE);
    }
}

#define CURAND_CALL(err)    _curandErrChk(err, __FILE__, __LINE__)


// =================================== CUSOLVER ERROR HANDLING =================================================

static inline const char* cusolverGetErrorString(curandStatus_t error) {
  switch (error)
  {
  case CUSOLVER_STATUS_SUCCESS:
    return "CUSOLVER_STATUS_SUCCESS";
  }
  return "Unknown curand status";
}

inline void _cusolverErrChk(cusolverStatus_t err, const char *file, const int line)
{
    if( CUSOLVER_STATUS_SUCCESS != err) {
        fprintf( stderr, "\n");
        fprintf( stderr, "CUSOLVER ERROR in file=`%s`, line=%d\n",file,line);
        fprintf( stderr, "CUSOLVER ERROR %d: %s\n", err, cusolverGetErrorString((curandStatus_t)(err)) );
        fprintf( stderr, "\n");
        cudaDeviceReset();
        exit(EXIT_FAILURE);
    }
}

#define CUSOLVER_CALL(err)    _cusolverErrChk(err, __FILE__, __LINE__)


// =================================== CUSPARSE ERROR HANDLING =================================================

#endif /* CUERRORS_H_ */
