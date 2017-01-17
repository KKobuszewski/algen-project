#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
// #include <omp.h>
#include <assert.h>

#include <algorithm>
#include <random>

#include <cuda_runtime.h>
#include <cusolverDn.h>

#include <eigensolver/Eigensolver.hpp>


inline void printMatrix(int m, int n, const std::complex<double>*A, int lda, const char* name)
{
    printf("%s:\n",name);
    for (int row = 0 ; row < m ; row++)
    {
        for (int col = 0 ; col < n ; col++)
        {
            double Are = A[row + col*lda].real();
            double Aim = A[row + col*lda].imag();
            if      (Are >= 0 && Aim >= 0)   printf(" %.3e+%.3e ", fabs(Are), fabs(Aim));
            else if (Are >= 0 && Aim <  0)   printf(" %.3e-%.3e ", fabs(Are), fabs(Aim));
            else if (Are <  0 && Aim >= 0)   printf("-%.3e+%.3e ", fabs(Are), fabs(Aim));
            else                             printf("-%.3e-%.3e ", fabs(Are), fabs(Aim));
        }
        printf("\n");
    }
}

inline void printMatrix(int m, int n, const double*A, int lda, const char* name)
{
    printf("%s:\n",name);
    for (int row = 0 ; row < m ; row++)
    {
        for (int col = 0 ; col < n ; col++)
        {
            double Areg = A[row + col*lda];
            if (Areg >= 0)   printf(" %.3e ", Areg);
            else            printf( "%.3e ", Areg);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[])
{
    const unsigned m = 8; // matrix size
    const unsigned n = 64; // number matrices
    
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937_64 mersenne(seed);
    std::uniform_real_distribution<double> uniform_dist(1., 2.);
    std::normal_distribution<double> gaussian_dist(0.,0.1);
    
    
    std::complex<double> *H[n];
    std::complex<double> *S[n];
    std::complex<double> *V[n];
    double *E[n];
    
    for (unsigned ii=0; ii<n; ii++)
    {
        H[ii] = (std::complex<double>*) malloc( sizeof(std::complex<double>) * m*m );
        S[ii] = (std::complex<double>*) malloc( sizeof(std::complex<double>) * m*m );
        V[ii] = (std::complex<double>*) malloc( sizeof(std::complex<double>) * m*m );
        E[ii] = (double*) malloc( sizeof(double) * m );
    }
    
    // construct random hermitean matrices
    for (unsigned ii=0; ii<n; ii++)
    {
        for (unsigned jj=0; jj<m; jj++)
        {
            //for (unsigned kk=1; kk<jj; kk++)
            for (unsigned kk=0; kk<m; kk++)
            {
                if (kk == jj)
                {
                    std::complex<double> z(5*(jj+1.)*uniform_dist(mersenne),0);
                    H[ii][jj*m+kk] = z;
                }
                else
                {
                    std::complex<double> z(gaussian_dist(mersenne),gaussian_dist(mersenne));
                    H[ii][jj*m+kk] = z;
                    H[ii][kk*m+jj] = std::conj(z);
                }
            }
        }
    }
    for (unsigned ii=0; ii<n; ii++)
    {
        for (unsigned jj=0; jj<m; jj++)
        {
            for (unsigned kk=0; kk<m; kk++)
            {
                if (kk == jj)
                {
                    std::complex<double> z(1.,0);
                    S[ii][jj*m+kk] = z;
                }
                else
                {
                    //std::complex<double> z(gaussian_dist(mersenne),gaussian_dist(mersenne));
                    std::complex<double> z(0.,0.);
                    S[ii][jj*m+kk] = z;
                    S[ii][kk*m+jj] = std::conj(z);
                }
            }
        }
    }
    
    
    for (unsigned ii=0; ii<n; ii++)
    {
        printMatrix(m,m,H[ii],m,"H");
        printMatrix(m,m,S[ii],m,"S");
        printf("\n");
    }
    
    
    Eigensolver eig(m,n,n);
    eig.find_generalized_eigvals_batched(H,S,V,E);
    
    
    for (unsigned ii=0; ii<n; ii++)
    {
        printMatrix(m,1,E[ii],1,"eigenvalues");
        printMatrix(m,m,V[ii],m,"eigenvectors");
        printf("\n");
    }
    
    for (unsigned ii=0; ii<n; ii++)
    {
        free(H[ii]);
        free(S[ii]);
        free(V[ii]);
        free(E[ii]);
    }
    
    return EXIT_SUCCESS;
}