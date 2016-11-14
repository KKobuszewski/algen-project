#include <iostream>
#include <complex>
#include <fftw3.h>
#include <cmath>

#include "fftw_interface.hpp"
#define MY_NULLPTR nullptr

/////////////////////////////////////////////////////////////////////////
#define NX 16
#define DX 0.1f
#define M  100
#define N  16
#define D  1

/////////////////////////////////////////////////////////////////////////
typedef std::complex<double> Cplx;
typedef Cplx (*funcPtr)  (double* coordinates, double* params, unsigned int ncoords, unsigned int nparams);
typedef Cplx (*integPtr) (Cplx* cpsi, Cplx* psi);

/////////////////////////////////////////////////////////////////////////
Cplx psi(double* coordinates, double* params, unsigned int ncoords, unsigned int nparams) {
	Cplx psi = { 1.0f, 1.0f };
	return psi;
}

/////////////////////////////////////////////////////////////////////////
void eavluate_hamiltonian(funcPtr gen_psi, Cplx** H, Cplx** S, double*** params, const int param_size) {
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {

			///////////////////////////////////
			double coords[D];
			for (int d = 0; d < D; d++) {
				coords[d] = 0.0F;
			}

			///////////////////////////////////
			Cplx cpsi_j = std::conj( gen_psi(coords, params[i][j], D, param_size) );
			for (int k = 0; k < N; k++) {
				Cplx psi_k = gen_psi(coords, params[i][k], D, param_size);
				S[j][k] = 0.0f;
				H[j][k] = 0.0f;
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
	Cplx** H = new Cplx*[N];
	Cplx** S = new Cplx*[N];
	for (int i = 0; i < M; i++) {
		H[i] = new Cplx[N];
		S[i] = new Cplx[N];
	}

	///////////////////////////////////
	

	getchar();
}