#include <complex>
#include <fftw3.h>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fftw_interface.hpp"
#define MY_NULLPTR nullptr




///////////////////////////////////////////////////////////////////////////////////////////////
// MY TYPES:
///////////////////////////////////////////////////////////////////////////////////////////////
typedef std::complex<double> Cplx;
typedef double Data;

typedef Cplx (*funcPtr)  (Data* coordinates, double* params, unsigned int ncoords, unsigned int nparams);
typedef Cplx (*integPtr) (Cplx* conjPsi, Cplx* psi);

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// CONFIGURATION:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// DEFINED: 
#define NX 16   // Mesh
#define DX 0.1f // sta³a siaci
#define M  100  // liczba osobników
#define N  16   // liczba f. bazy
#define D  1    // liczba wymiarów

///////////////////////////////////
// D - DIMENSIONAL GRID:
int  g_NQ[D] = { NX };               // User input.
Data g_DQ[D] = { DX };               // User input.
Data g_X0[D] = { 0 };                // Output: run set_grid() function to calculate !!!

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION DEFINITIONS:
///////////////////////////////////////////////////////////////////////////////////////////////

void set_grid() {
#if (_GRID_FIX_)
#	define FIX_Q0 - 1
#else
#	define FIX_Q0
#endif

	for (int i = 0; i < D; i++) {
		g_X0[i] = (Data)((g_NQ[i] FIX_Q0) * g_DQ[i] / 2);
	}

#	undef FIX_Q0
}

Cplx psi(Data* coordinates, double* params, unsigned int ncoords, unsigned int nparams) {
	Cplx psi;
	
	psi.real(1.0f); // ??? real nie zwraca wartosci przez referencje, czy mi siê tylko wydaje ?? ... 
	psi.imag(1.0f); // ??? imag -||-

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
			Cplx conjPsi_j = std::conj( gen_psi(coords, params[i][j], D, param_size) );
			for (int k = 0; k < N; k++) {
				Cplx psi_k = gen_psi(coords, params[i][k], D, param_size);
				S[j][k] = 0.0f;
				H[j][k] = 0.0f;
			}
		}
	}
}

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// START:
///////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {
	
	////////////////////////////////////
	// MEM:
	Cplx** H = new Cplx*[M];
	Cplx** S = new Cplx*[M];
	for (int i = 0; i < M; i++) {
		H[i] = new Cplx[N];
		S[i] = new Cplx[N];
	}

	////////////////////////////////////


	////////////////////////////////////
	// CLEAN-UP:
	for (int i = 0; i < M; i++) {
		delete[] H[i];
		delete[] S[i];
	}
	delete[] H;
	delete[] S;

	////////////////////////////////////
	// EXIT:
	getchar();
	return 0;
}

///////////