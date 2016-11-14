#include <complex>
#include <fftw3.h>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fftw_interface.hpp"


///////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////
#define MY_NULLPTR nullptr

///////////






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
// PROGRAM MODES:
#define _GRID_FIX_  0              // Bad for fourier transform. Keep as false.

///////////////////////////////////
// DEFINED: 
#define NX 256                     // Mesh
#define DX 0.1f                    // sta³a siaci

#define M  100                     // liczba osobników
#define N  16                      // liczba f. bazy
#define D  1                       // liczba wymiarów

///////////////////////////////////
// D - DIMENSIONAL GRID:
int  g_nq[D] = { 0 };              // User input.
Data g_dq[D] = { 0 };              // User input.
Data g_q0[D] = { 0 };              // Output: run set_grid() function to calculate !!!

///////////






///////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION DEFINITIONS:
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////
// D - DIMENSIONAL GRID:
void set_grid() {
#if (_GRID_FIX_)
#	define FIX_Q0 - 1
#else
#	define FIX_Q0
#endif

	for (int i = 0; i < D; i++) {
		g_q0[i] = (Data)((g_nq[i] FIX_Q0) * g_dq[i] / 2);
	}

#	undef FIX_Q0
}

///////////////////////////////////
// QM RELATED:
Cplx psi(Data* coordinates, double* params, unsigned int ncoords, unsigned int nparams) {
	Cplx psi;
	
	psi.real(1.0f); // ??? real nie zwraca wartosci przez referencje, czy tylko mi siê wydaje ?? ... 
	psi.imag(1.0f); // ??? imag -||-

	return psi;
}
void eavluate_hamiltonian(funcPtr gen_psi, Cplx** H, Cplx** S, double*** params, const int param_size) {
	///////////////////////////////////
	// GEN COORDS:
	Data** coords_all = new Data*[D];
	for (int d = 0; d < D; d++) {
		// ALLOC MEM: 
		coords_all[d] = new Data[g_nq[d]];

		// SET VALUES:
		for (int i = 0; i < g_nq[d]; i++) {
			coords_all[d][i] = g_q0[d] + i * g_dq[d];
		}
	}

	///////////////////////////////////
	// DO YOUR WORST:
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {

			Data coords[D];
			for (int d = 0; d < D; d++) {
				for (int mesh = 0; mesh < g_nq[d]; mesh++) {
					
					///////////////////////////////////
					// GET RIGHT COORDS:
					for (int dim = 0; dim < D; dim++) {
						coords[dim] = coords_all[dim][mesh];
					}

					///////////////////////////////////
					// DO YOUR WORST:
					Cplx conjPsi_j = std::conj(gen_psi(coords, params[i][j], D, param_size));
					for (int k = 0; k < N; k++) {
						Cplx psi_k = gen_psi(coords, params[i][k], D, param_size);
						S[j][k] += conjPsi_j * psi_k;  // sumowanie s[j][k] jako wstêp do ca³kowania. // ca³kowanie i tak sprowadzi sie do sumy z tablicy values[ g_nq[0] + g_nq[1] + g_nq[2] + ... + g_nq[D-1] ];
						H[j][k] += 0.0f;
					}
				}
			}
		}
	}

	////////////////////////////////////
	// CLEAN-UP:
	for (int i = 0; i < D; i++) {
		delete[] coords_all[i];
	}
	delete[] coords_all;
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
	// WORK-SPACE:
	g_nq[0] = NX; // dowolny mechanizm wype³nienia sieci. Taki wygodny na pocz¹tek.
	g_dq[0] = DX; // -||-
	set_grid();   // Wygeneruj punkty sieci Q0;





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