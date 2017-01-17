#ifndef __SIMULATIONMANAGER_HPP_
#define __SIMULATIONMANAGER_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include <map>
#include <algorithm>
#include <complex>

#include <random/Random.hpp>
#include <basis.h>
#include <hamiltonian.h>
#include <integrator/Integrator.hpp>

#ifndef THREADS
#define THREADS 8
#endif

#ifndef MAX_POPULATION
#define MAX_POPULATION 1024
#endif

#ifndef MAX_NBASIS
#define MAX_NBASIS 1024
#endif


namespace SimulationLimits
{

const unsigned max_population = MAX_POPULATION;
const unsigned max_nbasis     = MAX_NBASIS;
const unsigned max_threads    = THREADS;

}


class SimulationManager
{
private:
    Random* rng;  // LICZBY LOSOWE
    Integrator** integs;                                                     // CAŁKOWANIE <- do Ritza
    ptr_func_t*  basis_assortment;  // zbiór wszystkich funkcji wybieranych do bazy (128)
    unsigned*    shuffle_assortment; // tablica liczb {0,1,...,128-1} <- indeksy funkcji z basis_assortment // gdzies tam dalej jest permutacja liczb z tej tablicy, do inicjalizacji funkcji bazy dla określonego osobnika
    // flags
    bool         indvs_copied;
public:
    ptr_func_t*  basis_mem;      // pamięć na funkcje bazy osobników, o wielkości 8 x 32(liczba osobnikow)
    unsigned*    population_mem; // pamięć na indeksy {0,1,...,128-1}, o wielkości 8 x 32(liczba osobnikow)
    double*      params_mem;
    std::complex<double> *hamiltonians_mem, *overlaps_mem, *eigvecs_mem; // pamięć na macierze H, o wiekości (8x8) x 32(liczba osobników)
    double*      energies_mem;   // pamięć na energie, o wielkości 8 x 32(liczba osobnikow)
    double*      adaptations_cpy; // pamięć na fukcje celu, o wielkości 1 x 32(liczba osobnikow)
    double*      energies_mem_cpy; // kopie
    ptr_func_t*  basis_mem_cpy;
    unsigned*    population_mem_cpy; 
    double*      params_mem_cpy;
    // accessors
    unsigned dims, tot_dims, npart, nbasis, mpopulation, nparams, size_assortment;
    std::complex<double> **hamiltonians, **overlaps, **eigvecs;      // hamiltonians[ii] to jest macierz H ii-tego osobnika
    ptr_func_t** basis;  // basis[ii] to wektor funkcji bazy ii-tego osobnika, basis[ii][jj] to jj-ta funkcja bazy ii-tego osobnika
    // basis[ii=0..32-1][jj=0..8-1]    ==    basis_mem[  ( ii = {0...32-1} ) * 8 + (jj = {0...8-1})  ]
    unsigned**   population;  // to samo co wyżej tylko na indeksach funkcji w basis_assortment
    double**     params_to_funcs;  // parametry funkcji bazowej
    double**     energies;
    double*      adaptations;
    
    // constructors
    SimulationManager(unsigned _d, unsigned _npart, unsigned _nbasis, unsigned _mpopulation, unsigned _nparams,
                      ptr_func_t* _basis_assortment, unsigned _size_assortment);
    
    // destructors
    ~SimulationManager();
    
    
    // methods
    void copy_individuals();
    unsigned* get_indiv_func_indices(unsigned individual_index);
    double*   get_indiv_params(unsigned individual_index);
    std::complex<double>* get_hamiltonian(unsigned individual_index);   // oblicza macierz H dla individual_index osobnika
    std::complex<double>* get_overlap_mat(unsigned individual_index);   // oblicza macierz S (przekryć) dla individual_index osobnika
    
    void evaluate_matrices(unsigned samples, double* params, double* xl, double* xu); // stosuje dwie powyższe funkcje do każdego osobnika z osobna
    
    void print_matrices(unsigned specimen);
    void print_population();
    void print_functions();
    void print_params();
    std::map<unsigned, unsigned> histogramize_basis_functions();
    void print_basis_function_histogram();
    
    
    // for genetic algorithms
    void reproduce_individual(unsigned parent_index, unsigned child_index);
    void swap_basis_func(unsigned indv1_idx, unsigned indv2_idx, unsigned func1_idx, unsigned func2_idx);
    void mutate_basis_func(unsigned individual_index, unsigned basis_index, unsigned new_wavefunction_index);
    void check_basis_func_different();
};

/* ***************************************************************************************************************** *
 *                                                                                                                   *
 * Constructor of class SimulationManager. Allocates memory for simulation and initialize with random data.          *
 *                                                                                                                   *
 * @param unsigned    _d                 - number of spatial dimensions                                              *
 * @param unsigned    _npart             - number of particles in system                                             *
 * @param unsigned    _nbasis            - number of basis wavefunctions used for computation                        *
 * @param unsigned    _mpopulation       - number of individuals                                                     *
 * @param unsigned    _nparams           - number of parameters being passed to wavefunction                         *
 * @param ptr_func_t* _basis_assortment  - array containing all possible bassis functions                            *
 * @param unsigned    _size_assortment   - number of all possible basis functions                                    *
 *                                                                                                                   *
 * ***************************************************************************************************************** */
SimulationManager::SimulationManager(
    unsigned _d, unsigned _npart, unsigned _nbasis, unsigned _mpopulation, unsigned _nparams,
    ptr_func_t* _basis_assortment, unsigned _size_assortment) : 
    dims(_d), npart(_npart), nbasis(_nbasis), mpopulation(_mpopulation), nparams(_nparams),
    basis_assortment(_basis_assortment), size_assortment(_size_assortment)
{
    // =================  prerequests  ================================================================================
    tot_dims = _d*_npart;
    rng = new Random();
    integs = new Integrator*[THREADS];
    for (unsigned ii=0; ii<THREADS; ii++) integs[ii] = new Integrator(); // different integrator for each OMP thread
    indvs_copied = false;
    
    
    // =================  allocate memory  ============================================================================
    // vectors for each individual
    population       = (unsigned**)   malloc( mpopulation * sizeof(unsigned*) );
    population_mem   = (unsigned*)    malloc( mpopulation * nbasis * sizeof(unsigned) );
    params_to_funcs  = (double**)     malloc( mpopulation * sizeof(double*) );
    params_mem       = (double*)      malloc( mpopulation * nbasis*nparams * sizeof(double) );
    basis            = (ptr_func_t**) malloc( mpopulation * sizeof(ptr_func_t*) );
    basis_mem        = (ptr_func_t*)  malloc( mpopulation * nbasis * sizeof(ptr_func_t) );
    adaptations      = (double*)      malloc( mpopulation * sizeof(double) );
    
    // matrices for each individual
    hamiltonians     = (std::complex<double>**) malloc( mpopulation * sizeof(std::complex<double>*) );
    hamiltonians_mem = (std::complex<double>*)  malloc( mpopulation * nbasis*nbasis * sizeof(std::complex<double>) );
    overlaps         = (std::complex<double>**) malloc( mpopulation * sizeof(std::complex<double>*) );
    overlaps_mem     = (std::complex<double>*)  malloc( mpopulation * nbasis*nbasis * sizeof(std::complex<double>) );
    eigvecs          = (std::complex<double>**) malloc( mpopulation * sizeof(std::complex<double>*) );
    eigvecs_mem      = (std::complex<double>*)  malloc( mpopulation * nbasis*nbasis * sizeof(std::complex<double>) );
    energies         = (double**) malloc( mpopulation * sizeof(double*) );
    energies_mem     = (double*)  malloc( mpopulation * nbasis*nbasis * sizeof(double) );
    
    // allocate memory for copies
    population_mem_cpy   = (unsigned*)    malloc( mpopulation * nbasis * sizeof(unsigned) );
    params_mem_cpy       = (double*)      malloc( mpopulation * nbasis*nparams * sizeof(double) );
    basis_mem_cpy        = (ptr_func_t*)  malloc( mpopulation * nbasis * sizeof(ptr_func_t) );
    energies_mem_cpy     = (double*)      malloc( mpopulation * nbasis * sizeof(double) );
    adaptations_cpy      = (double*)      malloc( mpopulation * sizeof(double) );
    
    
    // 
    if (size_assortment > nbasis)
    {
        shuffle_assortment = (unsigned*) malloc( size_assortment*sizeof(unsigned) );
        for (unsigned ii=0; ii<size_assortment; ii++)
            shuffle_assortment[ii] = ii;
    }
    else
        shuffle_assortment = NULL;
    
    
    // =================  fill in arrays of ptrs  =====================================================================
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        population[kk] = population_mem + kk*nbasis;
    }
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        params_to_funcs[kk] = params_mem + kk*nbasis;
    }
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        basis[kk] = basis_mem + kk*nbasis;    // nbasis to 8
    }
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        hamiltonians[kk] = hamiltonians_mem + kk*nbasis*nbasis;
    }
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        overlaps[kk] = overlaps_mem + kk*nbasis*nbasis;
    }
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        eigvecs[kk] = eigvecs_mem + kk*nbasis*nbasis;
    }
    #pragma omp simd
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        energies[kk] = energies_mem + kk*nbasis*nbasis;
    }
    
    // =================  initialize population with a random vectors =================================================
    //#pragma omp parallel for num_threads(8)
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        /*#pragma omp simd
        for (unsigned ii=0; ii < nbasis; ii++)
        {
            population_mem[nbasis*kk + ii] = rng->get_random_uniform();
        }*/
        
        
        for (unsigned ii=0; ii < nbasis; ii++)
        {
            #pragma omp simd
            for (unsigned jj=0; jj<nparams; jj++)
//                 params_mem[nbasis*nparams*kk + ii*nparams + jj] = rng->get_random_uniform();
                params_mem[nbasis*nparams*kk + ii*nparams + jj] = 1.;
        }
        
        
        if (size_assortment == nbasis)
        {
            #pragma omp simd
            for (unsigned ii=0; ii < nbasis; ii++)
            {
                basis[kk][ii] = basis_assortment[ii];
            }
            #pragma omp simd
            for (unsigned ii=0; ii < nbasis; ii++)
            {
                population[kk][ii] = ii;
            }
        }
        else if (size_assortment > nbasis)
        {
            // size_assortment allocation before
            rng->shuffle_array<unsigned>(shuffle_assortment,size_assortment);
            std::sort(shuffle_assortment,shuffle_assortment+nbasis);
            
            #pragma omp simd
            for (unsigned ii=0; ii < nbasis; ii++)
            {
                population[kk][ii] = shuffle_assortment[ii];
                //printf("%u/%u. function drawn\n",kk,shuffle_assortment[ii]);
            }
            
            #pragma omp simd
            for (unsigned ii=0; ii < nbasis; ii++)
            {
                basis[kk][ii] = basis_assortment[shuffle_assortment[ii]];
                //printf("%u/%u. function drawn\n",kk,shuffle_assortment[ii]);
            }
        }
        else
        {
            fprintf(stderr,"\nERROR! (SimulationManager): Too less basis functions in a assortment!!!\n");
            exit(EXIT_FAILURE);
        }
    }
}

SimulationManager::~SimulationManager()
{
    delete rng;
    for (unsigned ii=0; ii<THREADS; ii++) delete integs[ii];
    delete[] integs;
    //printf("Integrators deleted.\n");
    
//     if (shuffle_assortment)  free(shuffle_assortment);
//     
//     if (population_mem)      free(population_mem);
//     if (population)          free(population);
//     if (params_mem)          free(params_mem);
//     if (params_to_funcs)     free(params_to_funcs);
//     if (basis_mem)           free(basis_mem);
//     if (basis)               free(basis);
//     if (hamiltonians_mem)    free(hamiltonians_mem);
//     if (hamiltonians)        free(hamiltonians);
//     if (overlaps_mem)        free(overlaps_mem);
//     if (overlaps)            free(overlaps);
    if (shuffle_assortment)  free(shuffle_assortment);
#ifdef DEBUG
    printf("shuffle_assortment OK\n");
#endif
    free(population_mem);
    free(population);
#ifdef DEBUG
    printf("population OK\n");
#endif
    free(params_mem);
    free(params_to_funcs);
#ifdef DEBUG
    printf("params OK\n");
#endif
    free(basis_mem);
    free(basis);
#ifdef DEBUG
    printf("basis OK\n");
#endif
    free(hamiltonians_mem);
    free(hamiltonians);
#ifdef DEBUG
    printf("hamiltonians OK\n");
#endif
    free(overlaps_mem);
    free(overlaps);
#ifdef DEBUG
    printf("overlaps OK\n");
#endif
    free(eigvecs_mem);
    free(eigvecs);
#ifdef DEBUG
    printf("eigvecs OK\n");
#endif
    free(energies_mem);
    free(energies);
#ifdef DEBUG
    printf("energies OK\n");
#endif
    free(adaptations);
#ifdef DEBUG
    printf("adaptations OK\n");
#endif
    
    free(population_mem_cpy);
    free(params_mem_cpy);
    free(basis_mem_cpy);
    free(energies_mem_cpy);
    free(adaptations_cpy);
#ifdef DEBUG
    printf("copies OK\n");
#endif
    
    
    //sleep(5);
    
}

/*
 * Make copy of most important individual's atributes to additional arrays to allow overwriting
 * these atributes. This method should be used before reproducing and recombining individuals.
 * Indices, params and ptrs func of basis functions are copied. Also energies and adaptations.
 * 
 * NOTE: This method is automatically called by SimulationManager::evaluate_matrices!
 */
inline void SimulationManager::copy_individuals()
{
    // copy individuals
    memcpy( population_mem_cpy, population_mem, mpopulation*nbasis         * sizeof(unsigned)   );
    memcpy( params_mem_cpy,     params_mem,     mpopulation*nbasis*nparams * sizeof(double)     );
    memcpy( basis_mem_cpy,      basis_mem,      mpopulation*nbasis         * sizeof(ptr_func_t) );
    memcpy( energies_mem_cpy,   energies_mem,   mpopulation*nbasis         * sizeof(double)     );
    memcpy( adaptations_cpy,    adaptations,    mpopulation                * sizeof(double)     );
    
    indvs_copied = true;
}


inline unsigned* SimulationManager::get_indiv_func_indices(const unsigned individual_index)
{
    return population[individual_index];
}

inline double* SimulationManager::get_indiv_params(const unsigned individual_index)
{
    return params_to_funcs[individual_index];
}

inline std::complex<double>* SimulationManager::get_hamiltonian(const unsigned individual_index)
{
    return hamiltonians[individual_index];
}

inline std::complex<double>* SimulationManager::get_overlap_mat(const unsigned individual_index)
{
    return overlaps[individual_index];
}


inline void SimulationManager::evaluate_matrices(const unsigned samples, double* p, double* xl, double* xu)
{
    #pragma omp parallel for num_threads(THREADS)
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
#ifdef DEBUG
        printf("(omp_thread %u) Evaluating matrices for %u-th individual...\n",omp_get_thread_num(),kk);
#endif
        
        Integrator* integ = integs[omp_get_thread_num()];
        
        #pragma omp simd
        for (unsigned ii=0; ii < nbasis; ii++)
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            double* params = params_to_funcs[kk];
            hamiltonians[kk][ii*nbasis + jj] = integ->integrate(hamiltonian_elem, basis[kk][ii], basis[kk][jj],
                                                                params, dims, nparams, npart, samples, xl, xu);
            overlaps[kk][ii*nbasis + jj]     = integ->integrate(identity,         basis[kk][ii], basis[kk][jj],
                                                                params, dims, nparams, npart, samples, xl, xu);
        }
        
        // make matrices hermitian
        #pragma omp simd
        for (unsigned ii=0;  ii < nbasis; ii++)
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            // for improved accurancy?
            std::complex<double> mean = .5*( hamiltonians[kk][ii*nbasis + jj] + hamiltonians[kk][jj*nbasis + ii] );
            hamiltonians[kk][ii*nbasis + jj] = mean;
            hamiltonians[kk][jj*nbasis + ii] = std::conj(mean);
        }
        #pragma omp simd
        for (unsigned ii=0;  ii < nbasis; ii++)
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            // for improved accurancy?
            std::complex<double> mean = .5*( overlaps[kk][ii*nbasis + jj] + overlaps[kk][jj*nbasis + ii] );
            overlaps[kk][ii*nbasis + jj] = mean;
            overlaps[kk][jj*nbasis + ii] = std::conj(mean);
        }
    }
    
    // copy individuals
    this->copy_individuals();
    //memcpy( population_mem_cpy, population_mem, mpopulation*nbasis         * sizeof(unsigned) );
    //memcpy( params_mem_cpy,     params_mem,     mpopulation*nbasis*nparams * sizeof(double) );
    //memcpy( basis_mem_cpy,      basis_mem,      mpopulation*nbasis         * sizeof(ptr_func_t) );
    //memcpy( energies_mem_cpy,   energies_mem,   mpopulation*nbasis         * sizeof(double) );
    //memcpy( adaptations_cpy,    adaptations,    mpopulation                * sizeof(double) );
}



inline void SimulationManager::print_matrices(const unsigned specimen)
{
    printf("========================= INDIVIDUAL %u =======================================================\n",
           specimen);
    printf("|                                                                                             \n");
    printf("|   Basis function indices:                                                                   \n");
    printf("|  ");
    for (unsigned jj=0; jj < nbasis; jj++)
    {
        printf(" %u ",population[specimen][jj]);
    }
    printf("\n");
    printf("|                                                                                             \n");
    printf("|   Energies:                                                                                 \n");
    printf("|  ");
    for (unsigned jj=0; jj < nbasis; jj++)
    {
        double ei = energies[specimen][jj];
        if (ei >= 0) printf(" ");
        printf("%.3lf ",ei);
    }
    printf("\n");
    printf("|                                                                                             \n");
    printf("|   Eigenvecs:                                                                                \n");
    for (unsigned ii=0; ii < nbasis; ii++)
    {
        printf("|  ");
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            std::complex<double> vij = eigvecs[specimen][ii*nbasis + jj];
            if (vij.real() >= 0) printf(" ");
            if (vij.imag() >= 0) printf("%.1e+%.1ej ",vij.real(),fabs(vij.imag()));
            else                 printf("%.1e+%.1ej ",vij.real(),-1*vij.imag());
        }
        printf("\n");
    }
    printf("|                                                                                             \n");
    printf("|   Hamiltonian:                                                                              \n");
    for (unsigned ii=0; ii < nbasis; ii++)
    {
        printf("|  ");
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            std::complex<double> hij = hamiltonians[specimen][ii*nbasis + jj];
            if (hij.real() >= 0) printf(" ");
            if (hij.imag() >= 0) printf("%.1e+%.1ej ",hij.real(),fabs(hij.imag()));
            else                 printf("%.1e+%.1ej ",hij.real(),-1*hij.imag());
        }
        printf("\n");
    }
    printf("|                                                                                             \n");
    printf("|   Overlaps:                                                                                  \n");
    for (unsigned ii=0; ii < nbasis; ii++)
    {
        printf("|  ");
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            std::complex<double> sij = overlaps[specimen][ii*nbasis + jj];
            if (sij.real() >= 0) printf(" ");
            if (sij.imag() >= 0) printf("%.1e+%.1ej ",sij.real(),fabs(sij.imag()));
            else                 printf("%.1e+%.1ej ",sij.real(),-1*sij.imag());
        }
        printf("\n");
    }
    printf("|                                                                                             \n");
    printf("|_______________________________________________________________________________________________\n");
    printf("\n");
}

inline void SimulationManager::print_population()
{
    for (unsigned kk=0; kk < mpopulation; kk++)
    {
        SimulationManager::print_matrices(kk);
    }
}

// void SimulationManager::print_all_matrices()
// {
//     for (unsigned kk=0; kk < mpopulation; kk++)
//     {
//         SimulationManager::print_matrices(kk);
//     }
// }


inline void SimulationManager::print_functions()
{
    printf("\n# BASIS FUNCTION INDICES OF INDIVIDUALS\n\n");
    for (unsigned kk=0; kk<mpopulation; kk++)
    {
        printf("Individual %u:\t",kk);
        for (unsigned jj=0; jj < nbasis; jj++)
        {
            printf("%u ",population[kk][jj]);
        }
        printf("\n");
    }
}



inline void SimulationManager::print_params()
{
    printf("\n# PARAMETERS OF INDIVIDUALS\n\n");
    for (unsigned kk=0; kk<mpopulation; kk++)
    {
        printf("Individual %u:\t",kk);
        for (unsigned jj=0; jj < nparams*nbasis; jj++)
        {
            printf("%.3lf ",params_to_funcs[kk][jj]);
        }
        printf("\n");
    }
}

inline std::map<unsigned, unsigned> SimulationManager::histogramize_basis_functions()
{
    std::map<unsigned, unsigned> hist;
    for (unsigned n=0; n<mpopulation*nbasis; n++)
    {
        ++hist[population_mem[n]];
    }
    
    return hist;
}

inline void SimulationManager::print_basis_function_histogram()
{
    std::map<unsigned, unsigned> hist = SimulationManager::histogramize_basis_functions();
    
    printf("Basis function index:\tnumber of representants:\n");
    for(auto p : hist) 
    {
        printf("%10u.        \t%10u \n",p.first,p.second);
    }
}


/* *********************************************************************************************************************************************** *
 * Copies individual with given index and saves in another position.                                                                               *
 * @param unsigned parent_index   - index of individual to copy                                                                                    *
 * @param unsigned child_index    - index where copy will be stored                                                                                *
 *                                                                                                                                                 * 
 * *********************************************************************************************************************************************** */
inline void SimulationManager::reproduce_individual(unsigned parent_index, unsigned child_index)
{
    if ((parent_index > mpopulation)) {  printf("Erorr!!! Too large parent_index: %u \n",parent_index); exit(EXIT_FAILURE);  }
    if ((child_index > mpopulation))  {  printf("Erorr!!! Too large child_index:  %u \n",child_index ); exit(EXIT_FAILURE);  }
    
    // BLAD !!! POD(PIERDALAM SWOJE WLASNE DANE) -PISUJE SWOBIE DANE !!!!
    //this->adaptations[child_index] = adaptations_cpy[parent_index];
    
    parent_index = parent_index*nbasis;
    child_index  = child_index*nbasis;
    
    memcpy( population_mem + child_index,         population_mem_cpy + parent_index,         nbasis         * sizeof(unsigned) );
    memcpy( params_mem     + child_index*nparams, params_mem_cpy     + parent_index*nparams, nbasis*nparams * sizeof(double) );
    memcpy( basis_mem      + child_index,         basis_mem_cpy      + parent_index,         nbasis         * sizeof(ptr_func_t) );
    memcpy( energies_mem   + child_index,         energies_mem_cpy   + parent_index,         nbasis         * sizeof(double) );
    
    indvs_copied = false;
}



/* *********************************************************************************************************************************************** *
 * This function changes basis function (index in basis_assortment, params and pointer to func) of two individuals at given positions.             *
 * The purpose is to use it when crossover.                                                                                                        *
 *                                                                                                                                                 *
 * @param unsigned indv1_idx    - index of 1st individual                                                                                          *
 * @param unsigned indv2_idx    - index of 2nd individual                                                                                          *
 * @param unsigned func1_idx    - index of basis function in 1st idividual to be changed with basis function in 2nd individual                     *
 * @param unsigned func2_idx    - index of basis function in 2nd idividual to be changed with basis function in 1st individual                     *
 *                                                                                                                                                 *
 * *********************************************************************************************************************************************** */
inline void SimulationManager::swap_basis_func(const unsigned indv1_idx, const unsigned indv2_idx, unsigned func1_idx, unsigned func2_idx)
{
    func1_idx = indv1_idx*nbasis + func1_idx;
    func2_idx = indv2_idx*nbasis + func2_idx;
    
    // 2 -> 1
    memcpy( population_mem + func2_idx,         population_mem_cpy + func1_idx,         sizeof(unsigned)         );
    memcpy( params_mem     + func2_idx*nparams, params_mem_cpy     + func1_idx*nparams, sizeof(double) * nparams );
    memcpy( basis_mem      + func2_idx,         basis_mem_cpy      + func1_idx,         sizeof(ptr_func_t)       );
    
    // 1. -> 2.
    //memcpy( population_mem + func2_idx,         population_mem_cpy + func1_idx,         sizeof(unsigned)         );
    //memcpy( params_mem     + func2_idx*nparams, params_mem_cpy     + func1_idx*nparams, sizeof(double) * nparams );
    //memcpy( basis_mem      + func2_idx,         basis_mem_cpy      + func1_idx,         sizeof(ptr_func_t)       );
    
    indvs_copied = false;
}



/* *********************************************************************************************************************************************** *
 * Performs mutation on individula's basis function (changes it with any basis function in basis_assortment).                                      *
 *                                                                                                                                                 *
 * *********************************************************************************************************************************************** */
inline void SimulationManager::mutate_basis_func(const unsigned individual_index, const unsigned basis_index, const unsigned new_wavefunction_index)
{
    basis[individual_index][basis_index]      = basis_assortment[new_wavefunction_index];
    population[individual_index][basis_index] = new_wavefunction_index;
}



inline void SimulationManager::check_basis_func_different()
{
    bool err = false;
    for (unsigned ii=0; ii<mpopulation; ii++)
    {
        for (unsigned jj=0; jj<nbasis; jj++)
        {
            if (population[ii][jj] >= size_assortment) 
            {
                printf("Error! Individual %u.\tfunction %u has too big index %u!\n",ii,jj,population[ii][jj]);
                err = true;
            }
            for (unsigned kk=0; kk<jj; kk++)
                if (population[ii][jj] == population[ii][kk])  { printf("Error! Individual %u.\tfunctions %u and %u the same!\n",ii,jj,kk); err=true; }
        }
    }
    
    if (err) { this->print_functions(); exit(EXIT_FAILURE); }
}


#endif
