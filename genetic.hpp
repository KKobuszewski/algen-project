#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include <random>
#include <algorithm>
#include <unordered_set>


#include "SimulationManager.hpp"


#define GENETIC_SUCCESS ((genetic_result_t) 0)

typedef int genetic_result_t;

namespace Genetic
{

double tot_adapation;
double beta = 10.;
double impact_factor = 1.;
unsigned n_energies = 1;

Random genetic_rng;


inline void set_selection_pressure(const double _beta)
{
    beta = _beta;
}

inline void set_impact_factor(const double _alpha)
{
    impact_factor = _alpha;
}

inline void set_levels_for_adaptation(const unsigned _n)
{
    n_energies = _n;
}


// ===========================================  SOME FANCY UTILITIES  =====================================================================



template <typename T>
inline void unique_elements_indices(T* arr, T* values, unsigned *indices, const unsigned n_arr, const unsigned n_val )
{
    for (unsigned ii=0; ii< n_val; ii++)
        indices[ii] = std::distance(arr, std::find(arr, arr + n_arr, values[ii]));
}


/*
 * This function gives all distinct elements from two given arrays.
 * 
 */
template <typename T>
inline void unique_merge(T* a1, T* a2, T* result, unsigned *n1, unsigned *n2, const unsigned size1, unsigned size2 = 0)
{
    if (size2 ==0) size2 = size1;
    T *p = result;
    p = std::copy_if( a1, a1+size1, p, [&]( T x ) { return std::find( a2, a2+size2, x ) == a2+size2; } );
    
    *n1 = p - result;
    
    p = std::copy_if( a2, a2+size2, p, [&]( T x ) { return std::find( a1, a1+size1, x ) == a1+size1; } );
    
    *n2 = p - result - *n1;
}





// ==========================================  ADAPTATION DETERMINATION  ==================================================================


/*
 * E -     energies
 * beta -  selection preassure/temperature
 * n -     number of energy levels taken into accont
 */
inline double adaptation_function(double* E, const unsigned n, const double beta, const double impact_factor=1.)
{
    double adaptation = 0.;
    #pragma omp for simd
    for (unsigned ii=0; ii<n; ii++)
        adaptation += E[ii]/(impact_factor*ii+1);
    if ( (exp( -1.*beta*adaptation ) < 1e-5) || (exp( -1.*beta*adaptation ) > 1e+10) )
        printf("f: %.3lf, sum:%.3lf, E0: %.3lf, beta:%.3lf\n",exp( -1.*beta*adaptation ),adaptation,E[0],beta);
    return exp( -1.*beta*adaptation );
}

inline genetic_result_t get_probablities(SimulationManager *sim)
{
    double *ad = sim->adaptations;
    const unsigned m = sim->mpopulation;
    
    #pragma omp for simd
    for(unsigned ii=0; ii<m; ii++)
    {
        ad[ii] = adaptation_function(sim->energies[ii],n_energies,beta,impact_factor);
    }
    
    #pragma omp for simd reduction(+:tot_adapation)
    for(unsigned ii=0; ii<m; ii++)
    {
        tot_adapation += ad[ii];
    }
    
    #pragma omp for simd
    for(unsigned ii=0; ii<m; ii++)
    {
        ad[ii] /= tot_adapation;
    }
    
    return GENETIC_SUCCESS;
}

inline genetic_result_t get_probablities_cumulative(SimulationManager *sim)
{
    double *ad = sim->adaptations;
    const unsigned m = sim->mpopulation;
    
    
    ad[0] = adaptation_function(sim->energies[0],n_energies,beta,impact_factor);
    #pragma omp for simd
    for(unsigned ii=1; ii<m; ii++)
    {
        ad[ii] = ad[ii-1] + adaptation_function(sim->energies[ii],n_energies,beta,impact_factor);
    }
    
    tot_adapation = ad[m-1];
    
    #pragma omp for simd
    for(unsigned ii=0; ii<m; ii++)
    {
        sim->adaptations[ii] /= tot_adapation;
    }
    
    return GENETIC_SUCCESS;
}



// ==========================================  REPRODUCTION  ==================================================================


inline genetic_result_t random_selection(SimulationManager *sim)
{
    const unsigned m = sim->mpopulation;
    // get probabilities of 
    // get_probablities_cumulative(sim);
    // double p[SimulationLimits::max_population];
    //size_t p_int[SimulationLimits::max_population];
    
    
    
    
    
    
    /*for (unsigned ii=0; ii<m; ii++)
    {
        p_int[ii] = (size_t) 
    }
    
    std::discrete_distribution<> dist(sim->adaptations);
    
    printf("indices: \n");
    //#pragma omp parallel for
    for (unsigned ii=0; ii<m; ii++)
    {
        const unsigned index=dist(genetic_rng.mersenne);
        
        printf("index: %u  ",index);
        sim->reproduce_individual(index,ii);
        printf("\n");
    }*/
    printf("\n");
}

inline genetic_result_t random_selection_patryk(SimulationManager *sim)
{
    double *ad = sim->adaptations_cpy;
    const unsigned m = sim->mpopulation;
    double sum = 0.0f;
    
    // Obliczam funkcje celu:
    #pragma omp for simd
    for(unsigned int ii=0; ii<m; ii++)
    {
        ad[ii] = adaptation_function(sim->energies[ii],n_energies,beta,impact_factor);
        sum += ad[ii];
    }
    printf("min PDF: %.3lf\n",*std::min_element(ad,ad+m));
    printf("max PDF: %.3lf\n",*std::max_element(ad,ad+m));
    printf("sum PDF: %.3lf\n",sum);
    
    // Obliczam ulamki
    #pragma omp for simd
    for(unsigned int ii=0; ii<m; ii++)
    {
        ad[ii] /= sum;
    }

    // Obliczam dystrybuante;
    #pragma omp for simd
    for(unsigned int ii=1; ii<m; ii++)
    {
        ad[ii] += ad[ii - 1];
    }
    
    printf("selection/CDF: [");
    for (unsigned int ii=0; ii<m; ii++)  printf("%.3lf ",ad[ii]);
    printf("]\n");
    
    // 
    double rand = 0.0;
    #pragma omp for simd
    for(unsigned int ii=1; ii<m; ii++)
    {
        rand = genetic_rng.get_random_uniform();        
        //ad[ii] += ad[ii - 1];
        
        int jj;
        for ( jj = 0; jj < m; jj++) {
            if (rand < ad[jj]) break;
        }
        printf("reproducing individual %u. (%.3lf/%.3lf)\n",jj,rand,ad[jj]);
        sim->reproduce_individual(jj,ii);
    }
}


template <const unsigned size_tournament>
inline genetic_result_t tournament_selection(SimulationManager *sim)
{
    get_probablities(sim);
    unsigned indices[size_tournament];
    const unsigned m = sim->mpopulation;
    printf("tournament_selection for %u individuals\n",size_tournament);
    
    //#pragma omp parallel for
    for (unsigned ii=0; ii < m; ii++)
    {
        // draw size_tournament indices of individuals
        indices[0] = rand() % m;
        for (unsigned jj=1; jj < size_tournament; jj++)
        {
            indices[jj] = (rand() % m);
            for (unsigned kk=0; kk < jj; kk++)
                while (indices[jj] == indices[kk]) {  indices[jj] = rand() % m;  }
        }
        
        
        // check which drawn specimen is best
        unsigned index_chosen = 0;
        #pragma omp for simd
        for (unsigned jj = 1; jj < size_tournament; jj++)
        {
            if ( sim->adaptations[indices[jj]] > sim->adaptations[index_chosen] )
                index_chosen = indices[jj];
        }
        //printf("index: %u\n",index_chosen);
        
        // copy
        sim->reproduce_individual(index_chosen,ii);
        
        //if (ii==0)
        //{
        //    sim->print_matrices(index_chosen);
        //    sim->print_matrices(ii);
        //}
    }
    return GENETIC_SUCCESS;
}


// ==========================================  RECOMBINATION  ==================================================================


#define MAX_ANCESTORS 10

typedef unsigned crossover_type_t;

#define CROSSOVER_MODE_1 ((Genetic::crossover_type_t) 1)
#define CROSSOVER_MODE_2 ((Genetic::crossover_type_t) 2)


unsigned indices_basis_functions[SimulationLimits::max_nbasis * MAX_ANCESTORS];
unsigned indices_in_populations [SimulationLimits::max_nbasis * MAX_ANCESTORS];
unsigned population_perm[SimulationLimits::max_nbasis * MAX_ANCESTORS];
unsigned random_perm[SimulationLimits::max_nbasis];

/*
ancestor_func_indices[2 * MAX_ANCESTORS * SimulationLimits::max_nbasis * SimulationLimits::max_threads];
unsigned ancestor_indices[MAX_ANCESTORS];

inline void crossover_mode1(SimulationManager *sim, const unsigned ancestor_num, const double pc)
{
    const unsigned nbasis = sim->nbasis;
    const unsigned m      = sim->mpopulation;
    
    //#pragma omp parallel for
    for (unsigned ii=0; ii < m; ii++)
    {
        if ( pc < genetic_rng.get_random_uniform() )
        {
            // choose random indices of ancestors
            for (unsigned jj=0; jj<ancestor_num; jj++)
            {
                ancestor_indices[jj] = genetic_rng.random_integer(0,sim->size_assortment-1);
            }
            
            // rewrite indices of basis functions 
            for (unsigned jj=0; jj<ancestor_num; jj++)
            for (unsigned kk=0; kk<nbasis      ; kk++)
            {
                const unsigned index = ancestor_indices[jj];
                ancestor_func_indices[jj*nbasis+kk] = sim->population[index][kk];
            }
            
            // find unique indices in ancestors
            std::unordered_set<unsigned> set(ancestor_func_indices, ancestor_func_indices + nbasis*ancestor_num);
            const unsigned unique_indices_num = set.size();
            std::copy(set.begin(), set.end(), ancestor_func_indices);
            
            // shuffle indices
            genetic_rng.shuffle_array<unsigned>(ancestor_func_indices,unique_indices_num);
            std::sort(ancestor_func_indices,ancestor_func_indices+nbasis);
            
            // create individual from indices
            
        }
    }
    
    
}

inline void crossover_mode2(SimulationManager *sim, const unsigned ancestor_num, const double pc)
{
    const unsigned nbasis = sim->nbasis;
    const unsigned m      = sim->mpopulation;
    
    //#pragma omp parallel for
    for (unsigned ii=0; ii < m; ii++)
    {
        if ( pc < genetic_rng.get_random_uniform() )
        {
            // choose random indices of ancestors
            for (unsigned jj=0; jj<ancestor_num; jj++)
            {
                ancestor_indices[jj] = genetic_rng.random_integer(0,sim->size_assortment-1);
            }
            
            // rewrite indices of basis functions 
            for (unsigned jj=0; jj<ancestor_num; jj++)
            for (unsigned kk=0; kk<nbasis      ; kk++)
            {
                const unsigned index = ancestor_indices[jj];
                ancestor_func_indices[jj*nbasis+kk] = sim->population[index][kk];
            }
            
            
        }
    }
}
*/


/*
 * @param double pc  - probability of crossover for an individual => p_c = <C>/M, and we need <C>/(M(M-1)/2) = 2*p_c/(M-1)
 * 
 */
inline genetic_result_t crossover_mode1(SimulationManager *sim, const double pc)
{
    const unsigned nbasis = sim->nbasis;
    const unsigned m      = sim->mpopulation;
    
    for (unsigned ii=0; ii < nbasis; ii++)    random_perm[ii] = ii;
    //printf("pc = %.5lf\n",pc*2./(m-1.));
    
    
    // iterate over each pair of individuals
    for (unsigned ii=0; ii < m; ii++)
    for (unsigned jj=0; jj < ii; jj++)
    {
        
        unsigned* population1 = sim->population[ii];
        unsigned* population2 = sim->population[jj];
        
        if ( pc*2./(m-1.) > genetic_rng.get_random_uniform() )
        {
            printf("Crossver individuals %u & %u\t",ii,jj);
            
            // find distinct indices of basis functions in the individuals
            unsigned n1, n2;
            unique_merge<unsigned>( population1, population2, indices_basis_functions, &n1, &n2, nbasis );
            unique_elements_indices<unsigned>( population1, indices_basis_functions, indices_in_populations, nbasis, n1 );
            unique_elements_indices<unsigned>( population2, indices_basis_functions+n1, indices_in_populations+n1, nbasis, n2 );
            
            printf("n1:%u, n2:%u\t",n1,n2);
            
            // TODO: MAKE SURE THAT ALL FUNCTIONS INDICES IN INDIVIDUAL ARE UNIQUE
            for (unsigned kk=0; kk<n1+n2; kk++)  printf("%u ",indices_basis_functions[kk]);
            printf("\n");
            
            genetic_rng.shuffle_array<unsigned>(random_perm,nbasis);
            
            // swap indices from 1. individual to 2. individual
            for (unsigned kk=0; kk<n1; kk++)
            {
                //unsigned idx_to_swap = genetic_rng.random_integer(0,nbasis-1);
                unsigned idx_to_swap = random_perm[kk];
                sim->swap_basis_func(ii,jj,indices_in_populations[kk],idx_to_swap);
            }
            
            genetic_rng.shuffle_array<unsigned>(random_perm,nbasis);
            // swap indices from 2. individual to 1. individual
            for (unsigned kk=0; kk<n2; kk++)
            {
                //unsigned idx_to_swap = genetic_rng.random_integer(0,nbasis-1);
                unsigned idx_to_swap = random_perm[kk];
                sim->swap_basis_func(jj,ii,idx_to_swap,(indices_in_populations+n1)[kk]);
            }
            
            sim->print_matrices(ii);
            sim->print_matrices(jj);
        }
    }
    
    
    return GENETIC_SUCCESS;
}

/*
 * @param double pc  - probability of crossover for an individual => p_c = <C>/M, and we need <C>/(M(M-1)/2) = 2*p_c/(M-1)
 * 
 */
inline genetic_result_t crossover_mode2(SimulationManager *sim, const double pc)
{
    const unsigned nbasis = sim->nbasis;
    const unsigned m      = sim->mpopulation;
    
    for (unsigned ii=0; ii < nbasis; ii++)    random_perm[ii] = ii;
    
    for (unsigned ii=0; ii < m; ii++)         population_perm[ii] = ii;
    genetic_rng.shuffle_array<unsigned>(population_perm,m);
    
    // iterate over each pair of individuals
    for (unsigned ii=0; ii < m/2; ii++)
    {
        unsigned jj = m - ii - 1;
        unsigned* population1 = sim->population[ population_perm[ii] ];
        unsigned* population2 = sim->population[ population_perm[jj] ];
        
        if ( pc > genetic_rng.get_random_uniform() )                // shouldn't be pc*2.? -> rather not, because it could be grater than 1.
        {
            printf("Crossver individuals %u & %u\t",ii,jj);
            
            // find distinct indices of basis functions in the individuals
            unsigned n1 = nbasis + 1;
            unsigned n2 = nbasis + 1;
            unique_merge<unsigned>( population1, population2, indices_basis_functions, &n1, &n2, nbasis );
            unique_elements_indices<unsigned>( population1, indices_basis_functions, indices_in_populations, nbasis, n1 );
            unique_elements_indices<unsigned>( population2, indices_basis_functions+n1, indices_in_populations+n1, nbasis, n2 );
            
            //printf("n1:%u, n2:%u\t",n1,n2);
            
            // TODO: MAKE SURE THAT ALL FUNCTIONS INDICES IN INDIVIDUAL ARE UNIQUE
            for (unsigned kk=0; kk<n1+n2; kk++)  printf("%u ",indices_basis_functions[kk]);
            printf("\n");
            
            genetic_rng.shuffle_array<unsigned>(random_perm,nbasis);
            
            // swap indices from 1. individual to 2. individual
            for (unsigned kk=0; kk<n1; kk++)
            {
                //unsigned idx_to_swap = genetic_rng.random_integer(0,nbasis-1);
                unsigned idx_to_swap = random_perm[kk];
                sim->swap_basis_func(population_perm[ii],population_perm[jj],indices_in_populations[kk],idx_to_swap);
                //printf("crossover (%u,%u): %u -> %u\n",ii,jj,indices_in_populations[kk],idx_to_swap);
            }
            
            genetic_rng.shuffle_array<unsigned>(random_perm,nbasis);
            // swap indices from 2. individual to 1. individual
            for (unsigned kk=0; kk<n2; kk++)
            {
                //unsigned idx_to_swap = genetic_rng.random_integer(0,nbasis-1);
                unsigned idx_to_swap = random_perm[kk];
                sim->swap_basis_func(population_perm[jj],population_perm[ii],(indices_in_populations+n1)[kk],idx_to_swap);
                //printf("crossover (%u,%u): %u -> %u\n",jj,ii,idx_to_swap,(indices_in_populations+n1)[kk]);
            }
            
            //sim->print_matrices(ii);
            //sim->print_matrices(jj);
        }
    }
    
    return GENETIC_SUCCESS;
}



// , unsigned ancestor_num
template <crossover_type_t crossover_mode>
inline genetic_result_t recombine(SimulationManager *sim, const double pc)
{
    //if (ancestor_num > MAX_ANCESTORS) {  printf("To much ancestors given! (%u/%u)\n",ancestor_num,MAX_ANCESTORS); exit(EXIT_FAILURE);  }
    
    switch(crossover_mode)
    {
        case CROSSOVER_MODE_1:
        {
            crossover_mode1(sim, pc);
        }
        case CROSSOVER_MODE_2:
        {
            crossover_mode2(sim, pc);
        }
    }
    
    return GENETIC_SUCCESS;
}



// ==========================================  MUTATIONS  ==================================================================


inline genetic_result_t mutate(SimulationManager *sim, const double pm)
{
    double   rand;
    unsigned index, new_wavefunction_index;
    const unsigned m = sim->mpopulation;
    const unsigned n = sim->nbasis;
    const unsigned L = sim->size_assortment;
    
    for (unsigned ii=0; ii<m; ii++)
    {
        rand = genetic_rng.get_random_uniform();
        
        // change random basis function if drawn number is smaller than probability of mutation
        if ( rand < pm )
        {
            bool mutation_forbidden;
            
            do
            {
                // random index in individual basis
                new_wavefunction_index = genetic_rng.random_integer(0,L-1);
                
                // check if mutation is valid <- do not want to have the same wavefunctions
                mutation_forbidden = false;
                for (unsigned jj=0; jj<n; jj++)
                    if ( sim->population[ii][jj] == new_wavefunction_index ) {  mutation_forbidden = true; break;  }
            } while(mutation_forbidden);  // means mutation will change
            
            index = genetic_rng.random_integer(0,n-1);
            printf("(%u. individual)\tMutation: %u.(%u) - > %u \tgenotype before: [",ii,index,sim->population[ii][index],new_wavefunction_index);
            for (unsigned jj=0; jj<n; jj++) printf("%u,",sim->population[ii][jj]);
            
            printf("]\t");
            // made change in basis of individual's wavefunction
            sim->mutate_basis_func(ii,index,new_wavefunction_index);
            printf("genotype after: [");
            for (unsigned jj=0; jj<n; jj++) printf("%u,",sim->population[ii][jj]);
            printf("]\n");
        }
    }
    
    return GENETIC_SUCCESS;
}



}