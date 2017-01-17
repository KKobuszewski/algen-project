#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include <stdlib.h>
#include <math.h>

#include <complex>


typedef std::complex<double> (*ptr_integrand_t)(ptr_func_t,ptr_func_t,double*,double*, unsigned, unsigned, unsigned);

class Integrator
{
private:
    std::complex<double>* cplx_func_on_lattice;
    double*               dbl_func_on_lattice;
    unsigned old_samples;
public:
    // constructors
    Integrator();
    
    // destructors
    ~Integrator();
    
    // methods
    std::complex<double> integrate(ptr_integrand_t func, ptr_func_t psi_i,ptr_func_t psi_j, 
                                   double* params, unsigned ndims, unsigned nparams, unsigned npart,
                                   unsigned samples, double* xl, double* xu);
};

Integrator::Integrator()
{
    cplx_func_on_lattice = NULL;
    dbl_func_on_lattice = NULL;
}

Integrator::~Integrator()
{
    if (cplx_func_on_lattice) free(cplx_func_on_lattice);
    if (dbl_func_on_lattice)  free(dbl_func_on_lattice);
}



/*
 * 
 * @params
 * func          - functional of psi_i & psi_j to be integrated
 * psi_i, psi_j  - two basis functions, could be different
 * 
 * 
 * 
 */
inline std::complex<double> Integrator::integrate(
    ptr_integrand_t func, ptr_func_t psi_i,ptr_func_t psi_j, 
    double* params, unsigned ndims, unsigned nparams, unsigned npart,
    unsigned samples, double* xl, double* xu)
{
    const unsigned tot_dims = ndims*npart;
//     printf("samples:   %u\n",samples);
//     printf("ndims:     %u\n",ndims);
//     printf("npart:     %u\n",npart);
//     printf("nparams:   %u\n",nparams);
    
    if (tot_dims <= 3)  // if small number of phase-space dimension use typical integration
    {
        // check if needed memory is allocable
        if (pow(samples,tot_dims)*sizeof(std::complex<double>) > 1073741824LU)
        {
            printf("ERROR (Integrator::integrate):trying to allocate more that 1073741824 bytes! (%lf)\n",
                   pow(samples,tot_dims)*sizeof(std::complex<double>));
            exit(EXIT_FAILURE);
        }
        
        // allocate/reallocate memory if needed    -> could be helpful only if using FFT for derivatives!
        if (!cplx_func_on_lattice) 
        {
            cplx_func_on_lattice = (std::complex<double>*) malloc( pow(samples,tot_dims)*sizeof(std::complex<double>) );
            old_samples = samples;
        }
        else if (samples > old_samples)
        {
            cplx_func_on_lattice = (std::complex<double>*) realloc( cplx_func_on_lattice, pow(samples,tot_dims)*sizeof(std::complex<double>) );
            old_samples = samples;
        }
        else if (samples < old_samples)
        {
            cplx_func_on_lattice = (std::complex<double>*) realloc( cplx_func_on_lattice, pow(samples,tot_dims)*sizeof(std::complex<double>) );
        }
        
        //printf("tot_dims:  %u\n",tot_dims);
        
        switch (tot_dims)
        {
        case 1:
            {
                const double L = *xu - *xl;
                const double dx = L/((double) samples);
                
                for (unsigned ix = 0; ix < samples; ix++)
                {            
                    double x[1] = {ix*dx - L/2.};
                    cplx_func_on_lattice[ix] = func(psi_i, psi_j, x, params, ndims, nparams, npart);
                }
                
                std::complex<double> integ_result(0.,0.);
                for (unsigned ix = 0; ix < samples; ix++) integ_result += cplx_func_on_lattice[ix];
                return integ_result;
            }
        case 2:
            {
                const double Lx = xu[0] - xl[0];
                const double Ly = xu[1] - xl[1];
                const double dx = Lx/((double) samples);
                const double dy = Ly/((double) samples);
                
    //             printf("\n!!!!!!!!!!!!!!!!!!!!!!!  2D not implemented yet  !!!!!!!!!!!!!!!!!!!!!!!\n");
    //             exit(EXIT_FAILURE);
                for (unsigned ixy = 0; ixy < samples*samples; ixy++)
                {
                    unsigned ix = 0;
                    unsigned iy = 0;
                    double x[2] = {ix*dx - Lx/2.,iy*dy - Ly/2.};
                    cplx_func_on_lattice[ixy] = func(psi_i, psi_j, x, params, ndims, nparams, npart);
                }
                
                std::complex<double> integ_result(0.,0.);
                for (unsigned ixy = 0; ixy < samples*samples; ixy++) integ_result += cplx_func_on_lattice[ixy];
                return integ_result;
            }
        case 3:
            {
                const double Lx = xu[0] - xl[0];
                const double Ly = xu[1] - xl[1];
                const double Lz = xu[2] - xl[2];
                const double dx = Lx/((double) samples);
                const double dy = Ly/((double) samples);
                const double dz = Lz/((double) samples);
                
    //             printf("\n!!!!!!!!!!!!!!!!!!!!!!!  3D not implemented yet  !!!!!!!!!!!!!!!!!!!!!!!\n");
    //             exit(EXIT_FAILURE);
                for (unsigned ixyz = 0; ixyz < samples*samples*samples; ixyz++)
                {
                    unsigned ix = 0;
                    unsigned iy = 0;
                    unsigned iz = 0;
                    double x[3] = {ix*dx - Lx/2.,iy*dy - Ly/2., iz*dz - Lz/2.};
                    cplx_func_on_lattice[ixyz] = func(psi_i, psi_j, x, params, ndims, nparams, npart);
                }
                
                std::complex<double> integ_result(0.,0.);
                for (unsigned ixyz = 0; ixyz < samples*samples*samples; ixyz++) integ_result += cplx_func_on_lattice[ixyz];
                return integ_result;
            }
        }
    }
    else                         // use Monte-Carlo integration
    {
        printf("\n!!!!!!!!!!!!!!!!!!!!!!!  MC not implemented yet  !!!!!!!!!!!!!!!!!!!!!!!\n");
        exit(EXIT_FAILURE);
    }
    
//     printf("\n");
}




#endif