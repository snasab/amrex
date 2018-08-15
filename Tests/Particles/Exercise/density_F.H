#ifndef DENSITY_F_H_
#define DENSITY_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void init_fluid(const int* domain_lo, const int* domain_hi, amrex_real* fluidMF,
                    const int* dlo, const int* dhi); 
    void amrex_move_particles(void* particles, const int* Np,
		              const amrex::Real* dt,
		              const amrex::Real* domain_lo,
		              const amrex::Real* domain_hi);


#ifdef __cplusplus
}
#endif

#endif 
