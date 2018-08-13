#ifndef DENSITY_F_H_
#define DENSITY_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_move_particles(void* particles, const int* Np,
		              const amrex::Real* dt,
		              const amrex::Real* domain_lo,
		              const amrex::Real* domain_hi);


#ifdef __cplusplus
}
#endif

#endif 
