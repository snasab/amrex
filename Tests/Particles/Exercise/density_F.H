#ifndef DENSITY_F_H_
#define DENSITY_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void ux_init(const int* domain_lo, const int* domain_hi, amrex_real* ux_MultiFab,
                    const int* dlo, const int* dhi);

    void uy_init(const int* domain_lo, const int* domain_hi, amrex_real* uy_MultiFab,
                    const int* dlo, const int* dhi);

    void uz_init(const int* domain_lo, const int* domain_hi, amrex_real* uz_MultiFab,
                    const int* dlo, const int* dhi); 

    void amrex_init_fluid(const int* domain_lo, const int* domain_hi, amrex_real* U_MultiFab,
                    const int* dlo, const int* dhi,  int ncU);
    
    void amrex_push_momentum(void* particles, const int* Np, amrex::Real* ux, 
         amrex::Real* uy, amrex::Real* uz, const amrex::Real* dt);

    void amrex_move_particles(void* particles, const int* Np,
		              const amrex::Real* dt,
		              const amrex::Real* domain_lo,
		              const amrex::Real* domain_hi);


#ifdef __cplusplus
}
#endif

#endif 
