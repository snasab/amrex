#ifndef _DENSITY_PARTICLE_CONTAINER_H_
#define _DENSITY_PARTICLE_CONTAINER_H_

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Params.H>

namespace amrex{

class DensityParticleContainer
    : public ParticleContainer<1+BL_SPACEDIM> //For mass and velocities (2 or 3) 
    {	
public:
    using DParIter	= ParIter<1+BL_SPACEDIM>;

    DensityParticleContainer(const Geometry 		&geom,
                             const DistributionMapping	&dmap,
		             const BoxArray		&ba);
    
    void InitParticles(TestParams& parms);
    

    void moveParticles(const Real dt);    
    
    void writeParticles(int n);
    };

}



#endif // _DENSITY_PARTICLE_CONTAINER_H_
