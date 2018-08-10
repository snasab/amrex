#ifndef _DENSITY_PARTICLE_CONTAINER_H_
#define _DENSITY_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

namespace amrex{

class DensityParticleContainer
    : public ParticleContainer<1+BL_SPACEDIM> //For mass and velocities (2 or 3) 
    {	
public:
    using DParIter	= ParIter<1+BL_SPACEDIM>;

    DensityParticleContainer(const Geometry 		&geom,
                             const DistributionMapping	&dmap,
		             const BoxArray		&ba);
    
    //void InitParticles( ...what do I put here? 
    //moveParticles(...
    };

}

#endif // _DENSITY_PARTICLE_CONTAINER_H_
