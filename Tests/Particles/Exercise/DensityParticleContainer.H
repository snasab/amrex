#ifndef _DENSITY_PARTICLE_CONTAINER_H_
#define _DENSITY_PARTICLE_CONTAINER_H_

//#include <AMReX_Particles.H>

namespace amrex{

class DensityParticleContainer
    : public ParticleContainer<1+BL_SPACEDIM> //For mass and velocities (2 or 3) 
    {	
public:
//dont need this?    using DParIter	= ParIter<1+BL_SPACEDIM>;

    DensityParticleContainer(const Geometry 		&geom,
                             const DistributionMapping	&dmap,
		             const BoxArray		&ba);
    
    //void InitParticles( ...what do I put here? 
    

    void moveParticles(const Real dt);
    

    };

}

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  bool verbose;
};


#endif // _DENSITY_PARTICLE_CONTAINER_H_
