//add some header files here 
#include <iostream>
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include "DensityParticleContainer.H"

using namespace amrex; 

DensityParticleContainer:: 
DensityParticleContainer(const Geometry             &geom,
                         const DistributionMapping  &dmap,
                         const BoxArray             &ba)
                        : ParticleContainer<1+BL_SPACEDIM>(geom, dmap, ba)
{}

/*void DensityParticleContainer::InitParticles(...same as header inputs ..){
  BL_PROFILE("DensityParticleContainer::InitParticles");
  
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  IntVect domain_lo(0 , 0, 0);
  IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1);
  const Box domain(domain_lo, domain_hi);

  // This says we are using Cartesian coordinates
  int coord = 0;

  // This sets the boundary conditions to be doubly or triply periodic
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++)
    is_per[i] = 1;
  Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

  BoxArray ba(domain);
  ba.maxSize(parms.max_grid_size);
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "Number of boxes              : " << ba[0].size() << '\n' << '\n';
  }

  DistributionMapping dmap(ba);
  
  MultiFab density(ba, dmap, 1, 0);
  density.setVal(0.0);

  MultiFab partMF(ba, dmap, 1 + BL_SPACEDIM, 1);
  partMF.setVal(0.0);

  using MyParticleContainer = ParticleContainer<1+BL_SPACEDIM> ;
  MyParticleContainer myPC(geom, dmap, ba);
  myPC.SetVerbose(false);

  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';

  bool serialize = true;
  int iseed = 451;
  Real mass = 10.0;

  MyParticleContainer::ParticleInitData pdata = {mass, 1.0, 2.0, 3.0};
  myPC.InitRandom(num_particles, iseed, pdata, serialize);
  myPC.AssignCellDensitySingleLevel(0, partMF, 0, 4, 0);
  	
}*/

void DensityParticleContainer::moveParticles(const Real dt){
    BL_PROFILE("DensityParticleContainer::moveParticles");
    
    TestParams parms;
    
    const int lev = 0; 
    
    //Read in parms
    IntVect domain_lo(0 , 0, 0);
    IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1);

    for (ParIter pti(*this, lev); pti.isValid(); ++pti){
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();
        move_particles(particles.data(), &Np, &dt, domain_lo, domain_hi);
    }
}




























