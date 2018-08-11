#include <iostream>
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include "myfunc.H"
#include "AMReX_TracerParticles.H"

using namespace amrex;

void test_assign_density(TestParams& parms)
{
  const int lev = 0;
  Real dt;
  int max_step; 
//INITIALIZE
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
  //END INITIALIZE

  MultiFab::Copy(density, partMF, 0, 0, 1, 0);

  myPC.WriteAsciiFile("particles");

  myPC.MoveRandom(lev);
  myPC.WriteAsciiFile("particles1");

  myPC.MoveRandom(lev);
  myPC.WriteAsciiFile("particles2");
}
