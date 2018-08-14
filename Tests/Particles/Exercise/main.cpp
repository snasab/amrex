#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include "DensityParticleContainer.H"

using namespace amrex;

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
 
  ParmParse pp;
  
  TestParams parms;
    
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nppc", parms.nppc);
  pp.get("max_step", parms.max_step);
  pp.get("dt", parms.dt);
  pp.get("max_step", parms.max_step);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");
  
  parms.verbose = false;
  pp.query("verbose", parms.verbose);
  
  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }

  const int lev = 0;

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

// Density multifab does not require a ghost cell
  MultiFab density(ba, dmap, 1,0);
  density.setVal(0.0); 

  
// This multifab is necessary here -- ghost cell // 
  MultiFab partMF(ba, dmap, 1 + BL_SPACEDIM, 1);
  partMF.setVal(0.0);

//  using MyParticleContainer = DensityParticleContainer<1+BL_SPACEDIM> ;
  DensityParticleContainer myPC(geom, dmap, ba);
  myPC.SetVerbose(false);
  
  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
  if (ParallelDescriptor::IOProcessor())
    std::cout << "Total number of particles    : " << num_particles << '\n' << '\n';
    
  bool serialize = true;
  int iseed = 451; 
  Real mass = 10.0;

  myPC.InitParticles(parms);
 // myPC.InitRandom(num_particles, iseed, pdata, serialize);

  // Requires ghost cell for multifab 
  myPC.AssignCellDensitySingleLevel(0, partMF, lev, 4, lev);
  
  //Write particle data -- MAKE THIS AN OPTION (WRITE OR PLT files)  
  for (int i = 0; i < parms.max_step; i++) { 
    myPC.writeParticles(i);
    myPC.moveParticles(parms.dt);
    myPC.Redistribute();
  }
   
 
  amrex::Finalize();
}
