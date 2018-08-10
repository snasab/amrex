//add some header files here 

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
   	
}
*/
