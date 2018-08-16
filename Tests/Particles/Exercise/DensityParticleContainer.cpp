//add some header files here
#include <iostream> 
#include "DensityParticleContainer.H"
#include "density_F.H"
//No Params header here.. because other header uses it? 

using namespace amrex; 

DensityParticleContainer:: 
DensityParticleContainer(const Geometry             &geom,
                         const DistributionMapping  &dmap,
                         const BoxArray             &ba)
                        : ParticleContainer<1+BL_SPACEDIM>(geom, dmap, ba)
{}

void DensityParticleContainer::InitParticles(TestParams& parms){
  BL_PROFILE("DensityParticleContainer::InitParticles");
  
  //TestParams parms; 
  int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;

  const int lev = 0;
  const Geometry& geom = Geom(lev);
  const Real* dx = geom.CellSize();
  
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi){
//    //Declare no particles first
    auto& particles = GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];  
    for (int i = 0; i < num_particles; ++i){ 
           // std::cout << "ME: " << i << std::endl;
            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            p.pos(0) = float(i)/float(num_particles);
            p.pos(1) = float(i)/float(num_particles);
#if (BL_SPACEDIM == 3)
            p.pos(2) = float(i)/float(num_particles);
            p.rdata(3) = 3.0; //uz
#endif
            p.rdata(1) = 1.0;//ux //XX Maybe parm this? 
            p.rdata(2) = 2.0; //uy
            p.rdata(0) = 10.0;//mass
            
            particles.push_back(p);
     }


}}

void DensityParticleContainer::moveParticles(const Real dt){
    BL_PROFILE("DensityParticleContainer::moveParticles");
    
    TestParams parms;
    
    const int lev = 0; 
    const RealBox& prob_domain = Geom(lev).ProbDomain();  
 
    //Read in parms
    IntVect domain_lo(0 , 0, 0);
    IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1);

    for (DParIter pti(*this, lev); pti.isValid(); ++pti){
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();
        amrex_move_particles(particles.data(), &Np, &dt, prob_domain.lo(), prob_domain.hi());
    }

}

/*void DensityParticleContainer::push_momentum(const amrex::MultiFab&U,
                                             const amrex::Real    dt)
{
    BL_PROFILE("DensityParticleContainer::push_momentum");

    

void DensityParticleContainer::push_momentum(const amrex::MultiFab& ux,
                                             const amrex::MultiFab& uy,
                                             const amrex::MultiFab& uz,
                                                   amrex::Real      dt)
{
    BL_PROFILE("DensityParticleContainer::push_momentum");
    
    const int lev = 0;
    for (DParIter pti(*this, lev); pti.isValid(); ++pti){
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();
        amrex_push_momentum(particles.data(), &Np, ux.data(), uy.data(), uz.data(), &dt);
     integer(c_int), intent(in) :: Np
    }
} */



 

void DensityParticleContainer::writeParticles(int n){
  
    BL_PROFILE("DensityParticleContainer::writeParticles");

    const std::string& pltfile = amrex::Concatenate("particles",n,5);
    //Concatenate(plot_file_root,level_steps[0],file_name_digits);
    WriteAsciiFile(pltfile);
}




























