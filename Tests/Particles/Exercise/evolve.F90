module evolve_module 
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding, only: c_int

  implicit none 
  private 

  public particle_t
  
  type, bind(C) :: particle_t 
    real(amrex_particle_real) :: pos(3) 
    real(amrex_particle_real) :: mass 
    real(amrex_particle_real) :: vel(3)
    integer(c_int)            :: id 
    integer(c_int)            :: cpu
  end type particle_t
end module evolve_module

!contains

   subroutine amrex_push_momentum(particles, Np, ux, uy, uz, dt) bind(c, name='amrex_push_momentum') 
     use iso_c_binding, only: c_ptr, c_int, c_f_pointer 
     use amrex_fort_module, only: amrex_real
     use part_fort_module, only: particle_t 
  
     implicit none 
      
     integer(c_int), intent(in) :: Np
     type(particle_t), intent(inout), target :: particles(Np) 
     real(amrex_real), intent(inout) :: ux(Np), uy(Np), uz(Np)
     type(particle_t), pointer :: p 
     real(amrex_real), intent(in) :: dt
     real, parameter :: gravity=9.8   
     integer :: i
     
     do i = 1, Np
       p => particles(i) 
       !p%pos(1) = p%pos(1) +  ux(ip)*dt
       !p%pos(2) = p%pos(2) +  uy(ip)*dt
       !p%pos(3) = p%pos(3) +  uz(ip)*dt
     !velocity is up = u - Ws, where Ws = "mass" * gravity
       p%vel(1) = ux(i) 
       p%vel(2) = uy(i)
       p%vel(3) = uz(i) - p%mass*gravity
   
       p%pos(1) = p%vel(1)*dt
       p%pos(2) = p%vel(2)*dt
       p%pos(3) = p%vel(3)*dt
     enddo 
  end subroutine amrex_push_momentum
