module part_fort_module 
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding, only: c_int

  implicit none 
  private 

  public particle_t, move_particles
  
  type, bind(C) :: particle_t 
    real(amrex_particle_real) :: pos(3) 
    real(amrex_particle_real) :: mass 
    real(amrex_particle_real) :: vel(3)
    integer(c_int)            :: id 
    integer(c_int)            :: cpu
  end type particle_t

contains

   subroutine move_particles(particles, Np, domain_lo, domain_hi, dt) bind(c, name='move_particles') 
     use iso_c_binding, only: c_ptr, c_int, c_f_pointer 
     use amrex_fort_module, only: amrex_real
     use part_fort_module, only: particle_t !Should I declare this? 
  
     implicit none 
   
     type(particle_t), intent(inout), target :: particles(Np) 
     integer(c_int), intent(in) :: Np
     integer(c_int), intent(in) :: domain_lo(3), domain_hi(3) 
     real(amrex_real), intent(in) :: dt 
     type(particle_t), pointer :: p
     
     do i = 1, np
       p => particles(i)

     !Updating positions and velocities 
     !! For now the velocity is constant (but will edit later) 
     !p%vel(1) = ...
     !p%vel(2) = ...
     !p%vel(3) = ...

     p%pos(1) = p%pos(1) + p%vel(1) * dt 
     p%pos(2) = p%pos(2) + p%vel(2) * dt
     p%pos(3) = p%pos(3) + p%vel(3) * dt

     !Cooler things later ...
     enddo 
  end subroutine move_particles

end module part_fort_module
