module part_fort_module 
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
end module part_fort_module

!contains

   subroutine amrex_move_particles(particles, Np, dt, domain_lo, domain_hi) bind(c, name='amrex_move_particles') 
     use iso_c_binding, only: c_ptr, c_int, c_f_pointer 
     use amrex_fort_module, only: amrex_real
     use part_fort_module, only: particle_t 
  
     implicit none 
   
     type(particle_t), intent(inout), target :: particles(Np) 
     integer(c_int), intent(in) :: Np
     integer(c_int), intent(in) :: domain_lo(3), domain_hi(3) 
     real(amrex_real), intent(in) :: dt  !intent(in), value:: 
     type(particle_t), pointer :: p
     integer :: i
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
     enddo 
  end subroutine amrex_move_particles

!end module part_fort_module
