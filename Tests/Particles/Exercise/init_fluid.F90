subroutine init_fluid(domain_lo, domain_hi) bind(C, name="init_fluid")

  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_int

  implicit none

  integer(c_int)  , intent(in)    :: domain_lo(3), domain_hi(3) 
  real(amrex_real), intent(inout) :: u(domain_lo(1):domain_hi(1), domain_lo(2):domain_hi(2), domain_lo(3):domain_hi(3))
  double precision :: dx, dy, dz
  double precision :: xc, yc, zc  

  integer             :: i,j,k
  
  dx = .1 !(domain_hi(1) - domain_lo(1))/...
  dy = .1
  dz = .1
  
  do k = domain_lo(3), domain_hi(3)
     zc = k*dz 
     do j = domain_lo(2), domain_hi(2)
        yc = j*dy
        do i = domain_lo(1), domain_hi(1) 
           xc = i*dx 
           u(i,j,k) = 0.2*sin(yc) 
  enddo
     enddo
        enddo

end subroutine init_fluid
