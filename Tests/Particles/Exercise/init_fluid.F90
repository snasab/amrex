subroutine ux_init(domain_lo, domain_hi, ux, ux_lo, ux_hi) bind(C, name="ux_init")

  use amrex_fort_module, only: amrex_real

  implicit none
  
  integer, intent(in)    :: domain_lo(3), domain_hi(3), ux_lo(3), ux_hi(3)  
  real(amrex_real), intent(inout) :: ux(ux_lo(1):ux_hi(1), ux_lo(2):ux_hi(2), ux_lo(3):ux_hi(3))
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
           ux(i,j,k) = 5. 
  enddo
     enddo
        enddo

end subroutine ux_init

subroutine uy_init(domain_lo, domain_hi, uy, uy_lo, uy_hi) bind(C, name="uy_init")

  use amrex_fort_module, only: amrex_real

  implicit none

  integer, intent(in)    :: domain_lo(3), domain_hi(3), uy_lo(3), uy_hi(3)
  real(amrex_real), intent(inout) :: uy(uy_lo(1):uy_hi(1), uy_lo(2):uy_hi(2), uy_lo(3):uy_hi(3))
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
           uy(i,j,k) = 0.
  enddo
     enddo
        enddo

end subroutine uy_init

subroutine uz_init(domain_lo, domain_hi, uz, uz_lo, uz_hi) bind(C, name="uz_init")

  use amrex_fort_module, only: amrex_real

  implicit none

  integer, intent(in)    :: domain_lo(3), domain_hi(3), uz_lo(3), uz_hi(3)
  real(amrex_real), intent(inout) :: uz(uz_lo(1):uz_hi(1), uz_lo(2):uz_hi(2), uz_lo(3):uz_hi(3))
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
           uz(i,j,k) = 0.
  enddo
     enddo
        enddo

end subroutine uz_init


subroutine amrex_init_fluid(domain_lo, domain_hi, fluid, f_lo, f_hi, ncU) bind(C, name="amrex_init_fluid") 
  use amrex_fort_module, only: amrex_real

  implicit none
  integer, intent(in) :: domain_lo(3), domain_hi(3), f_lo(3), f_hi(3), ncU 
  real(amrex_real), intent(inout) :: fluid(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), ncU)
  integer :: i,j,k,n
   
  do n= 1, ncU 
    do     k = domain_lo(3), domain_hi(3)
      do   j = domain_lo(2), domain_hi(2)
        do i = domain_lo(1), domain_hi(1)
          fluid(i,j,k,ncU) = 1.0
        end do
      end do
    end do
  enddo

  

end subroutine amrex_init_fluid




