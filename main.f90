include "parallel.f90"
include "communication.f90"


module global_variables
  use mpi
  use parallel
  use communication
  implicit none
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physics parameters
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0

! material parameters
  real(8),parameter :: v_fermi = 1d0
  real(8),parameter :: T12=1d6, T23=1d6, T13=1d6
  real(8),parameter :: kx0_K=0d0,ky0_K=0d0

! system
  integer :: nkxy,nkz,nk
  integer :: nk_start, nk_end
  integer :: nk_average, nk_remainder
  complex(8),allocatable :: zrho_k(:,:,:)
  real(8),allocatable :: eps_k(:)
  complex(8),allocatable :: zdip_B(:), zdip_A(:)
  real(8),allocatable :: kx0(:),ky0(:),kz0(:)
  real(8).allocatable :: kx(:),ky(:),kz(:)


! time-propagation
  integer :: nt
  real(8) :: dt, Tprop

! laser
  real(8),allocatable :: At_IR(:,:), At_IR_dt2(:,:)
  real(8),allocatable :: Et_XUV(:), Et_XUV_dt2(:)



end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call init_parallel

  call input
  call initialize

  call time_propagation


  call fin_parallel

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  nkxy = 64
  nkz  = 64
  nk = nkxy*nkz

  Tprop = 40d0*fs
  dt = 0.01d0
  nt = aint(Tprop/dt)+1



end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none

  if(nk < comm_nproc_global) call error_finalize('Error: nk < # of MPI processes')
  nk_average = nk/comm_nproc_global
  nk_remainder = mod(nk,comm_nproc_global)
  if(comm_id_global+1 <= nk_remainder)then
    nk_start = 1 + comm_id_global*(nk_average+1)
    nk_end   = nk_start + (nk_average + 1) -1
  else
    nk_start = 1 + nk_remainder*(nk_average+1) + nk_average*(comm_id_global - nk_remainder)
    nk_end    = nk_start + nk_average  -1
  end if
  

  allocate(zrho_k(3,3,nk_start:nk_end))
  allocate(zdip_v(nk), zdip_c(nk))
  allocate(kx0(nk),ky0(nk),kz0(nkz))
  allocate(kx(nk),ky(nk),kz(nkz))

  zrho_k = 0d0
  zrho_k(1,1,:) = 1d0

  call init_k_grids
  call init_laser



end subroutine initialize
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it

  do it = 0, nt

    call dt_evolve(it)

  end do



end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it

  do ik = nk_start, nk_end
    call dt_evolve_k(it,ik)
  end do


end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine dt_evolve_k(it,ik)
  use global_variables
  implicit none
  integer,intent(in) :: it, ik
  real(8) :: kxt, kyt, kzt, eps_k
  real(8) :: At_IR_t(2), Et_XUV_t
  complex(8) :: zrho_t(3,3),zham(3,3), zLrho_rk(3,3,4)

! time, t
  At_IR_t(:) = At_IR(:,it)
  Et_XUV_t = Et_XUV(it)

  kxt = kx0(ik) + At_IR_t(1)
  kyt = ky0(ik) + At_IR_t(2)
  kzt = kz0(ik)

  eps_k = 0.5d0*( (kxt+kx0_K)**2 + (kyt+ky0_K)**2 + kzt**2)


  zham(1,1) = 0d0
  zham(2,1) = zdip_B(ik)*Et_XUV_t
  zham(3,1) = zdip_A(ik)*Et_XUV_t
  zham(1,2) = conjg(zham(2,1))
  zham(2,2) = eps_k
  zham(3,2) = v_fermi*(kxt+zi*kyt)
  zham(1,3) = conjg(zham(3,1))
  zham(2,3) = conjg(zham(3,2))
  zham(3,3) = eps_k

! RK1
  zrho_t(:,:) = zrho_k(:,:,ik)
  call zLrho_op(zrho_t,zham,zLrho_rk(:,:,1))


! time, t+dt/2
  At_IR_t(:) = At_IR_dt2(:,it)
  Et_XUV_t = Et_XUV_dt2(it)

  kxt = kx0(ik) + At_IR_t(1)
  kyt = ky0(ik) + At_IR_t(2)
  kzt = kz0(ik)

  eps_k = 0.5d0*( (kxt+kx0_K)**2 + (kyt+ky0_K)**2 + kzt**2)


  zham(1,1) = 0d0
  zham(2,1) = zdip_B(ik)*Et_XUV_t
  zham(3,1) = zdip_A(ik)*Et_XUV_t
  zham(1,2) = conjg(zham(2,1))
  zham(2,2) = eps_k
  zham(3,2) = v_fermi*(kxt+zi*kyt)
  zham(1,3) = conjg(zham(3,1))
  zham(2,3) = conjg(zham(3,2))
  zham(3,3) = eps_k


! RK2
  zrho_t(:,:) = zrho_k(:,:,ik) + 0.5d0*dt*zLrho_rk(:,:,1)
  call zLrho_op(zrho_t,zham,zLrho_rk(:,:,2))

! RK3
  zrho_t(:,:) = zrho_k(:,:,ik) + 0.5d0*dt*zLrho_rk(:,:,2)
  call zLrho_op(zrho_t,zham,zLrho_rk(:,:,3))

! time, t+dt
  At_IR_t(:) = At_IR(:,it+1)
  Et_XUV_t = Et_XUV(it+1)

  kxt = kx0(ik) + At_IR_t(1)
  kyt = ky0(ik) + At_IR_t(2)
  kzt = kz0(ik)

  eps_k = 0.5d0*( (kxt+kx0_K)**2 + (kyt+ky0_K)**2 + kzt**2)


  zham(1,1) = 0d0
  zham(2,1) = zdip_B(ik)*Et_XUV_t
  zham(3,1) = zdip_A(ik)*Et_XUV_t
  zham(1,2) = conjg(zham(2,1))
  zham(2,2) = eps_k
  zham(3,2) = v_fermi*(kxt+zi*kyt)
  zham(1,3) = conjg(zham(3,1))
  zham(2,3) = conjg(zham(3,2))
  zham(3,3) = eps_k

! RK4
  zrho_t(:,:) = zrho_k(:,:,ik) + dt*zLrho_rk(:,:,3)
  call zLrho_op(zrho_t,zham,zLrho_rk(:,:,4))

  
  zrho_k(:,:,ik) + zrho_k(:,:,ik) + dt/6d0*(zLrho_rk(:,:,1)   &
                                       +2d0*zLrho_rk(:,:,2)   &
                                       +2d0*zLrho_rk(:,:,3)   &
                                       +    zLrho_rk(:,:,4))


end subroutine dt_evolve_k
!-------------------------------------------------------------------------------
subroutine zLrho_op(zrho_in,zham,zLrho_out)
  use global_variables
  implicit none
  complex(8),intent(in) :: zrho_in(3,3), zham(3,3)
  complex(8),intent(out) :: zLrho_out(3,3)
  real(8) :: alpha_r, alpha_i, phi
  complex(8) :: zalpha, zeig(3,3)
  complex(8) :: zrho_col

  zLrho_out = -zi*(matmul(zham,zrho_in)-matmul(zrho_in,zham))

  alpha_r =  real(zham(3,2))
  alpha_i = aimag(zham(3,2))

  if(alpha_r > 0)then
    phi = atan(alpha_i/alpha_r)
  else if(alpha_r <0)then
    phi = atan(alpha_i/alpha_r) + pi
  else if(alpha_i >0)then
    phi = 0.5d0*pi
  else
    phi = -0.5d0*pi
  end if
  zalpha = exp(zi*phi)

  zeig(1,1) = 1d0
  zeig(2,1) = 0d0
  zeig(3,1) = 0d0

  zeig(1,2) = 0d0
  zeig(2,2) = 1d0/sqrt(2d0)
  zeig(3,2) = -zalpha/sqrt(2d0)

  zeig(1,3) = 0d0
  zeig(2,3) = 1d0/sqrt(2d0)
  zeig(3,3) = zalpha/sqrt(2d0)

  zrho_col = matmul(transpose(conjg(zeig)), matmul(zrho_in,zeig))
  zrho_col(1,1) = 0d0
  zrho_col(2,1) = -zrho_col(2,1)/T12
  zrho_col(3,1) = -zrho_col(3,1)/T13
  zrho_col(1,2) = conjg(zrho_col(1,2))
  zrho_col(2,2) = 0d0
  zrho_col(3,2) = -zrho_col(3,2)/T23
  zrho_col(1,3) = conjg(zrho_col(3,1))
  zrho_col(2,3) = conjg(zrho_col(2,3))
  zrho_col(3,3) = 0d0

  zrho_col = matmul( matmul(zeig, zrho_col), transpose(conjg(zeig)))

  
  zLrho_out = zLrho_out + zrho_col 

  zLrho_out(1,1) = 0d0

end subroutine zLrho_op
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
