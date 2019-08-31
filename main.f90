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
  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: angstrom = 1d0/0.52917721067d0
  real(8),parameter :: clight = 137.035999139d0

! material parameters
  real(8),parameter :: v_fermi = clight*1.12d6/299792458d0
  real(8),parameter :: T12=20d0*fs, T23=20d0*fs, T13=20d0*fs
  real(8),parameter :: kx0_K=4d0*pi/(3d0*2.46d0*angstrom),ky0_K=0d0
  real(8),parameter :: work_function = 4.5d0*ev

! system
  integer :: nkxy,nkz,nk
  integer :: nk_start, nk_end
  integer :: nk_average, nk_remainder
  integer,allocatable :: ikxy_table(:),ikz_table(:)
  complex(8),allocatable :: zrho_k(:,:,:)
  complex(8),allocatable :: zdip_B(:), zdip_A(:)
  real(8),allocatable :: kx0(:),ky0(:),kz0(:)
  real(8),allocatable :: kx(:),ky(:),kz(:)


! time-propagation
  integer :: nt
  real(8) :: dt, Tprop

! laser
  real(8) :: E0_IR, omega0_IR, tpulse_IR
  real(8) :: E0_XUV, omega0_XUV, tpulse_XUV
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

  Tprop = 380d0*fs*2.5d0
  dt = 0.1d0
  nt = aint(Tprop/dt)+1

  E0_IR = 0d0*ev/angstrom
  omega0_IR = 280d-3*ev
  tpulse_IR = 380d0*fs/0.364056664d0

  E0_XUV = 1d-6
  omega0_XUV = 22d0*ev
  tpulse_XUV = 120d0*fs/0.364056664d0


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
  allocate(zdip_B(nk), zdip_A(nk))
  allocate(kx0(nk),ky0(nk),kz0(nk))
  allocate(kx(nk),ky(nk),kz(nk))
  allocate(ikxy_table(nk),ikz_table(nk))

  zrho_k = 0d0
  zrho_k(1,1,:) = 1d0

  call init_k_grids
  call init_laser



end subroutine initialize
!-------------------------------------------------------------------------------
subroutine init_k_grids
  use global_variables
  implicit none
  integer :: ik, ikx, ikz
  real(8) :: kx_max,kx_min, dkx
  real(8) :: kz_max,kz_min, dkz

  kx_min = -0.5d0*0.5d0*ev/v_fermi
  kx_max =  0.5d0*0.5d0*ev/v_fermi
  dkx = (kx_max-kx_min)/(nkxy-1)

  ky0 = 0d0

  ik = 0
  do ikx = 1, nkxy
    do ikz = 1, nkz
      ik = ik + 1
      kx0(ik) = kx_min + dkx*(ikx-1)

      kz_min = sqrt(2d0*(omega0_XUV-0.5d0*((kx0_K+kx0(ik))**2+ky0_K**2)-work_function-0.5d0*ev))
      kz_max = sqrt(2d0*(omega0_XUV-0.5d0*((kx0_K+kx0(ik))**2+ky0_K**2)-work_function+0.5d0*ev))
      dkz = (kz_max-kz_min)/(nkz-1)


      kz0(ik) = kz_min + dkz*(ikz-1)
      ikxy_table(ik) = ikx
      ikz_table(ik)  = ikz
    end do
  end do


end subroutine init_k_grids
!-------------------------------------------------------------------------------
subroutine init_laser
  use global_variables
  implicit none
  real(8) :: tt, xx, tcenter
  integer :: it

  allocate(At_IR(2,0:nt+1), At_IR_dt2(2,0:nt+1))
  allocate(Et_XUV(0:nt+1), Et_XUV_dt2(0:nt+1))
  At_IR = 0d0; At_IR_dt2 = 0d0
  Et_XUV = 0d0; Et_XUV_dt2 = 0d0

  tcenter = max(0.5d0*tpulse_IR,0.5d0*tpulse_XUV)
  

  do it = 0, nt+1

    tt = dt*it
    xx = tt - tcenter
    if(abs(xx) < 0.5d0*tpulse_IR)then
      At_IR(1,it) = -E0_IR/omega0_IR*cos(pi*xx/tpulse_IR)**2*sin(omega0_IR*xx)
      At_IR(2,it) = 0d0
    end if

    if(abs(xx) < 0.5d0*tpulse_XUV)then
      Et_XUV(it) = -E0_XUV*cos(pi*xx/tpulse_XUV)**2*sin(omega0_XUV*xx)
    end if

    tt = dt*it + 0.5d0*dt
    xx = tt - tcenter
    if(abs(xx) < 0.5d0*tpulse_IR)then
      At_IR_dt2(1,it) = -E0_IR/omega0_IR*cos(pi*xx/tpulse_IR)**2*sin(omega0_IR*xx)
      At_IR_dt2(2,it) = 0d0
    end if

    if(abs(xx) < 0.5d0*tpulse_XUV)then
      Et_XUV_dt2(it) = -E0_XUV*cos(pi*xx/tpulse_XUV)**2*sin(omega0_XUV*xx)
    end if

  end do

  if(if_root_global)then
    open(20,file='laser.out')
    do it = 0, nt + 1
      write(20,"(999e16.6e3)")dt*it,At_IR(:,it),Et_XUV(it)
    end do
    close(20)
  end if


end subroutine init_laser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it, ik, ikx, ikz
  real(8) :: pop_k(nkxy,nkz),eps_k

  do it = 0, nt
    if(if_root_global .and. (mod(it,nt/200)==0))write(*,*)'it',it,nt
    call dt_evolve(it)

  end do

  pop_k = 0d0
  do ik = nk_start, nk_end
    pop_k(ikxy_table(ik),ikz_table(ik))=real(zrho_k(2,2,ik) + zrho_k(3,3,ik))
  end do

  call comm_allreduce(pop_k)
  if(if_root_global)then
    open(20,file='pes_spectrum.out')
    ik = 0
    do ikx = 1, nkxy
      do ikz = 1,nkz
        ik = ik + 1
        eps_k = 0.5d0*( (kx0(ik)+kx0_K)**2 + (ky0(ik)+ky0_K)**2 + kz0(ik)**2) + work_function
        write(20,"(999e26.16e3)")kx0(ik),ky0(ik),eps_k,pop_k(ikxy_table(ik),ikz_table(ik))
      end do
      write(20,*)
    end do
    close(20)
  end if

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik

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

  zdip_B(ik)= -1d0
  zdip_A(ik)=  1d0

! time, t
  At_IR_t(:) = At_IR(:,it)
  Et_XUV_t = Et_XUV(it)

  kxt = kx0(ik) + At_IR_t(1)
  kyt = ky0(ik) + At_IR_t(2)
  kzt = kz0(ik)

  eps_k = 0.5d0*( (kxt+kx0_K)**2 + (kyt+ky0_K)**2 + kzt**2) + work_function


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

  eps_k = 0.5d0*( (kxt+kx0_K)**2 + (kyt+ky0_K)**2 + kzt**2) + work_function


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

  eps_k = 0.5d0*( (kxt+kx0_K)**2 + (kyt+ky0_K)**2 + kzt**2) + work_function


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

  
  zrho_k(:,:,ik) = zrho_k(:,:,ik) + dt/6d0*(zLrho_rk(:,:,1)   &
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
  complex(8) :: zrho_col(3,3)

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
  zrho_col(1,2) = conjg(zrho_col(2,1))
  zrho_col(2,2) = 0d0
  zrho_col(3,2) = -zrho_col(3,2)/T23
  zrho_col(1,3) = conjg(zrho_col(3,1))
  zrho_col(2,3) = conjg(zrho_col(3,2))
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
