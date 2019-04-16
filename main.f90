module global_variables
  implicit none
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! material parameters
  real(8),parameter :: v_fermi = 1d0
  real(8),parameter :: T12=1d6, T23=1d6, T13=1d6
  real(8),parameter :: kx0_K=0d0,ky0_K=0d0

! system
  integer :: nk,nkz
  real(8) :: rho_11
  real(8),allocatable :: rho_22(:,:),rho_33(:,:)
  complex(8),allocatable :: zrho_12(:,:), zrho_13(:,:), zrho_23(:,:)
  real(8),allocatable :: eps_v(:), eps_c(:), eps_k(:,:)
  complex(8),allocatable :: zdip_v(:,:), zdip_c(:,:), zdip_vc(:,:,:)
  real(8),allocatable :: kx0(:),ky0(:),kz0(:)
  real(8).allocatable :: kx(:),ky(:),kz(:)


! time-propagation
  integer :: nt
  real(8) :: dt, Tprop

! laser
  real(8),allocatable :: At_IR(:,:), At_IR_dt2(:,:)
  real(8),allocatable :: Et_IR(:,:), Et_IR_dt2(:,:)
  real(8),allocatable :: Et_XUV(:), Et_XUV_dt2(:)



end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialize

  call time_propagation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  nk = 128
  nkz = 64


end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none

  rho_11 = 1d0
  allocate(rho_22(nk,nkz), rho_33(nk,nkz))
  rho_22 = 0d0; rho_33 = 0d0
  allocate(zrho_12(nk,nkz), zrho_13(nk,nkz), zrho_23(nk,nkz))
  zrho_12 = 0d0; zrho_12 = 0d0; zrho_23 = 0d0
  allocate(eps_v(nk), eps_c(nk), eps_k(nk,nkz) )
  allocate(zdip_v(nk,nkz), zdip_c(nk,nkz), zdip_vc(2,nk,nkz))
  allocate(kx0(nk),ky0(nk),kz0(nkz))
  allocate(kx(nk),ky(nk),kz(nkz))


  kxy0_v = 0d0 ! has to be K or K' point.



end subroutine initialize
!-------------------------------------------------------------------------------
subroutine drho_dt(rho_22_t, rho_33_t, zrho_12_t, zrho_13_t, zrho_23_t, &
  drho_22_t, drho_33_t, dzrho_12_t, dzrho_13_t, dzrho_23_t, EXUV_t, EIR_t)
  use global_variables
  implicit none
  real(8),intent(in) :: rho_22_t(nk,nkz),rho_33_t(nk,nkz)
  complex(8),intent(in) :: zrho_12_t(nk,nkz), zrho_13_t(nk,nkz), zrho_23_t(nk,nkz)
  real(8),intent(out) :: drho_22_t(nk,nkz), drho_33_t(nk,nkz)
  complex(8),intent(out) :: dzrho_12_t(nk,nkz),dzrho_13_t(nk,nkz),dzrho_23_t(nk,nkz)
  real(8),intent(in) :: EXU_t, EIR_t(2)

  integer :: ik, ikz

  do ik = 1, nk
    do ikz = 1, nkz

      drho_22_t(ik,ikz) = 2d0*real(-zi*(&
        EXUV_t*conjg(zdip_v(ik,ikz))*zrho_12_t(ik,ikz) &
       +(EIR_t(1)*zdip_vc(1,ik,ikz)+EIR_t(2)*zdip_vc(2,ik,ikz))&
       *conjg(zrho_23_t(ik,ikz))&
        ))

      drho_33_t(ik,ikz) = 2d0*real(-zi*(&
        EXUV_t*conjg(zdip_c(ik,ikz))*zrho_13_t(ik,ikz) &
       +conjg(EIR_t(1)*zdip_vc(1,ik,ikz)+EIR_t(2)*zdip_vc(2,ik,ikz))&
       *zrho_23_t(ik,ikz)&
       ))

      dzrho_12_t(ik,ikz) = zi*(eps_k(ik,ikz)-eps_v(ik))*zrho_12_t(ik,ikz) &
        -zi*EXUV_t*zdip_v(ik,ikz)* 1d0 &
        -zi*EXUV_t*zdip_c(ik,ikz)*conjg(zrho_23_t(ik,ikz)) &
        +zi*conjg(EIR_t(1)*zdip_vc(1,ik,ikz)+EIR_t(2)*zdip_vc(2,ik,ikz))&
        *zrho_13_t(ik,ikz) -zrho_12_t(ik,ikz)/T12

      dzrho_13_t(ik,ikz) = zi*(eps_k(ik,ikz)-eps_c(ik))*zrho_13_t(ik,ikz) &
        -zi*EXUV_t*zdip_c(ik,ikz)* 1d0 &
        -zi*EXUV_t*zdip_v(ik,ikz)*zrho_23_t(ik,ikz) &
        +zi*(EIR_t(1)*zdip_vc(1,ik,ikz)+EIR_t(2)*zdip_vc(2,ik,ikz))&
        *zrho_12_t(ik,ikz) -zrho_13_t(ik,ikz)/T13

      dzrho_23_t(ik,ikz) = -zi*(eps_c(ik)-eps_v(ik))*zrho_23_t(ik,ikz) &
        -zi*(EIR_t(1)*zdip_vc(1,ik,ikz)+EIR_t(2)*zdip_vc(2,ik,ikz)) &
        *(rho_33_t(ik,ikz)-rho_22_t(ik,ikz)) &
        -zi*EXUV_t*conjg(zdip_v(ik,ikz))*zrho_13_t(ik,ikz) &
        +zi*EXUV_t*zdip_c(ik,ikz)*conjg(zrho_12_t(ik,ikz)) -zrho_23_t(ik,ikz)/T23

    end do
  end do

end subroutine drho_dt
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it

  call init_laser

  do it = 0, nt

    call dt_evolve(it)

  end do



end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  real(8) :: At_IR_t(2), Et_IR_t(2), Et_XUV_t

! time, t
  At_IR_t(:) = At_IR(:,it)
  Et_IR_t(:) = Et_IR(:,it)
  Et_XUV_t = Et_XUV(it)

  call set_matrix_elements(At_IR_t)



end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine set_matrix_elements(At_IR_t)
  use global_variables
  implicit none
  real(8),intent(in) :: At_IR_t(2)
  integer :: ik, ikz

  kx(:) = kx0(:) + At_IR_t(1)
  ky(:) = ky0(:) + At_IR_t(2)
  kz = kz0

  do ik = 1, nk
    eps_c(ik) = v_fermi*sqrt(kx(ik)**2+ky(ik)**2)
    eps_v(ik) = -eps_c(ik)
  end do

  do ik = 1, nk
    do ikz = 1, nkz
      eps_k(ik,ikz) = 0.5d0*(&
        (kx(ik)+kx0_K)**2+(ky(ik)+ky0_K)**2+kz(ikz)**2)
    end do
  end do



end subroutine set_matrix_elements
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
