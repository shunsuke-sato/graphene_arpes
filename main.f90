module global_variables
  implicit none
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)


! system
  integer :: nk,nkz
  real(8) :: rho_11
  real(8),allocatable :: rho_22(:,:),rho_33(:,:)
  complex(8),allocatable :: zrho_12(:,:), zrho_13(:,:), zrho_23(:,:)
  real(8),allocatable :: eps_v(:,:), eps_c(:,:), eps_k(:,:)
  complex(8),allocatable :: zdip_v(:,:), zdip_c(:,:), zdip_vc(:,:)
  real(8) :: T12, T23, T13
  real(8) :: kx0(:),ky0(:),kz0(:)
  real(8) :: kx(:),ky(:),kz(:)




end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none
end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none

end subroutine initialize
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
