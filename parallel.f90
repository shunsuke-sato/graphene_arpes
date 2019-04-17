module parallel
  use mpi
  implicit none

  private
! MPI global
  integer, public :: comm_group_global, &
                     comm_id_global, &
                     comm_nproc_global
  logical, public :: if_root_global
                     
! OMP
  integer, public :: nthread_omp

! temporal arrays
  character(len=1024),public :: message(64)

  public :: init_parallel, &
            fin_parallel,  &
            error_finalize, &
            write_message

  interface write_message
     module procedure write_message_array
     module procedure write_message_scalar
  end interface write_message

  interface error_finalize
     module procedure error_finalize_array
     module procedure error_finalize_scalar
  end interface error_finalize

contains
!-------------------------------------------------------------------------------
  subroutine init_parallel
    implicit none
    integer :: ierr
!$ integer :: omp_get_max_threads  

    call MPI_init(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,comm_nproc_global,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,comm_id_global,ierr)

    comm_group_global = MPI_COMM_WORLD

    if(comm_id_global == 0)then
       if_root_global = .true.
    else
       if_root_global = .false.
    end if

    nthread_omp = 1
!$  nthread_omp=omp_get_max_threads()

  end subroutine init_parallel
!-------------------------------------------------------------------------------
  subroutine fin_parallel
    implicit none
    integer :: ierr

    call MPI_Finalize(ierr)

  end subroutine fin_parallel
!-------------------------------------------------------------------------------
  subroutine error_finalize_array(message_t)
    implicit none
    character(*),intent(in) :: message_t(:)
    integer :: ierr
    integer :: ndim, i

    if(if_root_global)then
      ndim = ubound(message_t, 1)
      do i = 1,ndim
        write(*,"(A)")trim(message_t(i))
      end do
    end if

    call MPI_Finalize(ierr)
    stop

  end subroutine error_finalize_array
!-------------------------------------------------------------------------------
  subroutine error_finalize_scalar(message_t)
    implicit none
    character(*),intent(in) :: message_t
    integer :: ierr

    if(if_root_global)write(*,"(A)")trim(message_t)
    call MPI_Finalize(ierr)
    stop

  end subroutine error_finalize_scalar
!-------------------------------------------------------------------------------
  subroutine write_message_array(message_t)
    implicit none
    character(*),intent(in) :: message_t(:)
    integer :: ndim, i

    if(if_root_global)then
      ndim = ubound(message_t, 1)
      do i = 1,ndim
        write(*,"(A)")trim(message_t(i))
      end do
    end if

  end subroutine write_message_array
!-------------------------------------------------------------------------------
  subroutine write_message_scalar(message_t)
    implicit none
    character(*),intent(in) :: message_t

    if(if_root_global)then
        write(*,"(A)")trim(message_t)
    end if

  end subroutine write_message_scalar
!-------------------------------------------------------------------------------
end module parallel
