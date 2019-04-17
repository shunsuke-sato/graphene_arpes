!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer(a_in, a_out, communicator, method)
  implicit none
  integer,intent(inout)       :: a_in
  integer,intent(in),optional :: a_out
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer :: a_t
  integer :: id_comm, method_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, 1, MPI_INTEGER, method_t, id_comm, ierr)
  else
    call MPI_Allreduce(a_in, a_t, 1, MPI_INTEGER, method_t, id_comm, ierr)
    a_in = a_t
  end if

end subroutine comm_allreduce_integer
!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer_1d(a_in, a_out, communicator, method, nsize)
  implicit none
  integer,intent(inout)       :: a_in(:)
  integer,intent(in),optional :: a_out(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer,allocatable :: a_t(:)
  integer :: lb1, ub1
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    allocate(a_t(lb1:ub1))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_integer_1d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer_2d(a_in, a_out, communicator, method, nsize)
  implicit none
  integer,intent(inout)       :: a_in(:,:)
  integer,intent(in),optional :: a_out(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer,allocatable :: a_t(:,:)
  integer :: lb1, ub1, lb2, ub2
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    lb2 = lbound(a_in,2); ub2 = ubound(a_in,2)
    allocate(a_t(lb1:ub1,lb2:ub2))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_integer_2d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_integer_3d(a_in, a_out, communicator, method, nsize)
  implicit none
  integer,intent(inout)       :: a_in(:,:,:)
  integer,intent(in),optional :: a_out(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  integer,allocatable :: a_t(:,:,:)
  integer :: lb1, ub1, lb2, ub2, lb3, ub3
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    lb2 = lbound(a_in,2); ub2 = ubound(a_in,2)
    lb3 = lbound(a_in,3); ub3 = ubound(a_in,3)
    allocate(a_t(lb1:ub1,lb2:ub2,lb3:ub3))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_INTEGER, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_integer_3d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8(a_in, a_out, communicator, method)
  implicit none
  real(8),intent(inout)       :: a_in
  real(8),intent(in),optional :: a_out
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  real(8) :: a_t
  integer :: id_comm, method_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, 1, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
  else
    call MPI_Allreduce(a_in, a_t, 1, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
    a_in = a_t
  end if

end subroutine comm_allreduce_real8
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8_1d(a_in, a_out, communicator, method, nsize)
  implicit none
  real(8),intent(inout)       :: a_in(:)
  real(8),intent(in),optional :: a_out(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  real(8),allocatable :: a_t(:)
  integer :: lb1, ub1
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    allocate(a_t(lb1:ub1))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_real8_1d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8_2d(a_in, a_out, communicator, method, nsize)
  implicit none
  real(8),intent(inout)       :: a_in(:,:)
  real(8),intent(in),optional :: a_out(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  real(8),allocatable :: a_t(:,:)
  integer :: lb1, ub1, lb2, ub2
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    lb2 = lbound(a_in,2); ub2 = ubound(a_in,2)
    allocate(a_t(lb1:ub1,lb2:ub2))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_real8_2d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_real8_3d(a_in, a_out, communicator, method, nsize)
  implicit none
  real(8),intent(inout)       :: a_in(:,:,:)
  real(8),intent(in),optional :: a_out(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  real(8),allocatable :: a_t(:,:,:)
  integer :: lb1, ub1, lb2, ub2, lb3, ub3
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    lb2 = lbound(a_in,2); ub2 = ubound(a_in,2)
    lb3 = lbound(a_in,3); ub3 = ubound(a_in,3)
    allocate(a_t(lb1:ub1,lb2:ub2,lb3:ub3))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_DOUBLE_PRECISION, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_real8_3d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8(a_in, a_out, communicator, method)
  implicit none
  complex(8),intent(inout)       :: a_in
  complex(8),intent(in),optional :: a_out
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  complex(8) :: a_t
  integer :: id_comm, method_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, 1, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
  else
    call MPI_Allreduce(a_in, a_t, 1, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
    a_in = a_t
  end if

end subroutine comm_allreduce_complex8
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8_1d(a_in, a_out, communicator, method, nsize)
  implicit none
  complex(8),intent(inout)       :: a_in(:)
  complex(8),intent(in),optional :: a_out(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  complex(8),allocatable :: a_t(:)
  integer :: lb1, ub1
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    allocate(a_t(lb1:ub1))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_complex8_1d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8_2d(a_in, a_out, communicator, method, nsize)
  implicit none
  complex(8),intent(inout)       :: a_in(:,:)
  complex(8),intent(in),optional :: a_out(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  complex(8),allocatable :: a_t(:,:)
  integer :: lb1, ub1, lb2, ub2
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    lb2 = lbound(a_in,2); ub2 = ubound(a_in,2)
    allocate(a_t(lb1:ub1,lb2:ub2))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_complex8_2d
!-------------------------------------------------------------------------------
subroutine comm_allreduce_complex8_3d(a_in, a_out, communicator, method, nsize)
  implicit none
  complex(8),intent(inout)       :: a_in(:,:,:)
  complex(8),intent(in),optional :: a_out(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: method
  integer,intent(in),optional :: nsize
  complex(8),allocatable :: a_t(:,:,:)
  integer :: lb1, ub1, lb2, ub2, lb3, ub3
  integer :: id_comm, method_t, nsize_t, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)

  method_t = MPI_SUM
  if(present(method)) method_t = method

  nsize_t = size(a_in)
  if(present(nsize)) nsize_t = nsize

  if(present(a_out))then
    call MPI_Allreduce(a_in, a_out, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
  else
    lb1 = lbound(a_in,1); ub1 = ubound(a_in,1)
    lb2 = lbound(a_in,2); ub2 = ubound(a_in,2)
    lb3 = lbound(a_in,3); ub3 = ubound(a_in,3)
    allocate(a_t(lb1:ub1,lb2:ub2,lb3:ub3))
    call MPI_Allreduce(a_in, a_t, nsize_t, MPI_DOUBLE_COMPLEX, method_t, id_comm, ierr)
    a_in = a_t
    deallocate(a_t)
  end if

end subroutine comm_allreduce_complex8_3d


