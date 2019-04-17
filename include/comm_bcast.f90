!-------------------------------------------------------------------------------
subroutine comm_bcast_integer(a, communicator, root)
  implicit none
  integer,intent(inout) :: a
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, 1, MPI_INTEGER, id_root, id_comm, ierr)

end subroutine comm_bcast_integer
!-------------------------------------------------------------------------------
subroutine comm_bcast_integer_1d(a, communicator, root)
  implicit none
  integer,intent(inout) :: a(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, size(a), MPI_INTEGER, id_root, id_comm, ierr)

end subroutine comm_bcast_integer_1d
!-------------------------------------------------------------------------------
subroutine comm_bcast_integer_2d(a, communicator, root)
  implicit none
  integer,intent(inout) :: a(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_INTEGER, id_root, id_comm, ierr)
  
end subroutine comm_bcast_integer_2d
!-------------------------------------------------------------------------------
subroutine comm_bcast_integer_3d(a, communicator, root)
  implicit none
  integer,intent(inout) :: a(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_INTEGER, id_root, id_comm, ierr)

end subroutine comm_bcast_integer_3d
!-------------------------------------------------------------------------------
subroutine comm_bcast_real8(a, communicator, root)
  implicit none
  real(8),intent(inout) :: a
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, id_root, id_comm, ierr)

end subroutine comm_bcast_real8
!-------------------------------------------------------------------------------
subroutine comm_bcast_real8_1d(a, communicator, root)
  implicit none
  real(8),intent(inout) :: a(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, id_root, id_comm, ierr)
  
end subroutine comm_bcast_real8_1d
!-------------------------------------------------------------------------------
subroutine comm_bcast_real8_2d(a, communicator, root)
  implicit none
  real(8),intent(inout) :: a(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, id_root, id_comm, ierr)

end subroutine comm_bcast_real8_2d
!-------------------------------------------------------------------------------
subroutine comm_bcast_real8_3d(a, communicator, root)
  implicit none
  real(8),intent(inout) :: a(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, id_root, id_comm, ierr)
  
end subroutine comm_bcast_real8_3d
!-------------------------------------------------------------------------------
subroutine comm_bcast_real8_4d(a, communicator, root)
  implicit none
  real(8),intent(inout) :: a(:,:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, id_root, id_comm, ierr)
  
end subroutine comm_bcast_real8_4d
!-------------------------------------------------------------------------------
subroutine comm_bcast_complex8(a, communicator, root)
  implicit none
  complex(8),intent(inout) :: a
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, 1, MPI_DOUBLE_COMPLEX, id_root, id_comm, ierr)
  
end subroutine comm_bcast_complex8
!-------------------------------------------------------------------------------
subroutine comm_bcast_complex8_1d(a, communicator, root)
  implicit none
  complex(8),intent(inout) :: a(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)  
  call MPI_Bcast(a, size(a), MPI_DOUBLE_COMPLEX, id_root, id_comm, ierr)
  
end subroutine comm_bcast_complex8_1d
!-------------------------------------------------------------------------------
subroutine comm_bcast_complex8_2d(a, communicator, root)
  implicit none
  complex(8),intent(inout) :: a(:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, size(a), MPI_DOUBLE_COMPLEX, id_root, id_comm, ierr)
  
end subroutine comm_bcast_complex8_2d
!-------------------------------------------------------------------------------
subroutine comm_bcast_complex8_3d(a, communicator, root)
  implicit none
  complex(8),intent(inout) :: a(:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, size(a), MPI_DOUBLE_COMPLEX, id_root, id_comm, ierr)
  
end subroutine comm_bcast_complex8_3d
!-------------------------------------------------------------------------------
subroutine comm_bcast_complex8_4d(a, communicator, root)
  implicit none
  complex(8),intent(inout) :: a(:,:,:,:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, size(a), MPI_DOUBLE_COMPLEX, id_root, id_comm, ierr)
  
end subroutine comm_bcast_complex8_4d
!-------------------------------------------------------------------------------
subroutine comm_bcast_character(a, communicator, root)
  implicit none
  character(*),intent(inout) :: a
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, len(a), MPI_CHARACTER, id_root, id_comm, ierr)
  
end subroutine comm_bcast_character
!-------------------------------------------------------------------------------
subroutine comm_bcast_character_1d(a, communicator, root)
  implicit none
  character(*),intent(inout) :: a(:)
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, size(a)*len(a), MPI_CHARACTER, id_root, id_comm, ierr)
  
end subroutine comm_bcast_character_1d
!-------------------------------------------------------------------------------
subroutine comm_bcast_logical(a, communicator, root)
  implicit none
  logical,intent(inout) :: a
  integer,intent(in),optional :: communicator
  integer,intent(in),optional :: root
  integer :: id_comm, id_root, ierr

  id_comm = int_switch(present(communicator), communicator, comm_group_global)
  id_root = int_switch(present(root), root, 0)
  call MPI_Bcast(a, 1, MPI_LOGICAL, id_root, id_comm, ierr)
  
end subroutine comm_bcast_logical
!-------------------------------------------------------------------------------
