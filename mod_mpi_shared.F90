MODULE MPI_SHARED

#if(USE_MPI)

#if(MPI2)
  USE MPI
#else
  INCLUDE 'mpif.h'
#endif
!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  INTEGER, PARAMETER :: ndims = 2

  CHARACTER(LEN=40)  :: nam_proc

  LOGICAL        :: iamroot, periods(ndims)
  INTEGER        :: my_N
  INTEGER        :: comm, comm_cart, errcode, myid, root, numprocs,    &
                    status(MPI_STATUS_SIZE), request

  INTEGER        :: len, tag, tagq, OVERLAP

  INTEGER, DIMENSION(ndims)     :: dims, my_coord, coord

  INTEGER, ALLOCATABLE, DIMENSION(:) :: X_MIN_DA , X_MAX_DA ,          &
                                        Z_MIN_DA , Z_MAX_DA ,          &
                                        X_MIN_A  , X_MAX_A  ,          &
                                        Z_MIN_A  , Z_MAX_A

#ENDIF

END MODULE MPI_SHARED