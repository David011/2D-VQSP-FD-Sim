!=======================================================================

SUBROUTINE READ_MODEL

  USE PRECISION   , ONLY:                                              &
                          IK
  USE CONTROL_DATA, ONLY:                                              &
                          MX,     MZ, TPML,                            &
                          MO_FILE_NAME, JMH_FILE_NAME,                 &
                          STRESS_IMAGING
                          
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE AUXIL       , ONLY :                                             &
                          PMLL , PMLXH ,         PMLZH,                &
                          X_MIN_D , X_MAX_D,                           &
                          Z_MIN_D , Z_MAX_D

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=3) :: my_N_CH

  INTEGER :: K, L, J
  INTEGER(IK), ALLOCATABLE, DIMENSION(:,:) :: JMT

!-----------------------------------------------------------------------

  ALLOCATE ( JM (X_MIN_D:X_MAX_D,                  Z_MIN_D:Z_MAX_D) )

#if(USE_MPI)
  my_N = ( my_coord(2) )*dims(1) + my_coord(1)
  
  WRITE (my_N_CH, '(I3.3)') my_N
  
  OPEN (13,FILE = TRIM(JMH_FILE_NAME)//my_N_CH, FORM = 'UNFORMATTED',&
                                                      STATUS='OLD' )
  
  DO L = Z_MIN_D, Z_MAX_D
      READ ( 13 )  JM ( :, L)
  END DO
  
  OPEN ( 14, FILE =  TRIM(MO_FILE_NAME)//my_N_CH, FORM = 'UNFORMATTED',&
                                                        STATUS='OLD' )
  
  
#else

  ALLOCATE ( JMT (PMLL:PMLXH,0:PMLZH) )

  IF (STRESS_IMAGING) THEN
    DO L= 2, MZ

      READ  ( 13 )  JMT ( 1:MX, L)

      JMT (:,1) = JMT (:,2)  
      JMT (:,0) = JMT (:,3)
      DO J = 0, TPML
        JMT (1 -J, L) = JMT ( 1, L)
        JMT (MX+J, L) = JMT (MX, L)
      END DO

    END DO
  ELSE
    DO L= 0, MZ

      READ ( 13 )  JMT ( 1:MX, L)

      DO J = 0, TPML
        JMT (1 -J, L) = JMT ( 1, L)
        JMT (MX+J, L) = JMT (MX, L)
      END DO

    END DO
  END IF

  DO L= MZ+1, PMLZH
    JMT (:,L) = JMT (:,MZ)
  END DO

  JM = JMT

  DEALLOCATE(JMT)

#endif

END SUBROUTINE READ_MODEL
