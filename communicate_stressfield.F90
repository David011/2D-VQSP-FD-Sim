SUBROUTINE COMMUNICATE_STRESSFIELD

  USE STRESSFIELD , ONLY: TXX , TXZ , TZZ, PRES

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  INTERFACE
    SUBROUTINE     COMMUNICATE (VARIABLE)
    USE PRECISION   , ONLY : PP
    USE AUXIL       , ONLY : X_MIN_D, X_MAX_D,                         &
                             Z_MIN_D, Z_MAX_D
    IMPLICIT NONE
    REAL(PP), INTENT(INOUT),                                           &
      DIMENSION(X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) ::  VARIABLE
    END SUBROUTINE COMMUNICATE

    SUBROUTINE     COMMUNICATE_PLANE (VARIABLE)
    USE PRECISION   , ONLY : PP
    USE AUXIL       , ONLY : X_MIN_D, X_MAX_D
    IMPLICIT NONE
    REAL(PP), INTENT(INOUT),                                           &
               DIMENSION(X_MIN_D:X_MAX_D) :: VARIABLE
    END SUBROUTINE COMMUNICATE_PLANE

  END INTERFACE

!-----------------------------------------------------------------------

  CALL COMMUNICATE   ( TXX  )
  CALL COMMUNICATE   ( TZZ  )
  CALL COMMUNICATE   ( TXZ  )
  CALL COMMUNICATE   ( PRES  )


END SUBROUTINE COMMUNICATE_STRESSFIELD
