SUBROUTINE COMMUNICATE_WAVEFIELD

  USE WAVEFIELD   , ONLY: UM, WM, QX, QZ

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

  END INTERFACE

!-----------------------------------------------------------------------

  CALL COMMUNICATE   ( QX  )
  CALL COMMUNICATE   ( UM  )
  CALL COMMUNICATE   ( QZ  )
  CALL COMMUNICATE   ( WM  )
  
END SUBROUTINE COMMUNICATE_WAVEFIELD
