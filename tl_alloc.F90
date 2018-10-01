!=======================================================================

SUBROUTINE  TL_ALLOC

  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE AUXIL       , ONLY:                                              &
                          X_MIN_D, X_MAX_D,                            &
                          Z_MIN_D, Z_MAX_D                              
                                                
  USE STRESSFIELD , ONLY:                                              &
                  TXX ,       TZZ , TXZ, PRES
  USE WAVEFIELD   , ONLY:                                              &
                    UM,     WM,    QX,    QZ,    QXP,    QZP

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ALLOSTAT

!-----------------------------------------------------------------------

  ALLOCATE (  UM  (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              WM  (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              QX  (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              QZ  (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              QXP (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              QZP (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              !XXX ( N_FREQ,         X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              !XZZ ( N_FREQ,         X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              !XXZ ( N_FREQ,         X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              TXX (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              TZZ (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              PRES(                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
              TXZ (                 X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D), &
                                                       STAT=ALLOSTAT )
  
    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR T_ALLOC_02, STAT= ', ALLOSTAT
      STOP
    END IF

  UM      (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  WM      (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  QX      (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  QZ      (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  QXP     (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  QZP     (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  !XXX     (1:N_FREQ, X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.            
  !XZZ     (1:N_FREQ, X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  !XXZ     (1:N_FREQ, X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TXX     (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TZZ     (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  PRES    (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.
  TXZ     (          X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) = 0.


END  SUBROUTINE  TL_ALLOC
