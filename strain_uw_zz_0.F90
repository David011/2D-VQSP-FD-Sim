!=======================================================================

SUBROUTINE  STRAIN_UW_ZZ_0 ( MXL, MXH, MXT,                MZL, MZH,      &
                          X_MIN, X_MAX,               TPML, C,         &
                          WTM, EZZT                                    &
                          !,XZZT
                           )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                          MXL, MXH, MXT,                MZL, MZH, TPML,&
                          X_MIN, X_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                                     C      !1/H

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EZZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM

  !REAL   (PP), DIMENSION (N_FREQ, MXL:MXH,  MZL:MZH), INTENT(INOUT) :: &
                                                                  !XZZT

!-----------------------------------------------------------------------

  INTEGER       :: I, L, IO, LP1, MXT1

!-----------------------------------------------------------------------

L    = 0
LP1  = 1

MXT1  = MXT - 1  !MXT1 = MX - 1

!================================== EZZT ===============================
!----- ALL POINTS COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    EZZT(I) = ( C*( WTM (I  ,     LP1) - WTM (I  ,     L  ) ) ) 
                                                                              
  END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
#if(AFDA_ZZ)
  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    EZZT(I)   = C*( - 11._PP/12._PP * WTM (I,  0)                      & 
                    + 17._PP/24._PP * WTM (I,  1)                      &
                    +  3._PP/ 8._PP * WTM (I,  2)                      &
                    -  5._PP/24._PP * WTM (I,  3)                      &
                    +  1._PP/24._PP * WTM (I,  4) )

  END DO
#endif

END SUBROUTINE  STRAIN_UW_ZZ_0
