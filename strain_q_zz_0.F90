SUBROUTINE  STRAIN_Q_ZZ_0 ( MXL, MXH, MXT,                MZL, MZH,      &
                          X_MIN, X_MAX,               TPML, C,           &
                          QZT, QZZT)                                    
                           

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                          MXL, MXH, MXT,                MZL, MZH, TPML,&
                          X_MIN, X_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                                     C      !1/H

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  QZZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT

 
!-----------------------------------------------------------------------

  INTEGER       :: I, L, IO, LP1, MXT1

!-----------------------------------------------------------------------

L    = 0
LP1  = 1

MXT1  = MXT - 1

!================================== EZZT ===============================
!----- ALL POINTS COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    QZZT(I) = ( C*( QZT (I  ,     LP1) - QZT (I  ,     L  ) ) )
                                                                              
  END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
#if(AFDA_ZZ)
  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    QZZT(I)   = C*( - 11._PP/12._PP * QZT (I,  0)                      & 
                    + 17._PP/24._PP * QZT (I,  1)                      & 
                    +  3._PP/ 8._PP * QZT (I,  2)                      &
                    -  5._PP/24._PP * QZT (I,  3)                      &
                    +  1._PP/24._PP * QZT (I,  4) )

  END DO
#endif


END SUBROUTINE  STRAIN_Q_ZZ_0