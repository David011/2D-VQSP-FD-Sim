SUBROUTINE  STRAIN_Q_ZZ( L, MXL, MXH, MXT,                MZL, MZH,    &
                          X_MIN, X_MAX,                                &
                          TPML  ,A, B, C,                              &
                          QZT,  QZZT )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                                 L, MXL, MXH, MXT, TPML,               &
                                 MZL, MZH,                             &
                                 X_MIN, X_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  QZZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT


!-----------------------------------------------------------------------

  INTEGER       :: I,    IO, LP1, LP2, LM1, MXT1

  REAL (PP)     :: BC

!-----------------------------------------------------------------------

LP1  = L + 1
LP2  = L + 2
LM1  = L - 1

MXT1  = MXT - 1

BC = B/C

!================================== EZZT ===============================
!----- ALL POINTS COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)
                                                                                 
    QZZT(I  ) = ( C*( QZT (I  ,     LP1) - QZT (I  ,     L  ) ) )                 
                                                                                 
  END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    QZZT(I  ) = A  *( QZT (I  ,     LP2) - QZT (I  ,     LM1) )        &         
              + BC *  QZZT(I  )                                                  

  END DO


END SUBROUTINE  STRAIN_Q_ZZ