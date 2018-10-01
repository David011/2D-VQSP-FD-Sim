SUBROUTINE  STRAIN_UW_ZZ( L, MXL, MXH, MXT,                MZL, MZH,   &
                          X_MIN, X_MAX,                                &
                          TPML, A, B, C,                               &
                          !XZZT,                                        &
                          WTM,  EZZT )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                           INTENT (IN)    :: &
                                 MXL, MXH, MXT,                TPML, L,&
                                 MZL, MZH,                             &
                                 X_MIN, X_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EZZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM

  !REAL   (PP), DIMENSION (N_FREQ, MXL:MXH,  MZL:MZH), INTENT(INOUT) :: &
                                                                  !XZZT

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
                                                                  
    EZZT(I  ) = ( C*( WTM (I  ,     LP1) - WTM (I  ,     L  ) ) )  
                                                                   
  END DO

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(2,X_MIN), MIN(MXT1,X_MAX)

    EZZT(I  ) = A  *( WTM (I  ,     LP2) - WTM (I  ,     LM1) )        &    
              + BC *  EZZT(I  )                                             
  END DO

END SUBROUTINE  STRAIN_UW_ZZ
