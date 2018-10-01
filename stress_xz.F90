SUBROUTINE  STRESS_XZ (L    , MXL, MXH, MXT,                MZL, MZH,  &
                                          X_MIN, X_MAX,                &
                                          TXZT, EXZT,                  &
                                          !XXZT,                       &
                                          JMT, TPML )

  USE PRECISION   , ONLY: PP, IK
  USE CONTROL_DATA, ONLY: DT, TAU_EPS
  USE GRID_MEDIUM , ONLY: MU

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                       MXL, MXH, MXT,                MZL, MZH, L, TPML,&
                       X_MIN, X_MAX

  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT

  !REAL   (PP), DIMENSION (N_FREQ, MXL:MXH,  MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  TXZT

!-----------------------------------------------------------------------

  INTEGER       :: I, JM1, IO

  REAL (KIND=PP):: TXZ

!-----------------------------------------------------------------------

  DO  I = MAX(2-TPML,X_MIN), MIN(MXT+TPML,X_MAX)

    JM1     = JMT (I, L )

    TXZ =  MU(JM1) * EXZT(I)

    TXZT(I, L ) = TXZT(I, L ) + ( DT + TAU_EPS ) * TXZ
  END DO


  END SUBROUTINE  STRESS_XZ
