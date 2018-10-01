!=======================================================================

SUBROUTINE  INH_V_W_INT_0 ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT, QZPT,                                 &
                            REL1_WT, REL2_WT, REL3_WT, REL3_QZT,       &
                            WT, QZT,                                   &
                            TZZT,  PREST )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                            INTENT(IN)    :: &
                                MXL, MXH, MXT,                MZL, MZH,&
                                X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                            TZZT, PREST
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                          QZPT, WT, QZT

  INTEGER :: I, L, JM1, IO

!-----------------------------------------------------------------------

L    = 0

  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)                  

    JM1 = JMT (I  ,L  )

    WT (I,  L ) = WT (I,  L )                                                                       &
                  + (REL2_WT(JM1)/REL3_QZT(JM1)) * (1._PP - EXP(- REL3_QZT(JM1))) * QZPT  (I,  L )  &
                              + REL1_WT(JM1) *                                                      &   
                 ( 35._PP/ 8._PP * TZZT  (I  ,      0  )                                            &  
                 - 35._PP/24._PP * TZZT  (I  ,      1  )                                            &
                 + 21._PP/40._PP * TZZT  (I  ,      2  )                                            &
                 -  5._PP/56._PP * TZZT  (I  ,      3  ) )                                          &
                              + REL3_WT(JM1) *                                                      &
                 ( 35._PP/ 8._PP * PREST (I  ,      0  )                                            & 
                 - 35._PP/24._PP * PREST (I  ,      1  )                                            &
                 + 21._PP/40._PP * PREST (I  ,      2  )                                            &
                 -  5._PP/56._PP * PREST (I  ,      3  ) )
    
  END DO

END SUBROUTINE  INH_V_W_INT_0
