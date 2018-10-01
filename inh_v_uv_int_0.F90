!=======================================================================

SUBROUTINE  INH_V_UV_INT_0( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT, QXPT,                                 &
                            REL1_UT, REL2_UT, REL3_UT, REL3_QXT,       &
                            UT, QXT,                                   &
                            TXXT, TXZT, PREST  )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM
  USE AUXIL       , ONLY:                                              &
                          A, B
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                            INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                          QXPT, UT, QXT
  INTEGER              :: I, L, JM1, IO

!-----------------------------------------------------------------------

L    = 0

  DO  I = MAX(4,X_MIN), MIN(MXT-2,X_MAX)

    JM1 = JMT (I  , L  )
                                                                       
#if(SECOND)
    UT (I,  L ) = UT (I,  L )                                                                       &
                  + (REL2_UT(JM1)/REL3_QXT(JM1)) * (1._PP - EXP(- REL3_QXT(JM1))) * QXPT  (I,  L )  &
                              + REL1_UT(JM1) *                                                      & 
                 (   ( TXXT  (I  ,      L  ) - TXXT  (I-1,      L  ) )                              &
                 +     TXZT  (I  ,      1  ) )                                                      &
                              + REL3_UT(JM1) *                                                      &
                 (   ( PREST (I  ,      L  ) - PREST (I-1,      L  ) )
    
#else
    UT (I,  L ) = UT (I,  L )                                                                       &
                  + (REL2_UT(JM1)/REL3_QXT(JM1)) * (1._PP - EXP(- REL3_QXT(JM1))) * QXPT  (I,  L )  &
                              + REL1_UT(JM1) *                                                      & 
                 ( A*( TXXT  (I+1,      L  ) - TXXT  (I-2,      L  )   )                            &
                 + B*( TXXT  (I  ,      L  ) - TXXT  (I-1,      L  )   )                            &
                 + 17._PP/24._PP * TXZT (I  ,      1  )                                             &
                 +  3._PP/ 8._PP * TXZT (I  ,      2  )                                             &
                 -  5._PP/24._PP * TXZT (I  ,      3  )                                             &
                 +  1._PP/24._PP * TXZT (I  ,      4  )                )                            &
                              + REL3_UT(JM1) *                                                      &
                 ( A*( PREST (I+1,      L  ) - PREST (I-2,      L  )   )                            &
                 + B*( PREST (I  ,      L  ) - PREST (I-1,      L  )  ) )
    
#endif
  END DO


END SUBROUTINE  INH_V_UV_INT_0
