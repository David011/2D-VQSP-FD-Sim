!=======================================================================

SUBROUTINE  INH_V_QZ_INT_1 ( MXL, MXH, MXT,                MZL, MZH,   &
                            X_MIN, X_MAX,                              &
                            JMT,                                       &
                            REL1_QZT, REL2_QZT, REL3_QZT,              &
                            QZT,                                       &
                            TZZT, TXZT, PREST )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM
  USE AUXIL       , ONLY:                                              &
                          A, B

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXZT, TZZT, PREST
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    QZT

  INTEGER   :: I, L, IN1, IN2, IN3, IN4, JM1, IO

!-----------------------------------------------------------------------

L    = 1

IN2  = L-1
IN3  = L
IN4  = L+1
IN1  = L+2


  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    JM1 = JMT (I  ,L   )

#if(SEDOND)
    QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                   &
                         - REL1_QZT(JM1)  *                              &
                 (     TXZT (I+1,    L    ) - TXZT (I  ,    L    )       &
                 +     TZZT (I  ,     IN3 ) - TZZT (I  ,    IN2  )   )   &
                              - REL2_QZT(JM1) *                          &
                 (     PREST (I  ,     IN3 ) - PREST (I  ,    IN2  ) ) 

        
#else
    QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                   &
                       - REL1_QZT(JM1)  *                                &
                 ( A*( TXZT (I+2,    L    ) - TXZT (I-1,    L    )    )  &  
                 + B*( TXZT (I+1,    L    ) - TXZT (I  ,    L    )    )  &
                 - 31._PP/ 24._PP * TZZT  (I  ,    IN2  )                &   
                 + 29._PP/ 24._PP * TZZT  (I  ,    IN3  )                &   
                 -  3._PP/ 40._PP * TZZT  (I  ,    IN4  )                &   
                 +  1._PP/168._PP * TZZT  (I  ,    IN1  )              ) &
                              - REL2_QZT(JM1) *                          &
                (- 31._PP/ 24._PP * PREST (I  ,    IN2  )                &   
                 + 29._PP/ 24._PP * PREST (I  ,    IN3  )                &   
                 -  3._PP/ 40._PP * PREST (I  ,    IN4  )                &   
                 +  1._PP/168._PP * PREST (I  ,    IN1  )              )
    
#endif
  END DO


END SUBROUTINE  INH_V_QZ_INT_1