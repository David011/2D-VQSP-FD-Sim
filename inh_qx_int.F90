SUBROUTINE INH_QX_INT( L, MXL, MXH, MXT, MZL, MZH,                                  &
                       X_MIN, X_MAX,                                                &
                       JMT,                                                         &
                       REL1_QXT, REL2_QXT, REL3_QXT,                                &
                       QXT,                                                         &
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
                             L, MXL, MXH, MXT,                MZL, MZH,&
                             X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                           REL1_QXT, REL2_QXT, REL3_QXT             
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                           TXXT,    TXZT,   PREST
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    QXT

  INTEGER              :: I, IN1, IN2, IN3, IN4, JM1, IO

!-----------------------------------------------------------------------
  
IN1  = L-1
IN2  = L+0
IN3  = L+1
IN4  = L+2

  DO  I = MAX(4,X_MIN), MIN(MXT-2,X_MAX) 

    JM1 = JMT (I  ,L   )

    QXT (I,  L ) = EXP(- REL3_QXT(JM1)) * QXT (I,  L )                 &
                       - REL1_QXT(JM1)  *                              &
                 ( A*( TXXT (I+1,    L    ) - TXXT (I-2,    L    )     & 
                     + TXZT (I  ,    IN4  ) - TXZT (I  ,    IN1  ) )   & 
                 + B*( TXXT (I  ,    L    ) - TXXT (I-1,    L    )     &
                     + TXZT (I  ,    IN3  ) - TXZT (I  ,    IN2  ) ) ) &
                                                - REL2_QXT(JM1) *      &
                ( A*( PREST (I+1,    L    ) - PREST (I-2,    L    ) )  &
                + B*( PREST (I  ,    L    ) - PREST (I-1,    L    ) ) )

  END DO


END SUBROUTINE  INH_QX_INT