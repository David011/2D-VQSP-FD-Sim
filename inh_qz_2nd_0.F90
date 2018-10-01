SUBROUTINE  INH_QZ_2ND_0( MXL, MXH,                                    &
                           MXT,                                        &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT,                                        &
                           REL1_QZT, REL2_QZT, REL3_QZT,               &
                           QZT,                                        &
                           TZZT, PREST)

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                    MXL, MXH,           MXT,      TPML,&
                                    MZL, MZH,                          &
                                    X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                            TZZT, PREST
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    QZT

  INTEGER              :: I, L, JM1, MXT1, IO

!-----------------------------------------------------------------------
  L    = 0

  MXT1 = MXT - 1

!_____________________________________________________________ LEFT ____
  DO  I = MAX(2 - TPML,X_MIN), MIN(2,X_MAX)
    JM1 = JMT (I  , L  )

     QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                 &
                        - REL1_QZT(JM1)  *                              &  
                 ( 3._PP         * TZZT  (I  ,      0  )                &   
                 - 1._PP/3._PP   * TZZT  (I  ,      1  ) )              &  
!                 ( 35._PP/ 8._PP * TZZT (I  ,      0  )                &
!                 - 35._PP/24._PP * TZZT (I  ,      1  )                &
!                 + 21._PP/40._PP * TZZT (I  ,      2  )                &
!                 -  5._PP/56._PP * TZZT (I  ,      3  ) )              &
                              - REL2_QZT(JM1) *                         &
                 ( 3._PP         * PREST (I  ,      0  )                &
                 - 1._PP/3._PP   * PREST (I  ,      1  ) )   
!                 ( 35._PP/ 8._PP * PREST (I  ,      0  )                &
!                 - 35._PP/24._PP * PREST (I  ,      1  )                &
!                 + 21._PP/40._PP * PREST (I  ,      2  )                &
!                 -  5._PP/56._PP * PREST (I  ,      3  ) )              &
    
  END DO

!____________________________________________________________ RIGHT ____
  DO  I = MAX(MXT1,X_MIN) , MIN(MXT1  + TPML,X_MAX)
    JM1 = JMT (I  , L  )

     QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                 &
                       - REL1_QZT(JM1)  *                               &  
                 ( 3._PP         * TZZT  (I  ,      0  )                &   
                 - 1._PP/3._PP   * TZZT  (I  ,      1  ) )              &  
!                 ( 35._PP/ 8._PP * TZZT (I  ,      0  )                &
!                 - 35._PP/24._PP * TZZT (I  ,      1  )                &
!                 + 21._PP/40._PP * TZZT (I  ,      2  )                &
!                 -  5._PP/56._PP * TZZT (I  ,      3  ) )              &
                              - REL2_QZT(JM1) *                         &
                 ( 3._PP         * PREST (I  ,      0  )                &
                 - 1._PP/3._PP   * PREST (I  ,      1  ) )   
!                 ( 35._PP/ 8._PP * PREST (I  ,      0  )                &
!                 - 35._PP/24._PP * PREST (I  ,      1  )                &
!                 + 21._PP/40._PP * PREST (I  ,      2  )                &
!                 -  5._PP/56._PP * PREST (I  ,      3  ) )              &
    
  END DO


END SUBROUTINE  INH_QZ_2ND_0
