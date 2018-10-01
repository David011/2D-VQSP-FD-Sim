!=======================================================================

SUBROUTINE  INH_W_QZ_INT_1_SI ( MXL, MXH, MXT,                MZL, MZH, &
                            X_MIN, X_MAX,               TPML,          &
                            JMT, UT, WT, QXT, QZT )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          AUXIL_W, AUXIL_QZ

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX,               TPML
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                       UT, WT, QXT, QZT

  INTEGER   :: I,    L, JM1

!-----------------------------------------------------------------------

L    = 1

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-2+TPML,X_MAX)

    JM1 = JMT (I  ,0  )
                                                                               
    WT (I,  L) = WT  (I,  3)                                                  &    
               + AUXIL_W(JM1) * ( UT  (I+1,    1) - UT  (I  ,    1)           &
                              +   UT  (I+1,    2) - UT  (I  ,    2) )
    
    QZT(I,  L) = QZT  (I,  3)                                                 &
               + AUXIL_QZ(JM1)* ( UT  (I+1,    1) - UT  (I  ,    1)           &
                              +   UT  (I+1,    2) - UT  (I  ,    2) )         &
               +                ( QXT (I+1,    1) - QXT (I  ,    1)           &
                              +   QXT (I+1,    2) - QXT (I  ,    2) ) 

  END DO


END SUBROUTINE  INH_W_QZ_INT_1_SI
