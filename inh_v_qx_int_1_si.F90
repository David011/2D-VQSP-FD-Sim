!=======================================================================

SUBROUTINE  INH_V_QX_INT_1_SI                                          &
                          ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,               TPML,          &
                            UT, WT, QXT, QZT )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,                                            INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX,               TPML
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                     UT, WT, QXT, QZT

  INTEGER              :: I,    L

!-----------------------------------------------------------------------

L    = 1

  DO  I = MAX(3-TPML,X_MIN), MIN(MXT-1+TPML,X_MAX)

    UT (I,L) = UT (I,2) + WT (I,2) - WT (I-1,2)   
    QXT(I,L) = QXT(I,2) + QZT(I,2) - QZT(I-1,2)   
    
  END DO


END SUBROUTINE  INH_V_QX_INT_1_SI
