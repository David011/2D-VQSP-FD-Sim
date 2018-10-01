SUBROUTINE  STO_QXP                                                    &
                ( L,       MXL, MXH,                                   &
                           MXT,  MX1,                                  &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           QXT, QXPT)
                
  USE PRECISION   , ONLY:                                              &
                          PP, IK
  !----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                        L, MXL, MXH, MXT, MX1,         &
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
  
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                             QXPT, QXT
  
  INTEGER   :: I
  
   DO  I = MAX(3 - TPML,X_MIN), MIN(MX1  + TPML,X_MAX)
       
       QXPT (I, L) = QXT (I, L)
       
   END DO
   
END SUBROUTINE STO_QXP 