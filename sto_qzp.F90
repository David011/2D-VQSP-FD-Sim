SUBROUTINE  STO_QZP                                                    &
                ( L,       MXL, MXH,                                   &
                           MXT,  MX1,                                  &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           QZT, QZPT)
                
  USE PRECISION   , ONLY:                                              &
                          PP, IK
  !----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                        L, MXL, MXH, MXT, MX1,         &
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
  
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                             QZPT, QZT
  
  INTEGER   :: I
  
   DO  I = MAX(2 - TPML,X_MIN), MIN(MX1  + TPML,X_MAX)
       
       QZPT (I, L) = QZT (I, L)
       
   END DO
   
END SUBROUTINE STO_QZP