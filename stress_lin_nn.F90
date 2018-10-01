SUBROUTINE  STRESS_LIN_NN (L    , MXL, MXH, MXT,                MZL, MZH,  &
                     X_MIN, X_MAX,                                     &
                     TXXT,                   TZZT,        PREST,       &
                     EXXT,                   EZZT,                     &
                    !XXXT,              XZZT,                          &
                     QXXT,                   QZZT,                     &
                     JMT,                    TPML,  SOURS,  SOURF )

  USE PRECISION   , ONLY: PP, IK
  USE CONTROL_DATA, ONLY: DT, TAU_EPS, LIN_X, SOURSP, SOURFP 
  USE GRID_MEDIUM , ONLY:  SIG_XX_1,       SIG_XX_2,       SIG_XX_3,   &
                           SIG_ZZ_1,       SIG_ZZ_2,       SIG_ZZ_3,   &
                           PRES_1,         PRES_2,         PRES_3

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                       MXL, MXH, MXT,                MZL, MZH, L, TPML,&
                       X_MIN, X_MAX

  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                      EXXT,       EZZT
  
  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                      QXXT,       QZZT

  !REAL   (PP), DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                      !XXXT,       XZZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT,       TZZT,   PREST
  
  REAL   (PP),                                        INTENT(IN) :: SOURS, SOURF

!-----------------------------------------------------------------------

  INTEGER       :: I, IO, JM1

  REAL (KIND=PP):: TXX, TZZ, PRES

!-----------------------------------------------------------------------

    DO  I = MAX(2 - TPML,X_MIN), MIN(MXT-1 + TPML,X_MAX)

      JM1     = JMT (I, L)                                  

     TXX  =    SIG_XX_1(JM1)     *   EXXT(I)   +  SIG_XX_2(JM1)    *  EZZT(I)   +  SIG_XX_3(JM1)  *  QXXT(I)  +  SIG_XX_3(JM1)  *  QZZT(I)
     TZZ  =    SIG_ZZ_1(JM1)     *   EXXT(I)   +  SIG_ZZ_2(JM1)    *  EZZT(I)   +  SIG_ZZ_3(JM1)  *  QXXT(I)  +  SIG_ZZ_3(JM1)  *  QZZT(I)
     PRES =    PRES_1  (JM1)     *   EXXT(I)   +  PRES_2  (JM1)    *  EZZT(I)   +  PRES_3(JM1)    *  QXXT(I)  +  PRES_3(JM1)    *  QZZT(I)
      
      
      IF (I == LIN_X) THEN
          TXXT (I,L) = TXXT (I,L) + ( DT + TAU_EPS) * (TXX  - ((SOURS - SOURSP)/DT) )                      
          TZZT (I,L) = TZZT (I,L) + ( DT + TAU_EPS) * (TZZ  - ((SOURS - SOURSP)/DT) )
          PREST(I,L) = PREST(I,L) + ( DT + TAU_EPS) * (PRES + ((SOURF - SOURFP)/DT) )
          SOURSP = SOURS
          SOURFP = SOURF
      ELSE
          TXXT (I,L) = TXXT (I,L) + ( DT + TAU_EPS) *  TXX                    
          TZZT (I,L) = TZZT (I,L) + ( DT + TAU_EPS) *  TZZ
          PREST(I,L) = PREST(I,L) + ( DT + TAU_EPS) *  PRES
      END IF
      
    END DO


END SUBROUTINE  STRESS_LIN_NN