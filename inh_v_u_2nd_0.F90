
SUBROUTINE  INH_V_U_2ND_0 (MXL, MXH, PL,                               &
                           MXT, MX1, PXH,                              &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT, QXPT,                                  &
                           REL1_UT, REL2_UT, REL3_UT, REL3_QXT,        &
                           UT, QXT,                                    &
                           TXXT, TXZT, PREST,                          &
                           RXXLT, RXXHT, SXXLT, SXXHT,                 &
     ODA2_XL,ODB2_XL,ODC2_XL,ODD2_XL, ODA2_XH,ODB2_XH,ODC2_XH,ODD2_XH  & 
                          )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                           MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST 
  REAL   (PP), DIMENSION (PL :3  ,          MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXLT
  REAL   (PP), DIMENSION (MX1:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  RXXHT
  REAL   (PP), DIMENSION (PL :3  ,          MZL:MZH), INTENT(INOUT) :: &
                                                                  SXXLT
  REAL   (PP), DIMENSION (MX1:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  SXXHT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                          QXPT, UT, QXT
  REAL   (PP), DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

  INTEGER   :: L, I,      IN3, JM1, IO
  REAL (PP) :: DXXXA, PXXXA

!----------------------------------------------------------------------
L    = 0

IN3  = L+1

!_____________________________________________________________ LEFT ____
  DO  I = MAX(3-TPML,X_MIN), MIN(3,X_MAX)

    JM1 = JMT (I  , L  )
                                                                          
    DXXXA         = TXXT (I,     L    ) - TXXT (I-1,    L    )         
    PXXXA         = PREST (I,   L    ) - PREST (I-1,  L    )
    
    UT (I,  L ) = UT (I,  L )                                                                       &
                  + (REL2_UT(JM1)/REL3_QXT(JM1)) * (1._PP - EXP(- REL3_QXT(JM1))) * QXPT  (I,  L )  &
                              + REL1_UT(JM1) *                                                      &
                  ( ODA2_XL(3-I  ,  L) * ( ODD2_XL(3-I  ,  L)*DXXXA                                 &
                                         + RXXLT(I,  L)     )                                       &
                  + TXZT (I,     IN3  )                        )                                    &
                                               + REL3_UT(JM1)*                                      &
                  ( ODA2_XL(3-I  ,  L) * ( ODD2_XL(3-I  ,  L)*PXXXA                                 &
                                         + SXXLT(I,  L)     )      )

    RXXLT (I,  L )= ODB2_XL(3-I  ,  L) * RXXLT(I,  L)                                               &
                  + ODC2_XL(3-I  ,  L) * DXXXA
    
    SXXLT (I,  L )= ODB2_XL(3-I  ,  L) * SXXLT(I,  L)                                               &
                  + ODC2_XL(3-I  ,  L) * PXXXA

  END DO

!____________________________________________________________ RIGHT ____
  DO  I = MAX(MX1 ,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  , L  )

    DXXXA         = TXXT (I,     L    ) - TXXT (I-1,    L    )
    PXXXA         = PREST (I,   L    ) - PREST (I-1,  L    )
    
    UT (I,  L ) = UT (I,  L )                                                                       &
                  + (REL2_UT(JM1)/REL3_QXT(JM1)) * (1._PP - EXP(- REL3_QXT(JM1))) * QXPT  (I,  L )  &
                              + REL1_UT(JM1) *                                                      &
                  ( ODA2_XH(I-MX1,  L) * ( ODD2_XH(I-MX1,  L)* DXXXA                                &
                                         + RXXHT(I,  L)              )                              &
                  + TXZT (I,     IN3  )                        )                                    &
                              + REL3_UT(JM1)*                                                       &
                  ( ODA2_XH(I-MX1,  L) * ( ODD2_XH(I-MX1,  L)* PXXXA                                &
                                         + SXXHT(I,  L)     )          )

    RXXHT (I,  L )= ODB2_XH(I-MX1  ,L) * RXXHT(I,  L)                                               &
                  + ODC2_XH(I-MX1  ,L) * DXXXA

    SXXHT (I,  L )= ODB2_XH(I-MX1  ,L) * SXXHT(I,  L)                                               &
                  + ODC2_XH(I-MX1  ,L) * PXXXA

  END DO


END SUBROUTINE  INH_V_U_2ND_0
