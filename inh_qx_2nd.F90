SUBROUTINE  INH_QX_2ND(L, MXL, MXH, PL,                                &
                           MXT, MX1, PXH,                              &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT,                                        &
                           REL1_QXT, REL2_QXT, REL3_QXT,               &
                           QXT,                                        &
                           TXXT, TXZT, PREST,                          &
                           GXXLT, GXXHT, HXXLT, HXXHT,                 &
     ODA2_XL,ODB2_XL,ODC2_XL,ODD2_XL, ODA2_XH,ODB2_XH,ODC2_XH,ODD2_XH  & 
                       )

  USE PRECISION   , ONLY:                                              &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                              &
                          JMNUM

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                           REL1_QXT, REL2_QXT, REL3_QXT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
  REAL   (PP), DIMENSION (PL :3  ,          MZL:MZH), INTENT(INOUT) :: &
                                                                  GXXLT
  REAL   (PP), DIMENSION (MX1:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  GXXHT
  REAL   (PP), DIMENSION (PL :3  ,          MZL:MZH), INTENT(INOUT) :: &
                                                                  HXXLT
  REAL   (PP), DIMENSION (MX1:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  HXXHT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    QXT
  REAL   (PP), DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

  INTEGER   :: I, IN2, IN3, JM1, IO
  REAL (PP) :: DXXXA, PXXXA

!----------------------------------------------------------------------

IN2  = L+0
IN3  = L+1

!_____________________________________________________________ LEFT ____
  DO  I = MAX(3-TPML,X_MIN), MIN(3,X_MAX)                                  

    JM1 = JMT (I  ,L   )
    
    DXXXA         = TXXT  (I,   L    ) - TXXT  (I-1,  L    )                  
    PXXXA         = PREST (I,   L    ) - PREST (I-1,  L    )

    QXT (I,  L ) =   EXP(- REL3_QXT(JM1)) * QXT (I,  L )                  & 
                         - REL1_QXT(JM1)  *                               &
                   ( ODA2_XL(3-I  ,  L) * ( ODD2_XL(3-I  ,  L)*DXXXA      &
                                          + GXXLT(I,  L)     )            &
                     + TXZT (I,     IN3  ) - TXZT (I  ,    IN2  )    )    &
                                                  - REL2_QXT(JM1) *       &
                  ( ODA2_XL(3-I  ,  L) * ( ODD2_XL(3-I  ,  L)*PXXXA       &
                                         + HXXLT(I,  L)     )      )
  
    GXXLT (I,  L )= ODB2_XL(3-I  ,  L) * GXXLT(I,  L)                     &
                  + ODC2_XL(3-I  ,  L) * DXXXA
    
    HXXLT (I,  L )= ODB2_XL(3-I  ,  L) * HXXLT(I,  L)                     &
                  + ODC2_XL(3-I  ,  L) * PXXXA
    
  END DO

!____________________________________________________________ RIGHT ____
  DO  I = MAX(MX1 ,X_MIN), MIN(MX1  + TPML,X_MAX)                           

    JM1 = JMT (I  ,L   )

    DXXXA         = TXXT  (I,   L    ) - TXXT  (I-1,  L    )
    PXXXA         = PREST (I,   L    ) - PREST (I-1,  L    )
    
    QXT (I,  L ) =   EXP(- REL3_QXT(JM1)) * QXT (I,  L )                  &
                         - REL1_QXT(JM1)  *                               &
                   ( ODA2_XH(I-MX1,  L) * ( ODD2_XH(I-MX1,  L)*DXXXA      &
                                          + GXXHT(I,  L)     )            &
                   + TXZT (I,     IN3  ) - TXZT (I  ,    IN2  )    )      &
                                                  - REL2_QXT(JM1) *       &
                   ( ODA2_XH(I-MX1,  L) * ( ODD2_XH(I-MX1,  L)*PXXXA      &
                                          + HXXHT(I,  L)     )     )
  
    GXXHT (I,  L )= ODB2_XH(I-MX1  ,L) * GXXHT(I,  L)                     &
                  + ODC2_XH(I-MX1  ,L) * DXXXA
    
    HXXHT (I,  L )= ODB2_XH(I-MX1  ,L) * HXXHT(I,  L)                     &
                  + ODC2_XH(I-MX1  ,L) * PXXXA
     
   

  END DO


END SUBROUTINE  INH_QX_2ND