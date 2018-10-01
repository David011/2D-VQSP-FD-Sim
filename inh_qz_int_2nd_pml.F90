SUBROUTINE  INH_QZ_INT_2ND_PML                                          &
                   ( L, LL, MXL, MXH, PL, MXT, MX1, PXH,                &
                           MZL, MZH,  TPML,                             &
                           X_MIN, X_MAX,                                &
                           JMT,                                         &
                           REL1_QZT, REL2_QZT, REL3_QZT,                &
                           QZT,                                         &
                           TZZT, TXZT, PREST,                           &
                           GZZHT, HZZHT,                                &
                           GZXLT, GZXHT,                                &
ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,  &
                                    ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH  &
                   )
                       
  USE PRECISION   , ONLY:                                               &
                          PP, IK
  USE GRID_MEDIUM , ONLY:                                               &
                          JMNUM

!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER    ,                                        INTENT(IN)    :: &
                                    L, LL, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
  INTEGER(IK), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
  REAL   (PP), DIMENSION (1:JMNUM                  ), INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TZZT, TXZT, PREST
  REAL   (PP), DIMENSION (PL :2  ,          MZL:MZH), INTENT(INOUT) :: &
                                                                  GZXLT
  REAL   (PP), DIMENSION (MX1:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  GZXHT
  REAL   (PP), DIMENSION (MXL:MXH,           0:TPML), INTENT(INOUT) :: &
                                                                  GZZHT
  REAL   (PP), DIMENSION (MXL:MXH,           0:TPML), INTENT(INOUT) :: &
                                                                  HZZHT
  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                    QZT
  REAL   (PP), DIMENSION ( 0:TPML,          MZL:MZH), INTENT(INOUT) :: &
  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
  REAL   (PP), DIMENSION(0:TPML+1,          MXL:MXH), INTENT(INOUT) :: &
                                     ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
  INTEGER   :: I, IN2, IN3, JM1, IO
  REAL (PP) :: DXZXA, DZZZA, PZZZA

!----------------------------------------------------------------------

IN2  = L-1
IN3  = L

!_____________________________________________________________ LEFT ____
  DO  I = MAX(2 - TPML,X_MIN), MIN(2,X_MAX)

    JM1 = JMT (I  ,L   )

    DXZXA         = TXZT  (I+1,    L    ) - TXZT  (I  ,    L    )
    DZZZA         = TZZT  (I  ,    IN3  ) - TZZT  (I  ,    IN2  )
    PZZZA         = PREST (I  ,    IN3  ) - PREST (I  ,    IN2  )  
    
    QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                    &
                       - REL1_QZT(JM1)  *                                 &
                  ( ODA2_ZH(LL   ,I  ) * ( ODD2_ZH(LL   ,I  )*DZZZA       &
                                        + GZZHT(I,  LL)            ) )    &
                              - REL1_QZT(JM1) *                           &
                  ( ODA_XL (2-I  ,  L) * ( ODD_XL (2-I  ,  L)*DXZXA       &
                                        + GZXLT(I,  L)             ) )    &
                              - REL2_QZT(JM1) *                           &
                  ( ODA2_ZH(LL   ,I  ) * ( ODD2_ZH(LL   ,I  )*PZZZA       &
                                        + HZZHT(I,  LL)            ) )


    GZXLT (I,  L )= ODB_XL (2-I,  L) * GZXLT(I,  L)                       &
                  + ODC_XL (2-I,  L) * DXZXA

    GZZHT (I,  LL)= ODB2_ZH(LL ,I  ) * GZZHT(I,  LL)                      &
                  + ODC2_ZH(LL ,I  ) * DZZZA
    
    HZZHT (I,  LL)= ODB2_ZH(LL ,I  ) * HZZHT(I,  LL)                      &
                  + ODC2_ZH(LL ,I  ) * PZZZA

  END DO
  
!_________________________________________________________________ INTERIOR____
  
  DO  I = MAX(3,X_MIN), MIN(MX1 -1,X_MAX)
    JM1 = JMT (I  ,L   )

    DXZXA         = TXZT  (I+1,    L    ) - TXZT  (I  ,    L    )
    DZZZA         = TZZT  (I  ,    IN3  ) - TZZT  (I  ,    IN2  )
    PZZZA         = PREST (I  ,    IN3  ) - PREST (I  ,    IN2  )
    
    
    QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                    &
                       - REL1_QZT(JM1)  *                                 &
                  ( ODA2_ZH(LL   ,I  ) * ( ODD2_ZH(LL   ,I  )*DZZZA       &
                                        + GZZHT(I,  LL)             )     &
                                                            + DXZXA )     &         
                              - REL2_QZT(JM1) *                           &
                  ( ODA2_ZH(LL   ,I  ) * ( ODD2_ZH(LL   ,I  )*PZZZA       &
                                        + HZZHT(I,  LL)            ) )     

    GZZHT (I,  LL)= ODB2_ZH(LL ,I  ) * GZZHT(I,  LL)                      &
                  + ODC2_ZH(LL ,I  ) * DZZZA
    
    HZZHT (I,  LL)= ODB2_ZH(LL ,I  ) * HZZHT(I,  LL)                      &
                  + ODC2_ZH(LL ,I  ) * PZZZA
  END DO
  
!____________________________________________________________ RIGHT ____

   DO  I = MAX(MX1 ,X_MIN), MIN(MX1  + TPML,X_MAX)

    JM1 = JMT (I  ,L   )

    DXZXA         = TXZT  (I+1,    L    ) - TXZT  (I  ,    L    )
    DZZZA         = TZZT  (I  ,    IN3  ) - TZZT  (I  ,    IN2  )
    PZZZA         = PREST (I  ,    IN3  ) - PREST (I  ,    IN2  )
    
    QZT (I,  L ) = EXP(- REL3_QZT(JM1)) * QZT (I,  L )                    &
                       - REL1_QZT(JM1)  *                                 &
                  ( ODA2_ZH(LL   ,I  ) * ( ODD2_ZH(LL   ,I  )*DZZZA       &
                                        + GZZHT(I,  LL)            ) )    &
                              - REL1_QZT(JM1) *                           &
                  ( ODA_XH (I-MX1  ,  L) * ( ODD_XH (I-MX1  ,  L)*DXZXA   &
                                        + GZXHT(I,  L)             ) )    &
                              - REL2_QZT(JM1) *                           &
                  ( ODA2_ZH(LL   ,I  ) * ( ODD2_ZH(LL   ,I  )*PZZZA       &
                                        + HZZHT(I,  LL)            ) )
    
    GZXHT (I,  L )= ODB_XH (I-MX1,  L) * GZXHT(I,  L)                     &
                  + ODC_XH (I-MX1,  L) * DXZXA

    GZZHT (I,  LL)= ODB2_ZH(LL ,I  ) * GZZHT(I,  LL)                      &
                  + ODC2_ZH(LL ,I  ) * DZZZA
    
    HZZHT (I,  LL)= ODB2_ZH(LL ,I  ) * HZZHT(I,  LL)                      &
                  + ODC2_ZH(LL ,I  ) * PZZZA
    
  END DO
  
  END SUBROUTINE INH_QZ_INT_2ND_PML