MODULE  INTERFACES

  INTERFACE


    SUBROUTINE      ALLOC_NR
    END SUBROUTINE  ALLOC_NR

    SUBROUTINE      ALLOC_T_SEIS
    END SUBROUTINE  ALLOC_T_SEIS

    SUBROUTINE      AUX1
    END SUBROUTINE  AUX1
    
    SUBROUTINE      COMMUNICATE_STRESSFIELD
    END SUBROUTINE  COMMUNICATE_STRESSFIELD

    SUBROUTINE      COMMUNICATE_WAVEFIELD
    END SUBROUTINE  COMMUNICATE_WAVEFIELD

    SUBROUTINE      DISTRIBUTE
    END SUBROUTINE  DISTRIBUTE
    
    SUBROUTINE      DOMAIN_DECOMPOSITION
    END SUBROUTINE  DOMAIN_DECOMPOSITION

    SUBROUTINE      HPLANE_0_STRESS
    END SUBROUTINE  HPLANE_0_STRESS

    SUBROUTINE      HPLANE_0_VELOCITY
    END SUBROUTINE  HPLANE_0_VELOCITY

    SUBROUTINE      HPLANE_1_STRESS
    END SUBROUTINE  HPLANE_1_STRESS

    SUBROUTINE      HPLANE_1_VELOCITY
    END SUBROUTINE  HPLANE_1_VELOCITY

    SUBROUTINE      HPLANE_IC_STRESS ( L )
      INTEGER, INTENT(IN) ::           L
    END SUBROUTINE  HPLANE_IC_STRESS

    SUBROUTINE      HPLANE_IC_VELOCITY ( L )
      INTEGER, INTENT(IN) ::             L
    END SUBROUTINE  HPLANE_IC_VELOCITY

    SUBROUTINE      HPLANE_LIN_STRESS (SOURTF_LIN, SOURTF_LOU)
      USE PRECISION, ONLY: PP
      REAL(PP), INTENT(IN) ::          SOURTF_LIN, SOURTF_LOU
    END SUBROUTINE  HPLANE_LIN_STRESS

    SUBROUTINE      HPLANE_LIN_VELOCITY
    END SUBROUTINE  HPLANE_LIN_VELOCITY

    SUBROUTINE      HPLANE_LINPM1_STRESS ( L )
      INTEGER, INTENT(IN) ::               L
    END SUBROUTINE  HPLANE_LINPM1_STRESS

    SUBROUTINE      HPLANE_LINPM1_VELOCITY ( L )
      INTEGER, INTENT(IN) ::                 L
    END SUBROUTINE  HPLANE_LINPM1_VELOCITY

    SUBROUTINE      HPLANE_MZ1_STRESS
    END SUBROUTINE  HPLANE_MZ1_STRESS

    SUBROUTINE      HPLANE_PML_STRESS ( L, LL )
      INTEGER, INTENT(IN) ::            L, LL
    END SUBROUTINE  HPLANE_PML_STRESS

    SUBROUTINE      HPLANE_PML_VELOCITY ( L, LL )
      INTEGER, INTENT(IN) ::              L, LL
    END SUBROUTINE  HPLANE_PML_VELOCITY

    SUBROUTINE      HPLANE_PMLZH0_STRESS
    END SUBROUTINE  HPLANE_PMLZH0_STRESS

    SUBROUTINE      HPLANE_PMLZH0_VELOCITY
    END SUBROUTINE  HPLANE_PMLZH0_VELOCITY

!-------------------------------------------------------------------------------INH_2ND
    SUBROUTINE      INH_V_U_2ND                                        &
                       (L, MXL, MXH, PL,                               &
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
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                             JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
      REAL   (PP),DIMENSION(PL :3  ,        MZL:MZH), INTENT(INOUT) :: &
                                                           RXXLT, SXXLT
      REAL   (PP),DIMENSION(MX1:PXH,        MZL:MZH), INTENT(INOUT) :: &
                                                           RXXHT, SXXHT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                          UT, QXT, QXPT
      REAL   (PP),DIMENSION(0:TPML+1,        MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_V_U_2ND
                 
    SUBROUTINE      INH_QX_2ND                                         &
                       (L, MXL, MXH, PL,                               &
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
    USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                          REL1_QXT, REL2_QXT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                     TXXT, TXZT, PREST
      REAL   (PP), DIMENSION (PL :3  ,      MZL:MZH), INTENT(INOUT) :: &
                                                          GXXLT, HXXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                          GXXHT, HXXHT
      REAL   (PP), DIMENSION (MXL:MXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                   QXT
      REAL   (PP), DIMENSION (0:TPML+1,     MZL:MZH), INTENT(INOUT) :: &
  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_QX_2ND
 
    SUBROUTINE      INH_V_W_2ND                                        &
                       (L, MXL, MXH, PL,                               &
                           MXT, MX1, PXH,                              &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT, QZPT,                                  &
                           REL1_WT, REL2_WT, REL3_WT, REL3_QZT,        &
                           WT, QZT,                                    &
                           TZZT, TXZT, PREST,                          &
                           RZXLT, RZXHT,                               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH                    &
                       )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                             L, MXL, MXH, PL, MXT, PXH,&
                                                MZL, MZH,  TPML,       &
                                                MX1,                   &
                                             X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TXZT, TZZT, PREST
      REAL   (PP),DIMENSION(PL :2  ,         MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXLT
      REAL   (PP),DIMENSION(MX1:PXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QZPT, WT, QZT
      REAL   (PP),DIMENSION( 0:TPML,         MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
    END SUBROUTINE  INH_V_W_2ND
    
    SUBROUTINE      INH_QZ_2ND                                         &
                       (L, MXL, MXH, PL,                               &
                           MXT, MX1, PXH,                              &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT,                                        &
                           REL1_QZT, REL2_QZT, REL3_QZT,               &
                           QZT,                                        &
                           TZZT, TXZT, PREST,                          &
                           GZXLT, GZXHT,                               &
ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH  &
                       )
    USE PRECISION   , ONLY: PP, IK
    USE GRID_MEDIUM , ONLY: JMNUM
      
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                        L, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                                   JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                     TZZT, TXZT, PREST
      REAL   (PP), DIMENSION (PL :2  ,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GZXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GZXHT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT
      REAL   (PP), DIMENSION ( 0:TPML,      MZL:MZH), INTENT(INOUT) :: &
  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
    END SUBROUTINE  INH_QZ_2ND
                      
 
!-------------------------------------------------------------------------------INH_INT_2ND_PML                         
                    
    SUBROUTINE      INH_V_U_INT_2ND_PML                                &
                  ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT, QXPT,                                  &
                           REL1_UT, REL2_UT, REL3_UT, REL3_QXT,        &
                           UT, QXT,                                    &
                           TXXT, TXZT, PREST,                          &
                           RXXLT, RXXHT, SXXLT, SXXHT,                 &
                           RXZHT,                                      &
                 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                   &
                 ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,                   &
                 ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH                     & 
                   )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    L, LL, MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
      REAL   (PP),DIMENSION(PL :3  ,         MZL:MZH),INTENT(INOUT) :: &
                                                           RXXLT, SXXLT 
      REAL   (PP),DIMENSION(MX1:PXH,         MZL:MZH),INTENT(INOUT) :: &
                                                           RXXHT, SXXHT
      REAL   (PP),DIMENSION(MXL:MXH,          0:TPML),INTENT(INOUT) :: &
                                                                  RXZHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QXPT, UT, QXT
      REAL   (PP),DIMENSION(0:TPML+1,        MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
      REAL   (PP),DIMENSION( 0:TPML,         MXL:MXH),INTENT(INOUT) :: &
                                    ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH
    END SUBROUTINE  INH_V_U_INT_2ND_PML
                  
    SUBROUTINE      INH_V_W_INT_2ND_PML                                &
                  ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT, QZPT,                                  &
                           REL1_WT, REL2_WT, REL3_WT, REL3_QZT,        &
                           WT, QZT,                                    &
                           TZZT, TXZT, PREST,                          &
                           RZZHT, SZZHT,                               &
                           RZXLT, RZXHT,                               &
                  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                   &
                  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,                   &
                 ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH                    &
                   )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                         L, LL, MXL, MXH, PL, MXT, PXH,&
                                                MZL, MZH,  TPML,       &
                                                MX1,                   &
                                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TZZT, TXZT, PREST
      REAL   (PP),DIMENSION(PL :2  ,         MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXLT
      REAL   (PP),DIMENSION(MX1:PXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                  RZXHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QZPT, WT, QZT
      REAL   (PP),DIMENSION(MXL:MXH,          0:TPML),INTENT(INOUT) :: &
                                                           RZZHT, SZZHT
      REAL   (PP),DIMENSION( 0:TPML,         MZL:MZH),INTENT(INOUT) :: &
                    ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,                 &
                    ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
      REAL   (PP),DIMENSION(0:TPML+1,MXL:MXH         ),INTENT(INOUT):: &
                   ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
      END SUBROUTINE  INH_V_W_INT_2ND_PML
                           
      SUBROUTINE INH_QX_INT_2ND_PML                                    &
                       ( L, LL, MXL, MXH, PL, MXT, PXH, MX1,           &
                                MZL, MZH,  TPML,                       &
                                X_MIN, X_MAX,                          &
                                JMT,                                   &
                                REL1_QXT, REL2_QXT, REL3_QXT,          &
                                QXT,                                   &
                                TXXT, TXZT, PREST,                     &
                                GXXLT, GXXHT, HXXLT, HXXHT,            &
                                GXZHT,                                 &
                                ODA2_XL, ODB2_XL,                      &
                                ODC2_XL, ODD2_XL,                      &
                                ODA2_XH, ODB2_XH,                      &
                                ODC2_XH, ODD2_XH,                      &
                                ODA_ZH,  ODB_ZH,                       &
                                ODC_ZH,  ODD_ZH                        &
                       )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    L, LL, MXL, MXH, PL, MXT ,MX1, PXH,&
                                                MZL, MZH,  TPML,       &
                                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                           REL1_QXT, REL2_QXT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
      REAL   (PP), DIMENSION (PL :3  ,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GXXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GXXHT
      REAL   (PP), DIMENSION (PL :3  ,      MZL:MZH), INTENT(INOUT) :: &
                                                                  HXXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                  HXXHT
      REAL   (PP), DIMENSION (MXL:MXH,      0:TPML), INTENT(INOUT) :: &
                                                                  GXZHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    QXT
      REAL   (PP), DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
      REAL   (PP), DIMENSION ( 0:TPML,          MXL:MXH), INTENT(INOUT) :: &
                                      ODA_ZH,  ODB_ZH,  ODC_ZH,  ODD_ZH
      END SUBROUTINE INH_QX_INT_2ND_PML
                       
      SUBROUTINE INH_QZ_INT_2ND_PML                                    &
        ( L, LL, MXL, MXH, PL, MXT, MX1, PXH,                          &
                                MZL, MZH,  TPML,                       &
                                X_MIN, X_MAX,                          &
                                JMT,                                   &
                                REL1_QZT, REL2_QZT, REL3_QZT,          &
                                QZT,                                   &
                                TZZT, TXZT, PREST,                     &
                                GZZHT, HZZHT,                          &
                                GZXLT, GZXHT,                          &
                                ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,     &
                                ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH,     &
                                ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH     &
                       )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                     L, LL, MXL, MXH, PL, MXT, MX1,PXH,&
                                                MZL, MZH,  TPML,       &
                                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP), DIMENSION (1:JMNUM              ), INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TZZT, TXZT, PREST
      REAL   (PP), DIMENSION (PL :2  ,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GZXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GZXHT
      REAL   (PP), DIMENSION (MXL:MXH,       0:TPML), INTENT(INOUT) :: &
                                                                  GZZHT
      REAL   (PP), DIMENSION (MXL:MXH,       0:TPML), INTENT(INOUT) :: &
                                                                  HZZHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    QZT
      REAL   (PP), DIMENSION ( 0:TPML,      MZL:MZH), INTENT(INOUT) :: &
  ODA_XL,  ODB_XL,  ODC_XL,  ODD_XL,  ODA_XH,  ODB_XH,  ODC_XH,  ODD_XH
  REAL   (PP), DIMENSION(0:TPML+1,          MXL:MXH), INTENT(INOUT) :: &
                                     ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
      END SUBROUTINE INH_QZ_INT_2ND_PML
!-------------------------------------------------------------------------------INH_INT                  
                 
    SUBROUTINE      INH_V_UV_INT                                       &
                        (L, MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT, QXPT,                                 &
                            REL1_UT, REL2_UT, REL3_UT, REL3_QXT,       &
                            UT,  QXT,                                  &
                            TXXT,   TXZT,    PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                             L, MXL, MXH, MXT,                MZL, MZH,&
                             X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                           TXXT,    TXZT,   PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                    UT, QXPT, QXT
    END SUBROUTINE  INH_V_UV_INT

    SUBROUTINE      INH_QX_INT                                       &
                        (L, MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT,                                       &
                            REL1_QXT, REL2_QXT, REL3_QXT,              &
                            QXT,                                       &
                            TXXT,       TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      USE CONTROL_DATA, ONLY: N_FREQ
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                             L, MXL, MXH, MXT,                MZL, MZH,&
                             X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                           REL1_QXT, REL2_QXT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                           TXXT,    TXZT,    PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                               QXT
    END SUBROUTINE  INH_QX_INT    
                            
    SUBROUTINE      INH_V_W_INT                                        &
                       ( L, MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT,  QZPT,                                &
                            REL1_WT, REL2_WT, REL3_WT, REL3_QZT,       &
                            WT, QZT,                                   &
                            TZZT,   TXZT,    PREST )

      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                            L, MXL, MXH, MXT,                MZL, MZH, &
                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                    TXZT,   TZZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          WT, QZT, QZPT
    END SUBROUTINE  INH_V_W_INT
                            
    SUBROUTINE      INH_QZ_INT                                       &
                        (L, MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT,                                       &
                            REL1_QZT, REL2_QZT, REL3_QZT,              &
                            QZT,                                       &
                            TZZT,       TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                             L, MXL, MXH, MXT,                MZL, MZH,&
                             X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                           TZZT,    TXZT,    PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    QZT
    END SUBROUTINE  INH_QZ_INT                        
     
!-------------------------------------------------------------------------------INH_INT_0                       
                            
    SUBROUTINE      INH_V_UV_INT_0                                     &
                          ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT, QXPT,                                 &
                            REL1_UT, REL2_UT, REL3_UT, REL3_QXT,       &
                            UT, QXT,                                   &
                            TXXT , TXZT, PREST  )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                     TXXT , TXZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QXPT, UT, QXT
    END SUBROUTINE  INH_V_UV_INT_0
                            
    SUBROUTINE      INH_QX_INT_0                                     &
                          ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT,                                       &
                            REL1_QXT, REL2_QXT, REL3_QXT,              &
                            QXT,                                       &
                            TXXT, TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                           REL1_QXT, REL2_QXT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    QXT
   END SUBROUTINE  INH_QX_INT_0
                            
   SUBROUTINE      INH_V_W_INT_0                                       &
                       ( MXL, MXH, MXT,                MZL, MZH,       &
                            X_MIN, X_MAX,                              &
                            JMT, QZPT,                                 &
                            REL1_WT, REL2_WT, REL3_WT, REL3_QZT,       &
                            WT, QZT,                                   &
                            TZZT , TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                     TZZT , TXZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QZPT, WT, QZT
   END SUBROUTINE  INH_V_W_INT_0
                            
   SUBROUTINE      INH_QZ_INT_0                                     &
                          ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT,                                       &
                            REL1_QZT, REL2_QZT, REL3_QZT,              &
                            QZT,                                       &
                            TZZT, TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TZZT, TXZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    QZT
    END SUBROUTINE  INH_QZ_INT_0

!-------------------------------------------------------------------------------INH_INT_1_SI                              
                            
    SUBROUTINE      INH_V_QX_INT_1_SI                                  &
                          ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,               TPML,          &
                            UT, WT, QXT, QZT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX,               TPML
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                     UT, WT, QXT, QZT
    END SUBROUTINE  INH_V_QX_INT_1_SI
                            
    SUBROUTINE      INH_W_QZ_INT_1_SI                                  &
                          ( MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,               TPML,          &
                            JMT,                                       &
                            UT, WT, QXT, QZT )
                            
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                               X_MIN, X_MAX,               TPML
      INTEGER(IK),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                       UT, WT, QXT, QZT
     END SUBROUTINE  INH_W_QZ_INT_1_SI

!-------------------------------------------------------------------------------INH_2ND_0,1 
                            
     SUBROUTINE      INH_V_U_2ND_0                                      &
                       (   MXL, MXH, PL, MXT, PXH, MX1,                &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT, QXPT,                                  &
                           REL1_UT, REL2_UT, REL3_UT, REL3_QXT,        &
                           UT, QXT,                                    &
                           TXXT, TXZT, PREST,                          &
                           RXXLT, RXXHT, SXXLT, SXXHT,                 &
                 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                   &
                 ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH                    &
                        )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                           MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH        ,MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                    REL1_UT, REL2_UT, REL3_UT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
      REAL   (PP),DIMENSION(PL :3  ,        MZL:MZH), INTENT(INOUT) :: &
                                                           RXXLT, SXXLT 
      REAL   (PP),DIMENSION(MX1:PXH,        MZL:MZH), INTENT(INOUT) :: &
                                                           RXXHT, SXXHT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                          QXPT, UT, QXT
      REAL   (PP),DIMENSION(0:TPML+1,        MZL:MZH),INTENT(INOUT) :: &
                                   ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, &
                                   ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_V_U_2ND_0
                        
    SUBROUTINE      INH_QX_2ND_0                                           &
                           (   MXL, MXH, PL, MXT, PXH, MX1,                &
                           MZL, MZH,  TPML,                                &
                           X_MIN, X_MAX,                                   &
                           JMT,                                            &
                           REL1_QXT, REL2_QXT, REL3_QXT,                   &
                           QXT,                                            &
                           TXXT, TXZT, PREST,                              &
                           GXXLT, GXXHT, HXXLT, HXXHT,                     &
     ODA2_XL,ODB2_XL,ODC2_XL,ODD2_XL, ODA2_XH,ODB2_XH,ODC2_XH,ODD2_XH   )
                       
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                           MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH        ,MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                           REL1_QXT, REL2_QXT, REL3_QXT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TXZT, PREST
      REAL   (PP), DIMENSION (PL :3  ,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GXXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                  GXXHT
      REAL   (PP), DIMENSION (PL :3  ,      MZL:MZH), INTENT(INOUT) :: &
                                                                  HXXLT
      REAL   (PP), DIMENSION (MX1:PXH,      MZL:MZH), INTENT(INOUT) :: &
                                                                  HXXHT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                                    QXT
      REAL   (PP), DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  INH_QX_2ND_0
                           
    SUBROUTINE      INH_QZ_2ND_0                                           &
                           (   MXL, MXH, PL, MXT, PXH, MX1,                &
                           MZL, MZH,  TPML,                                &
                           X_MIN, X_MAX,                                   &
                           JMT,                                            &
                           REL1_QZT, REL2_QZT, REL3_QZT,                   &
                           QZT,                                            &
                           TZZT, PREST)
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                           MXL, MXH, PL, MXT, MX1, PXH,&
                                           MZL, MZH,  TPML,            &
                                           X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH        ,MZL:MZH), INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                ), INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                      TZZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,        MZL:MZH), INTENT(INOUT) :: &
                                                                    QZT
    END SUBROUTINE  INH_QZ_2ND_0
    
    SUBROUTINE      INH_V_W_2ND_0                                      &
                         ( MXL, MXH,     MXT,                          &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           JMT, QZPT,                                  &
                           REL1_WT, REL2_WT, REL3_WT, REL3_QZT,        &
                           WT, QZT,                                    &
                           TZZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      USE CONTROL_DATA, ONLY: N_FREQ
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                                    MXL, MXH,           MXT,      TPML,&
                                    MZL, MZH,                          &
                                    X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TZZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QZPT, WT, QZT
    END SUBROUTINE  INH_V_W_2ND_0

                           
    SUBROUTINE      INH_V_QZ_INT_1                                     &
                       (    MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT,                                       & 
                            REL1_QZT, REL2_QZT, REL3_QZT,              &
                            QZT,                                       &
                            TZZT, TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                           REL1_QZT, REL2_QZT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TZZT, TXZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    QZT
    END SUBROUTINE  INH_V_QZ_INT_1
                            
    SUBROUTINE      INH_V_W_INT_1                                      &
                       (    MXL, MXH, MXT,                MZL, MZH,    &
                            X_MIN, X_MAX,                              &
                            JMT, QZPT,                                 & 
                            REL1_WT, REL2_WT, REL3_WT, REL3_QZT,       &
                            WT, QZT,                                   &
                            TZZT, TXZT, PREST )
      USE PRECISION   , ONLY: PP, IK
      USE GRID_MEDIUM , ONLY: JMNUM
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                               MXL, MXH, MXT,                MZL, MZH, &
                            X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                    JMT
      REAL   (PP),DIMENSION(1:JMNUM                 ),INTENT(INOUT) :: &
                                    REL1_WT, REL2_WT, REL3_WT, REL3_QZT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                      TZZT, TXZT, PREST
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                          QZPT, WT, QZT
    END SUBROUTINE  INH_V_W_INT_1

    SUBROUTINE  STO_QXP                                                &
                ( L,       MXL, MXH,                                   &
                           MXT,  MX1,                                  &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           QXT, QXPT)
                
      USE PRECISION   , ONLY:     PP, IK
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                            L, MXL, MXH, MXT, MX1, MZL, MZH,  TPML,    &
                            X_MIN, X_MAX
      REAL   (PP), DIMENSION (MXL:MXH,      MZL:MZH), INTENT(INOUT) :: &
                                                             QXPT, QXT
    END SUBROUTINE  STO_QXP

    SUBROUTINE  STO_QZP                                                &
                ( L,       MXL, MXH,                                   &
                           MXT,  MX1,                                  &
                           MZL, MZH,  TPML,                            &
                           X_MIN, X_MAX,                               &
                           QZT, QZPT)
                
      USE PRECISION   , ONLY:     PP, IK
      IMPLICIT NONE
      INTEGER    ,                                    INTENT(IN)    :: &
                            L, MXL, MXH, MXT, MX1, MZL, MZH,  TPML,    &
                            X_MIN, X_MAX
      REAL   (PP), DIMENSION (MXL:MXH,      MZL:MZH), INTENT(INOUT) :: &
                                                             QZPT, QZT
    END SUBROUTINE  STO_QZP
    
    SUBROUTINE      OPEN_SR_FILE
    END SUBROUTINE  OPEN_SR_FILE

    SUBROUTINE      PPML
    END SUBROUTINE  PPML

    SUBROUTINE      PRE_NR_COEF1
    END SUBROUTINE  PRE_NR_COEF1

    SUBROUTINE      READ_MODEL
    END SUBROUTINE  READ_MODEL

    SUBROUTINE      READ_MODEL_PAR
    END SUBROUTINE  READ_MODEL_PAR

    SUBROUTINE      READ_OPEN         ( WORK_NAME )
      IMPLICIT NONE
      CHARACTER (LEN=*), INTENT (IN) :: WORK_NAME
    END SUBROUTINE  READ_OPEN

    SUBROUTINE      STO_BOR                                            &
                    ( MXL, MXH,           PXH,      I1,                &
                      X_MIN, X_MAX,               T,                   &
                      SX_T,       SX_TP,        TPML )
      USE PRECISION   , ONLY: PP
      USE NONREF_BOUND, ONLY: NR_DATA_X
      IMPLICIT NONE
      INTEGER                                  , INTENT (IN)    ::     &
                                     MXL, MXH,           I1,     TPML, &
                                     PXH,                              &
                                     X_MIN, X_MAX
      TYPE(NR_DATA_X),                           INTENT (INOUT) ::     &
                                                            SX_T, SX_TP
      REAL(KIND=PP)  ,DIMENSION(MXL:MXH        ),INTENT (INOUT) :: T
    END SUBROUTINE  STO_BOR

    SUBROUTINE      STORE_SEISMOGRAMS( ITILE )
      INTEGER, INTENT(IN) :: ITILE
    END SUBROUTINE  STORE_SEISMOGRAMS

    SUBROUTINE      STORE_SEISMOGRAMS_PRESSURE( ITILE )
      INTEGER, INTENT(IN) :: ITILE
    END SUBROUTINE  STORE_SEISMOGRAMS_PRESSURE

    SUBROUTINE      STORE_SNAPSHOTS (ITIST)
      INTEGER, INTENT (IN) ::        ITIST
    END SUBROUTINE  STORE_SNAPSHOTS
!-------------------------------------------------------------------------------STRAINS
    SUBROUTINE      STRAIN_UW_XX                                          &
                      (L, MXL, MXH, PL, MXT, PXH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          A, B, C, TPML,                               &
                          UTM, EXXT,                                   &
                          !XXXT,                                       &
                          PXXLT, PXXHT,                                &
                          ODA_XL, ODB_XL, ODC_XL, ODD_XL,              &
                          ODA_XH, ODB_XH, ODC_XH, ODD_XH               &
                      )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, TPML,      &
                                           X_MIN, X_MAX
      REAL   (PP),                                    INTENT (IN)   :: &
                                                               A, B, C
      REAL   (PP),DIMENSION(MXL:MXH                 ),INTENT(INOUT) :: &
                                                                  EXXT
      REAL   (PP),DIMENSION(PL :  1,         MZL:MZH),INTENT(INOUT) :: &
                                                                 PXXLT
      REAL   (PP),DIMENSION(MXT:PXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                 PXXHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                   UTM
      !REAL   (PP),DIMENSION(N_FREQ, MXL:MXH, MZL:MZH),INTENT(INOUT) :: &
                                                                  !XXXT
      REAL   (PP),DIMENSION( 0:TPML,         MZL:MZH),INTENT(INOUT) :: &
        ODA_XL, ODB_XL, ODC_XL, ODD_XL, ODA_XH, ODB_XH, ODC_XH, ODD_XH
    END SUBROUTINE  STRAIN_UW_XX
                          
    SUBROUTINE      STRAIN_Q_XX                                        &
                      (L, MXL, MXH, PL, MXT, PXH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          A, B, C, TPML,                               &
                          QXT, QXXT,                                   & 
                          OXXLT, OXXHT,                                &
                          ODA_XL, ODB_XL, ODC_XL, ODD_XL,              &
                          ODA_XH, ODB_XH, ODC_XH, ODD_XH               &
                          )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, TPML,      &
                                           X_MIN, X_MAX
      REAL   (PP),                                    INTENT (IN)   :: &
                                                               A, B, C
      REAL   (PP),DIMENSION(MXL:MXH                 ),INTENT(INOUT) :: &
                                                                  QXXT
      REAL   (PP), DIMENSION (PL :  1,       MZL:MZH),INTENT(INOUT) :: &
                                                                 OXXLT
      REAL   (PP), DIMENSION (MXT:PXH,       MZL:MZH),INTENT(INOUT) :: &
                                                                 OXXHT
      REAL   (PP),DIMENSION(MXL:MXH,         MZL:MZH),INTENT(INOUT) :: &
                                                                   QXT
      REAL   (PP), DIMENSION ( 0:TPML,       MZL:MZH),INTENT(INOUT) :: &
        ODA_XL, ODB_XL, ODC_XL, ODD_XL, ODA_XH, ODB_XH, ODC_XH, ODD_XH
    
    END SUBROUTINE  STRAIN_Q_XX                      

    SUBROUTINE      STRAIN_UW_XZ                                       &
                     ( L, MXL, MXH, PL, MXT, PXH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          A, B, C, TPML,                               &
                          UTM, WTM, EXZT,                              &
                          !XXZT,                                        &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH           &
                          )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT
      REAL(PP),DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  STRAIN_UW_XZ

    SUBROUTINE      STRAIN_UW_XZ_1                                        &
                     (    MXL, MXH, PL, MXT, PXH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          A, B, C, TPML,                               &
                          UTM, WTM, EXZT,                              &
                          !XXZT,                                        &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH           &
                     )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH,    TPML,      &
                                            X_MIN, X_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT
      REAL(PP),DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
    END SUBROUTINE  STRAIN_UW_XZ_1

    SUBROUTINE      STRAIN_UW_XZ_2ND_PML                               &
                             ( L, LL, C, TPML,                         &
                          MXL, MXH, PL, MXT, PXH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          UTM, WTM,  EXZT,                             &
                          !,XXZT,
                          PZXLT, PZXHT, PXZHT,                         &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,          &
                          ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH           &
                            )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, LL, TPML,  &
                                           X_MIN, X_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                               C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT
      REAL(PP),DIMENSION (PL :  2,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT
      REAL(PP),DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT
      REAL(PP),DIMENSION (MXL:MXH,         0:TPML+1), INTENT(INOUT) :: &
                                                                 PXZHT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT
      REAL(PP),DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
                  ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,                  &
                  ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH
      REAL(PP),DIMENSION(0:TPML+1, MXL:MXH         ), INTENT(INOUT) :: &
                  ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH
    END SUBROUTINE  STRAIN_UW_XZ_2ND_PML

    SUBROUTINE      STRAIN_UW_ZZ                                       &
                     ( L, MXL, MXH, MXT,                MZL, MZH,      &
                          X_MIN, X_MAX,                                &
                          TPML, A, B, C,                               &
                          !XZZT,                                        &
                          WTM, EZZT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                 MXL, MXH, MXT,                TPML, L,&
                                 MZL, MZH,                             &
                                 X_MIN, X_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XZZT
      END SUBROUTINE  STRAIN_UW_ZZ
                          
    SUBROUTINE      STRAIN_Q_ZZ                                        &
                     ( L, MXL, MXH, MXT,                MZL, MZH,      &
                          X_MIN, X_MAX,                                &
                          TPML, A, B, C,                               &
                          QZT, QZZT )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                                 L, MXL, MXH, MXT, TPML,               &
                                 MZL, MZH,                             &
                                 X_MIN, X_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                               A, B, C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  QZZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT
      
      END SUBROUTINE  STRAIN_Q_ZZ                  

    SUBROUTINE      STRAIN_UW_ZZ_0                                        &
                        ( MXL, MXH, MXT,                MZL, MZH,      &
                          X_MIN, X_MAX,               TPML, C,         &
                          WTM, EZZT                                    &
                         !,XZZT 
                        )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                          MXL, MXH, MXT,                MZL, MZH, TPML,&
                          X_MIN, X_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                                     C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH         , MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XZZT
     END SUBROUTINE   STRAIN_UW_ZZ_0
                        
     SUBROUTINE      STRAIN_Q_ZZ_0                                     &
                        ( MXL, MXH, MXT,                MZL, MZH,      &
                          X_MIN, X_MAX,               TPML, C,         &
                          QZT, QZZT)

      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                       INTENT (IN)    :: &
                          MXL, MXH, MXT,                MZL, MZH, TPML,&
                          X_MIN, X_MAX
      REAL(PP),                                      INTENT (IN)    :: &
                                                                     C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  QZZT
      REAL(PP),DIMENSION (MXL:MXH         , MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT

     END SUBROUTINE   STRAIN_Q_ZZ_0                   

    SUBROUTINE      STRAIN_UW_ZZ_2ND_PML                               &
                             ( L, LL, C, TPML,                         &
                          MXL, MXH,     MXT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          WTM,  EZZT,                                  &
                         !XZZT,                                        &                      
                          PZZHT,                                       &
                          ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH               &
                             )
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                            X_MIN, X_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                                     C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EZZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   WTM
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XZZT
      REAL(PP),DIMENSION (MXL:MXH,           0:TPML), INTENT(INOUT) :: &
                                                                 PZZHT
      REAL(PP),DIMENSION ( 0:TPML, MXL:MXH         ), INTENT(INOUT) :: &
                                        ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH
    END SUBROUTINE  STRAIN_UW_ZZ_2ND_PML
                          
    SUBROUTINE      STRAIN_Q_ZZ_2ND_PML                                  &
                             ( L, LL, C, TPML,                         &
                          MXL, MXH,     MXT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          QZT,  QZZT,                                  &
                          OZZHT,                                       &
                          ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH               &
                                )
                        
      USE PRECISION   , ONLY: PP
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                            X_MIN, X_MAX
      REAL(PP),                                       INTENT(IN)    :: &
                                                                     C
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  QZZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT
      REAL   (PP), DIMENSION (MXL:MXH,      0:TPML), INTENT(INOUT) ::  &
                                                                 OZZHT
  REAL   (PP), DIMENSION ( 0:TPML,          MXL:MXH), INTENT(INOUT) :: &
                                        ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH

    END SUBROUTINE  STRAIN_Q_ZZ_2ND_PML                      

    SUBROUTINE      STRESS_NN                                          &
                      (L    , MXL, MXH, MXT,                MZL, MZH,  &
                     X_MIN, X_MAX,                                     &
                      TXXT,          TZZT,          PREST,             &
                      EXXT,          EZZT,                             &
                      QXXT,          QZZT,                             &
                            JMT, TPML )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                       MXL, MXH, MXT,                MZL, MZH, L, TPML,&
                       X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,          MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                      EXXT, EZZT
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                      QXXT, QZZT
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                      !XXXT, XZZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TZZT, PREST
    END SUBROUTINE  STRESS_NN
                            
    SUBROUTINE      STRESS_LIN_NN                                      &
                      (L    , MXL, MXH, MXT,                MZL, MZH,  &
                     X_MIN, X_MAX,                                     &
                      TXXT,          TZZT,          PREST,             &
                      EXXT,          EZZT,                             &
                      QXXT,          QZZT,                             &
                      JMT,           TPML,     SOURS,     SOURF)
    
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                       MXL, MXH, MXT,                MZL, MZH, L, TPML,&
                       X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,          MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                      EXXT, EZZT
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                      QXXT, QZZT
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                      !XXXT, XZZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                      TXXT, TZZT, PREST
      REAL(PP),                                       INTENT(IN)    :: &
                                                      SOURS, SOURF
      
    END SUBROUTINE  STRESS_LIN_NN                      

    SUBROUTINE      STRESS_XZ                                          &
                      (L    , MXL, MXH, MXT,                MZL, MZH,  &
                                          X_MIN, X_MAX,                &
                                          TXZT, EXZT, JMT, TPML        &
                                          !,XXZT
                                           )
      USE PRECISION   , ONLY: PP, IK
      IMPLICIT NONE
      INTEGER,                                        INTENT(IN)    :: &
                       MXL, MXH, MXT,                MZL, MZH, L, TPML,&
                       X_MIN, X_MAX
      INTEGER(IK),DIMENSION(MXL:MXH,          MZL:MZH),INTENT(INOUT):: &
                                                                   JMT
      REAL(PP),DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT
      !REAL(PP),DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT
      REAL(PP),DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                  TXZT
    END SUBROUTINE  STRESS_XZ                     

    SUBROUTINE      TIME_LOOP
    END SUBROUTINE  TIME_LOOP
    
    SUBROUTINE CALCUL_VELOCITY ( SIG_1, P,    AINV, MU,  &
                                 REL1,  REL2, REL3,      &
                                 VPF ,   VPS, VS         )
      USE PRECISION   , ONLY:   PP
    
      IMPLICIT NONE
      
      REAL (KIND=PP), INTENT (IN) :: SIG_1, P, AINV, MU
      REAL (KIND=PP), INTENT (IN) :: REL1, REL2, REL3
      REAL (KIND=PP), INTENT (OUT) :: VPF, VPS, VS
      
    END SUBROUTINE CALCUL_VELOCITY

  END INTERFACE

END MODULE  INTERFACES
