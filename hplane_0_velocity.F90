!=======================================================================

SUBROUTINE  HPLANE_0_VELOCITY

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM ,                                         &
                          REL1_U,         REL2_U,         REL3_U,      &
                          REL1_W,         REL2_W,         REL3_W,      &
                          REL1_QX,        REL2_QX,        REL3_QX,     &
                          REL1_QZ,        REL2_QZ,        REL3_QZ
  USE CONTROL_DATA, ONLY:                                              &
                          MX, TPML, STRESS_IMAGING
  USE AUXIL       , ONLY:                                              &
                               MX1 ,                                   &
                          PMLL, PMLXH,                                 &
                           X_MIN  , X_MAX  ,                           &
                           X_MIN_D, X_MAX_D,                           &
                           Z_MIN_D, Z_MAX_D
  USE WAVEFIELD   , ONLY:                                              &
                          UM ,      WM ,                               &
                          QX,       QZ ,                               &
                          QXP,      QZP

  USE STRESSFIELD , ONLY:                                              &
                          TXX , TZZ , TXZ, PRES
  USE NONREF_BOUND, ONLY:                                              &
                          AXU , AXW ,                                  &
                                    SX_U   , SX_UP  ,                  &
                                    SX_W   , SX_WP  ,                  &
                                    SX_QX  , SX_QXP ,                  &
                                    SX_QZ  , SX_QZP 
  USE PML         , ONLY:                                              &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
                               RXXL  , RXXH  ,  SXXL,  SXXH,           &
                               RZXL  , RZXH  ,                         &
                               GXXL  , GXXH  ,  HXXL,  HXXH,           &
                               GZXL  , GZXH
  USE INTERFACES  , ONLY:                                              &
                          INH_QX_INT_0,     INH_QZ_INT_0   ,           &
                          INH_V_UV_INT_0,   INH_V_W_INT_0  ,           &
                          INH_QX_2ND_0 ,                               &           
                          INH_QZ_2ND_0,                                &
                          INH_V_U_2ND_0 ,   INH_V_W_2ND_0  ,           &
                          STO_QXP,          STO_QZP        ,           &
                          STO_BOR

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: L

!-----------------------------------------------------------------------

IF (STRESS_IMAGING) THEN
ELSE
  L   = 0
      
  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     PXH = PMLXH ,                                     &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                     SX_T = SX_U  (  L ),                              &
                     SX_TP= SX_UP       ,                              &
                     T = UM (:,  L ), I1 = 2,         TPML = TPML )

  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     PXH = PMLXH ,                                     &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                     SX_T = SX_W  (  L ),                              &
                     SX_TP= SX_WP       ,                              &
                     T = WM (:,  L ), I1 = 1,         TPML = TPML )
  
  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     PXH = PMLXH ,                                     &     
                        X_MIN = X_MIN , X_MAX = X_MAX ,                & 
                     SX_T = SX_QX  (  L ), SX_TP= SX_QXP,              &     
                     T = QX (:,L ), I1 = 2, TPML = TPML )                
                                                                           
  CALL  STO_BOR    ( MXL = X_MIN_D, MXH = X_MAX_D,                     &
                     PXH = PMLXH ,                                     &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                     SX_T = SX_QZ  (  L ), SX_TP= SX_QZP,              &
                     T = QZ (:,L ), I1 = 1, TPML = TPML )

  CALL  INH_QX_INT_0                                                                &
                   (   MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,                     &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                                &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                              &
                       JMT = JM ,                                                   &
                       REL1_QXT = REL1_QX, REL2_QXT = REL2_QX, REL3_QXT = REL3_QX,  &
                       QXT = QX,                                                    &
                       TXXT = TXX , TXZT = TXZ, PREST = PRES )

  CALL  INH_V_UV_INT_0                                                        &
                   (   MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,               &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                          &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                        &
                       JMT = JM , QXPT = QXP,                                 &
                       REL1_UT = REL1_U, REL2_UT = REL2_U, REL3_UT = REL3_U,  &
                       REL3_QXT = REL3_QX,                                    &
                       UT = UM, QXT = QX,                                     &
                       TXXT = TXX , TXZT = TXZ, PREST = PRES  )

  CALL  INH_QZ_INT_0                                                                &
                   (   MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,                     &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                                &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                              &
                       JMT = JM ,                                                   &
                       REL1_QZT = REL1_QZ, REL2_QZT = REL2_QZ, REL3_QZT = REL3_QZ,  &
                       QZT = QZ,                                                    &
                       TZZT = TZZ, TXZT = TXZ , PREST = PRES  )
  
  CALL  INH_V_W_INT_0                                                         &
                   (   MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,               &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                          &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                       &
                       JMT = JM , QZPT = QZP,                                 &
                       REL1_WT = REL1_W, REL2_WT = REL2_W, REL3_WT = REL3_W,  &
                       REL3_QZT = REL3_QZ,                                    &
                       WT = WM, QZT = QZ,                                     &
                       TZZT = TZZ , TXZT = TXZ, PREST = PRES )

  CALL  INH_QX_2ND_0                                                                 &
               (       MXL = X_MIN_D, MXH = X_MAX_D, PL  = PMLL,                     &
                       MXT = MX, MX1 = MX1,          PXH = PMLXH,                    &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,                   &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                               &
                       JMT = JM ,                                                    &
                       REL1_QXT = REL1_QX, REL2_QXT = REL2_QX, REL3_QXT = REL3_QX,   &
                       QXT = QX,                                                     &
                       TXXT = TXX , TXZT = TXZ , PREST = PRES,                       &
                       GXXLT = GXXL , GXXHT = GXXH , HXXLT = HXXL , HXXHT = HXXH ,   &
           ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                                   &
           ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                                   &
           ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                                   &
           ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH                                     & 
              )
                   
  CALL  INH_V_U_2ND_0                                                                &
               (       MXL = X_MIN_D, MXH = X_MAX_D, PL  = PMLL,                     &
                       MXT = MX, MX1 = MX1,          PXH = PMLXH,                    &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,                   &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                               &
                       JMT = JM , QXPT = QXP,                                        &
                       REL1_UT = REL1_U, REL2_UT = REL2_U, REL3_UT = REL3_U,         &
                       REL3_QXT = REL3_QX,                                           &
                       UT = UM, QXT = QX,                                            &
                       TXXT = TXX , TXZT = TXZ , PREST = PRES,                       &
                       RXXLT = RXXL , RXXHT = RXXH , SXXLT = SXXL , SXXHT = SXXH ,   &
           ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                                   &
           ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                                   &
           ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                                   &
           ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH                                     &
               )
      
   CALL  INH_QZ_2ND_0                                                               &
               (       MXL = X_MIN_D, MXH = X_MAX_D, PL  = PMLL,                    &
                       MXT = MX, MX1 = MX1,          PXH = PMLXH,                   &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,                  &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                              &
                       JMT = JM ,                                                   &
                       REL1_QZT = REL1_QZ, REL2_QZT = REL2_QZ, REL3_QZT = REL3_QZ,  &
                       QZT = QZ,                                                    &
                       TZZT = TZZ , PREST = PRES )
      
  CALL  INH_V_W_2ND_0                                                         &
                   (   MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,               &
                       MZL = Z_MIN_D, MZH = Z_MAX_D, TPML = TPML,             &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                       &
                       JMT = JM , QZPT = QZP,                                 &
                       REL1_WT = REL1_W, REL2_WT = REL2_W, REL3_WT = REL3_W,  &
                       REL3_QZT = REL3_QZ,                                    &
                       WT = WM, QZT = QZ,                                     &
                       TZZT = TZZ, PREST = PRES  )
END IF

  CALL  STO_QXP                                                                &
              ( L = L, MXL = X_MIN_D, MXH = X_MAX_D,                           &
                       MXT = MX, MX1 = MX1,                                    &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,             &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                         &
                       QXT = QX, QXPT = QXP)
              
  CALL  STO_QZP                                                                &
              ( L = L, MXL = X_MIN_D, MXH = X_MAX_D,                           &
                       MXT = MX, MX1 = MX1,                                    &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,  TPML = TPML,             &
                       X_MIN = X_MIN , X_MAX = X_MAX ,                         &
                       QZT = QZ, QZPT = QZP)

END SUBROUTINE HPLANE_0_VELOCITY

