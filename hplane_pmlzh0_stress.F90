!=======================================================================

SUBROUTINE  HPLANE_PMLZH0_STRESS

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MZ, TPML
  USE AUXIL       , ONLY:                                              &
                          AH , BH , CH ,                               &
                          PMLL, PMLXH,                                 &
                           X_MIN  , X_MAX  ,                           &
                           X_MIN_D, X_MAX_D,                           &
                           Z_MIN_D, Z_MAX_D
  USE WAVEFIELD   , ONLY:                                              &
                          UM ,      WM ,                               &
                          QX,       QZ                             
  USE STRESSFIELD , ONLY:                                              &
                          TXX , TZZ , TXZ, PRES
  USE STRAINFIELD , ONLY:                                              &
                          EXX , EZZ , EXZ,                             &
                          QXX , QZZ
  USE PML         , ONLY:                                              &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                                         ODD2_ZH ,                     &
                                          PZXL  , PZXH  ,              &
                                                  PXZH
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_UW_XZ_2ND_PML, STRESS_XZ

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: L, LL

!-----------------------------------------------------------------------

  L  = MZ+TPML
  LL = TPML+1

  CALL STRAIN_UW_XZ_2ND_PML  ( L = L , LL = LL  , C = CH , TPML = TPML,&
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        UTM = UM , WTM = WM , EXZT = EXZ,              &
                        !XXZT = XXZ ,                                 &
                        PZXLT = PZXL , PZXHT = PZXH ,                  &
                                       PXZHT = PXZH ,                  &
       ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                         &
       ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                         &
       ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                         &
       ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH ,                         &
       ODA2_ZH = ODA2_ZH , ODB2_ZH = ODB2_ZH ,                         &
       ODC2_ZH = ODC2_ZH , ODD2_ZH = ODD2_ZH                           &
                             )

  CALL STRESS_XZ (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXZT = TXZ , EXZT = EXZ ,                       &
                       JMT = JM , TPML = TPML                          &
                         !,XXZT = XXZ                                  &
                 )


END SUBROUTINE HPLANE_PMLZH0_STRESS
