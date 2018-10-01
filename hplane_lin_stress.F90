!=======================================================================

SUBROUTINE  HPLANE_LIN_STRESS ( SOURTF_LIN, SOURTF_LOU )

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          MX, TPML, LIN_Z, LIN_X
  USE AUXIL       , ONLY:                                              &
                          AH , BH , CH ,                               &
                               MX1 ,                                   &
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
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
                                          PZXL  , PZXH  ,              &
                                          PXXL  , PXXH  ,              &
                                          OXXL  , OXXH
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_UW_XZ, STRAIN_UW_ZZ, STRAIN_UW_XX,    &
                                        STRAIN_Q_ZZ , STRAIN_Q_XX,     &
                          STRESS_XZ, STRESS_LIN_NN

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: L
  REAL(PP), INTENT(IN) :: SOURTF_LIN, SOURTF_LOU 
  REAL(PP), DIMENSION (X_MIN_D:X_MAX_D) :: TXZ_TEMP, TXX_TEMP, TZZ_TEMP
  !REAL(PP), DIMENSION (N_FREQ, X_MIN_D:X_MAX_D) :: XXZ_TEMP, XXX_TEMP, XZZ_TEMP
  REAL(PP) :: A, B

!-----------------------------------------------------------------------
  L = LIN_Z


   CALL STRAIN_UW_XZ(L = L ,                                           &
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
       A = AH, B = BH, C = CH, TPML = TPML,                            &
                        UTM = UM , WTM = WM  , EXZT = EXZ,             &
                        !XXZT = XXZ,                                    &
                        PZXLT = PZXL , PZXHT = PZXH ,                  &
       ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                         &
       ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                         &
       ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                         &
       ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH                           &
                    )


  CALL STRAIN_UW_ZZ(L =L , MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,    &
                        MZL = Z_MIN_D, MZH = Z_MAX_D,                  &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        TPML = TPML, A = AH , B = BH , C = CH ,        &
                        !XZZT = XZZ ,                                   &
                        WTM = WM ,             EZZT = EZZ  )
  
  CALL STRAIN_Q_ZZ(L =L , MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,     &
                        MZL = Z_MIN_D, MZH = Z_MAX_D,                  &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        TPML = TPML, A = AH , B = BH , C = CH ,        &
                        QZT = QZ ,             QZZT = QZZ  )

  CALL STRAIN_UW_XX(L = L ,                                            &
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
       A = AH, B = BH, C = CH, TPML = TPML,                            &
                        UTM = UM ,             EXXT = EXX,             &
                        !XXXT = XXX ,                                   &
                        PXXLT = PXXL , PXXHT = PXXH ,                  &
             ODA_XL = ODA_XL , ODB_XL = ODB_XL ,                       &
             ODC_XL = ODC_XL , ODD_XL = ODD_XL ,                       &
             ODA_XH = ODA_XH , ODB_XH = ODB_XH ,                       &
             ODC_XH = ODC_XH , ODD_XH = ODD_XH                         &
                   )
      
  CALL STRAIN_Q_XX(L = L ,                                             &
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
       A = AH, B = BH, C = CH, TPML = TPML,                            &
                        QXT = QX ,             QXXT = QXX,             &
                        OXXLT = OXXL , OXXHT = OXXH ,                  &
             ODA_XL = ODA_XL , ODB_XL = ODB_XL ,                       &
             ODC_XL = ODC_XL , ODD_XL = ODD_XL ,                       &
             ODA_XH = ODA_XH , ODB_XH = ODB_XH ,                       &
             ODC_XH = ODC_XH , ODD_XH = ODD_XH                         &
                   )


  CALL STRESS_XZ (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXZT = TXZ , EXZT = EXZ ,                       &
                       JMT = JM , TPML = TPML                          &
                         !,XXZT = XXZ                                  &
                         )

  CALL STRESS_LIN_NN (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,    &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXXT = TXX , TZZT = TZZ , PREST = PRES,         &
                       EXXT = EXX ,              EZZT = EZZ ,          &
                       !XXXT = XXX ,              XZZT = XZZ ,          &
                       QXXT = QXX ,              QZZT = QZZ ,          &
                       JMT = JM , TPML = TPML,                         &
                       SOURS = SOURTF_LIN, SOURF = SOURTF_LOU)


END SUBROUTINE HPLANE_LIN_STRESS
