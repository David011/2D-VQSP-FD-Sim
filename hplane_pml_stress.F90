SUBROUTINE  HPLANE_PML_STRESS ( L, LL )

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          MX, TPML
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
                                          ODA_ZH ,  ODB_ZH ,  ODC_ZH , &
                                          ODD_ZH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                                         ODD2_ZH ,                     &
                                          PZXL  , PZXH  ,              &
                                          PXXL  , PXXH  ,              &
                                                  PXZH  ,              &
                                          OXXL  , OXXH  ,              &
                                                  PZZH  ,              &
                                                  OZZH
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_UW_XZ_2ND_PML, STRAIN_UW_ZZ_2ND_PML,  &
                          STRAIN_Q_ZZ_2ND_PML,                         &
                          STRAIN_UW_XX, STRAIN_Q_XX,                   &
                          STRESS_XZ, STRESS_NN

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: L, LL

!-----------------------------------------------------------------------

  CALL STRAIN_UW_XZ_2ND_PML  ( L = L , LL = LL  , C = CH , TPML = TPML,&
       MXL = X_MIN_D, MXH = X_MAX_D, PL = PMLL, MXT = MX, PXH = PMLXH, &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        UTM = UM , WTM = WM , EXZT = EXZ,              &
                        !,XXZT = XXZ ,                                 &
                        PZXLT = PZXL , PZXHT = PZXH ,                  &
                                       PXZHT = PXZH ,                  &
       ODA2_XL = ODA2_XL , ODB2_XL = ODB2_XL ,                         &
       ODC2_XL = ODC2_XL , ODD2_XL = ODD2_XL ,                         &
       ODA2_XH = ODA2_XH , ODB2_XH = ODB2_XH ,                         &
       ODC2_XH = ODC2_XH , ODD2_XH = ODD2_XH ,                         &
       ODA2_ZH = ODA2_ZH , ODB2_ZH = ODB2_ZH ,                         &
       ODC2_ZH = ODC2_ZH , ODD2_ZH = ODD2_ZH                           &
                             )
      
      

  CALL STRAIN_UW_ZZ_2ND_PML  ( L = L , LL = LL  , C = CH , TPML = TPML,&
       MXL = X_MIN_D, MXH = X_MAX_D,            MXT = MX,              &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        WTM = WM ,             EZZT = EZZ,             &           
                      !,XZZT = XZZ ,
                                       PZZHT = PZZH ,                  &
                ODA_ZH = ODA_ZH, ODB_ZH = ODB_ZH,                      &
                ODC_ZH = ODC_ZH, ODD_ZH = ODD_ZH                       &
                            )
      
  CALL STRAIN_Q_ZZ_2ND_PML  ( L = L , LL = LL , C = CH , TPML = TPML,  &
       MXL = X_MIN_D, MXH = X_MAX_D,            MXT = MX,              &
       MZL = Z_MIN_D, MZH = Z_MAX_D,                                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                        QZT = QZ , QZZT = QZZ,                         &           
                        OZZHT = OZZH ,                                 &
                ODA_ZH = ODA_ZH, ODB_ZH = ODB_ZH,                      &
                ODC_ZH = ODC_ZH, ODD_ZH = ODD_ZH                       &
                            )

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

   CALL STRESS_XZ (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,       &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXZT = TXZ , EXZT = EXZ ,                       &
                       JMT = JM , TPML = TPML                          &
                         !,XXZT = XXZ                                  &
                         )

  CALL STRESS_NN (L=L ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXXT = TXX , TZZT = TZZ , PREST = PRES,         &
                       EXXT = EXX ,              EZZT = EZZ ,          &
                       !XXXT = XXX ,              XZZT = XZZ ,          &
                       QXXT = QXX ,              QZZT = QZZ ,          &
                       JMT = JM , TPML = TPML  )


END SUBROUTINE HPLANE_PML_STRESS
