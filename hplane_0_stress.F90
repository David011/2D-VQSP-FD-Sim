!=======================================================================

SUBROUTINE  HPLANE_0_STRESS

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          JM
  USE CONTROL_DATA, ONLY:                                              &
                          MX,     TPML, STRESS_IMAGING
  USE AUXIL       , ONLY:                                              &
                          AH , BH , CH , A, B,                         &
                          PMLL, PMLXH,                                 &
                           X_MIN  , X_MAX  ,                           &
                           X_MIN_D, X_MAX_D,                           &
                           Z_MIN_D, Z_MAX_D
 USE WAVEFIELD   , ONLY:                                               &
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
                                          PXXL  , PXXH  ,              &
                                          OXXL  , OXXH
  
  USE INTERFACES  , ONLY:                                              &
                          STRAIN_UW_XX,   STRAIN_Q_XX,                 &
                          STRESS_XZ,      STRESS_NN,                   &
                          STRAIN_UW_ZZ_0, STRAIN_Q_ZZ_0
                          

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER             :: L

!-----------------------------------------------------------------------

IF (STRESS_IMAGING) THEN                                                   
    
  TXZ (:,  0) = -TXZ (:,  4)
  TZZ (:,  0) = -TZZ (:,  3)                                      
  !TZZ (:,  0) = -TZZ (:,  3)-B/A*(-TZZ (:,  2)+1./3.*TZZ (:,  3))
  
  PRES (:,  0) = -PRES (:,  3)                                       
 !PRES (:,  0) = -PRES (:,  3)-B/A*(-PRES (:,  2)+1./3.*PRES (:,  3))
ELSE
  L   = 0

  EXZ (:    ) = 0.  
  !XXZ (:,:,L) = 0.                                                 !artefact

  CALL STRAIN_UW_ZZ_0    ( MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,    &    
                        MZL = Z_MIN_D, MZH = Z_MAX_D,                  &     
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &    
                        TPML = TPML,                   C = CH ,        &    
                        WTM = WM ,             EZZT = EZZ              &
                        !,XZZT = XZZ                                   
                         )
  CALL STRAIN_Q_ZZ_0     ( MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,    &    
                        MZL = Z_MIN_D, MZH = Z_MAX_D,                  &    
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &    
                        TPML = TPML,                   C = CH ,        &
                        QZT = QZ ,             QZZT = QZZ              &
                         )
   
   

 CALL STRAIN_UW_XX(L = 0 ,                                             &    
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

 CALL STRAIN_Q_XX(L = 0 ,                                              &
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
     
  CALL STRESS_XZ (L=0 ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXZT = TXZ , EXZT = EXZ ,                       &
                       JMT = JM , TPML = TPML                          &
                         !,XXZT = XXZ                                  &
                 )

  CALL STRESS_NN (L=0 ,MXL = X_MIN_D, MXH = X_MAX_D, MXT = MX ,        &
                       MZL = Z_MIN_D, MZH = Z_MAX_D,                   &
                        X_MIN = X_MIN , X_MAX = X_MAX ,                &
                       TXXT = TXX , TZZT = TZZ , PREST = PRES,         &
                       EXXT = EXX ,              EZZT = EZZ ,          &
                       !XXXT = XXX ,              XZZT = XZZ ,          &
                       QXXT = QXX ,              QZZT = QZZ ,          &
                       JMT = JM , TPML = TPML  )
END IF

END SUBROUTINE HPLANE_0_STRESS
