!=======================================================================

SUBROUTINE PPML                                                                

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MX,     MZ,  DT,                             &
                          TPML, R, TAU, POW, SF, WC, KTBO
  USE GRID_MEDIUM , ONLY:                                              &
                          JM,                                          &
                          MU,                                          &
                          SIG_XX_1, SIG_ZZ_2, PRES_1, PRES_2, PRES_3,  &
                          REL1_U,   REL2_QX,  REL3_U,                  &
                          REL1_W,   REL2_QZ,  REL3_W,                  &
                          H

  USE PML         , ONLY:                                              &
            ODA_XL ,  ODB_XL ,  ODC_XL ,  ODA_XH ,  ODB_XH ,  ODC_XH , &
            ODD_XL ,                      ODD_XH ,                     &
                                          ODA_ZH ,  ODB_ZH ,  ODC_ZH , &
                                          ODD_ZH ,                     &
           ODA2_XL , ODB2_XL , ODC2_XL , ODA2_XH , ODB2_XH , ODC2_XH , &
           ODD2_XL ,                     ODD2_XH ,                     &
                                         ODA2_ZH , ODB2_ZH , ODC2_ZH , &
                                         ODD2_ZH ,                     &
                                        PZXL , PZXH ,                  &   
                                        PXXL , PXXH ,                  &
                                        OXXL , OXXH ,                  &
                                               PXZH ,                  &
                                               PZZH ,                  &  
                                               OZZH ,                  &
                                          RXXL , RXXH ,  SXXL,  SXXH,  &     
                                          RZXL , RZXH ,                &
                                          GXXL , GXXH ,  HXXL,  HXXH,  &
                                          GZXL , GZXH ,                &
                                          RXZH ,                       &
                                          GXZH ,                       &
                                          RZZH , SZZH,                 &
                                          GZZH , HZZH

  USE AUXIL       , ONLY:                                              &
                          PMLL, PMLXH,        PMLZH,                   &
                           X_MIN_D, X_MAX_D, X_MIN, X_MAX,             &
                           Z_MIN_D, Z_MAX_D, Z_MIN, Z_MAX
  
  USE INTERFACES  , ONLY: CALCUL_VELOCITY

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: I, L, J, JM1

  REAL (KIND=PP) :: VEL, VELS, VELPF, VELPS, OMG, OMGPDT, AP, DENOM

!-----------------------------------------------------------------------

IF ( Z_MIN_D <= Z_MAX_D ) THEN
  ALLOCATE (ODA_XL   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODB_XL   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODC_XL   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODD_XL   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODA_XH   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODB_XH   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODC_XH   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODD_XH   (0:TPML   ,                    Z_MIN_D :Z_MAX_D ),&
            ODA_ZH   (0:TPML   , X_MIN_D :X_MAX_D                    ),&
            ODB_ZH   (0:TPML   , X_MIN_D :X_MAX_D                    ),&
            ODC_ZH   (0:TPML   , X_MIN_D :X_MAX_D                    ),&
            ODD_ZH   (0:TPML   , X_MIN_D :X_MAX_D                    ),&
            ODA2_XL  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODB2_XL  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODC2_XL  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODD2_XL  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODA2_XH  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODB2_XH  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODC2_XH  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODD2_XH  (0:TPML +1,                    Z_MIN_D :Z_MAX_D ),&
            ODA2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D                    ),&
            ODB2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D                    ),&
            ODC2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D                    ),&
            ODD2_ZH  (0:TPML +1, X_MIN_D :X_MAX_D                    ) )

  ALLOCATE (  PZXL (PMLL   :2      ,                  Z_MIN_D:Z_MAX_D),&
              PZXH (MX     :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
              PXXL (PMLL   :1      ,                  Z_MIN_D:Z_MAX_D),&
              OXXL (PMLL   :1      ,                  Z_MIN_D:Z_MAX_D),&
              PXXH (MX     :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
              OXXH (MX     :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
              PXZH (X_MIN_D:X_MAX_D,                        0:TPML+1 ),&
              PZZH (X_MIN_D:X_MAX_D,                        0:TPML   ),&
              OZZH (X_MIN_D:X_MAX_D,                        0:TPML   ) )

  ALLOCATE (   RXXL(PMLL   :3      ,                  Z_MIN_D:Z_MAX_D),&
               RXXH(MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
               SXXL(PMLL   :3      ,                  Z_MIN_D:Z_MAX_D),&
               SXXH(MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
               GXXL(PMLL   :3      ,                  Z_MIN_D:Z_MAX_D),&
               GXXH(MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
               HXXL(PMLL   :3      ,                  Z_MIN_D:Z_MAX_D),&
               HXXH(MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
               RZXL(PMLL   :2      ,                  Z_MIN_D:Z_MAX_D),&
               RZXH(MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
               GZXL(PMLL   :2      ,                  Z_MIN_D:Z_MAX_D),&
               GZXH(MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D),&
               RXZH(X_MIN_D:X_MAX_D,                  0      :TPML   ),&
               GXZH(X_MIN_D:X_MAX_D,                  0      :TPML   ),&
               RZZH(X_MIN_D:X_MAX_D,                  0      :TPML   ),&
               SZZH(X_MIN_D:X_MAX_D,                  0      :TPML   ),&
               GZZH(X_MIN_D:X_MAX_D,                  0      :TPML   ),&
               HZZH(X_MIN_D:X_MAX_D,                  0      :TPML   ) )

  !ZETA
  PZXL (PMLL   :2      ,                  Z_MIN_D:Z_MAX_D) = 0._PP 
  PZXH (MX     :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP 
  
  PXXL (PMLL   :1      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  PXXH (MX     :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  !KSI
  OXXL (PMLL   :1      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  OXXH (MX     :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  !ZETA
  PXZH (X_MIN_D:X_MAX_D,                  0      :TPML+1 ) = 0._PP 
  PZZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  !KSI
  OZZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP

  !THETA
  RXXL (PMLL   :3      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  RXXH (MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  !PSI
  SXXL (PMLL   :3      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  SXXH (MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  
  GXXL (PMLL   :3      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  GXXH (MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  !PHI
  HXXL (PMLL   :3      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  HXXH (MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  !THETA
  RZXL (PMLL   :2      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  RZXH (MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  
  GZXL (PMLL   :2      ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  GZXH (MX-1   :PMLXH  ,                  Z_MIN_D:Z_MAX_D) = 0._PP
  
  RXZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  GXZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  RZZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  SZZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  GZZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  HZZH (X_MIN_D:X_MAX_D,                  0      :TPML   ) = 0._PP
  

  IF ( TPML == 0 ) THEN
    OMG = 0._PP
  ELSE
    OMG = -TAU*LOG(R)*DT/REAL(TPML)/(REAL(TPML)**POW)/H 
  END IF

  DO  L = Z_MIN, Z_MAX

!______________________________________________________ LEFT TRANSPARENT

    IF ( X_MIN <= 1 ) THEN
        JM1 = JM ( X_MIN, L )
        CALL CALCUL_VELOCITY( SIG_XX_1(JM1), - PRES_1(JM1), - PRES_3(JM1), MU(JM1), H*REL1_U(JM1)/DT, H*REL2_QX(JM1)/DT, H*REL3_U(JM1)/DT, VELPF , VELPS, VELS )
        VEL  = VELS

        ODA_XL   (:,  L) = 1._PP
        ODA2_XL  (:,  L) = 1._PP

        ODB_XL   (:,  L) = 1._PP
        ODB2_XL  (:,  L) = 1._PP

        ODC_XL   (:,  L) = 0._PP
        ODC2_XL  (:,  L) = 0._PP

        DO J = 1, TPML
            
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91   
          DENOM = EXP(-OMGPDT/AP-WC*DT)     
          ODA_XL  (J,  L) = (DENOM+1.)/2./AP
          ODB_XL  (J,  L) =  DENOM
          ODC_XL  (J,  L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XL (J,  L) = (DENOM+1.)/2./AP
          ODB2_XL (J,  L) =  DENOM
          ODC2_XL (J,  L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XL (J,  L) = (DENOM+1.)/2./AP
          ODB2_XL (J,  L) =  DENOM
          ODC2_XL (J,  L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        ODD_XL   (:,  L) = (ODC_XL  (:,  L)+2.)/(ODB_XL  (:,  L)+1.)
        ODD2_XL  (:,  L) = (ODC2_XL (:,  L)+2.)/(ODB2_XL (:,  L)+1.)

    END IF

!_____________________________________________________ RIGHT TRANSPARENT

    IF ( X_MAX >= MX ) THEN
        JM1 = JM(X_MAX ,  L)
        CALL CALCUL_VELOCITY( SIG_XX_1(JM1), - PRES_1(JM1), - PRES_3(JM1), MU(JM1), H*REL1_U(JM1)/DT, H*REL2_QX(JM1)/DT, H*REL3_U(JM1)/DT, VELPF , VELPS, VELS )
        VEL  = VELS

        ODA_XH   (:,  L) = 1._PP
        ODA2_XH  (:,  L) = 1._PP

        ODB_XH   (:,  L) = 1._PP
        ODB2_XH  (:,  L) = 1._PP

        ODC_XH   (:,  L) = 0._PP
        ODC2_XH  (:,  L) = 0._PP

        DO J = 1, TPML
          OMGPDT = OMG*VEL* REAL(J)     **POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA_XH  (J,  L) = (DENOM+1.)/2./AP
          ODB_XH  (J,  L) =  DENOM
          ODC_XH  (J,  L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XH (J,  L) = (DENOM+1.)/2./AP
          ODB2_XH (J,  L) =  DENOM
          ODC2_XH (J,  L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        END DO

        J = TPML +1
          OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
          AP     = 1.+SF*OMGPDT*0.075/0.91
          DENOM = EXP(-OMGPDT/AP-WC*DT)
          ODA2_XH (J,  L) = (DENOM+1.)/2./AP
          ODB2_XH (J,  L) =  DENOM
          ODC2_XH (J,  L) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        ODD_XH   (:,  L) = (ODC_XH  (:,  L)+2.)/(ODB_XH  (:,  L)+1.)
        ODD2_XH  (:,  L) = (ODC2_XH (:,  L)+2.)/(ODB2_XH (:,  L)+1.)

    END IF

  END DO

!____________________________________________________ BOTTOM TRANSPARENT

  IF ( Z_MAX >= MZ ) THEN
    IF ( KTBO /= 0 ) OMG = 0.  
    DO I = X_MIN, X_MAX
      JM1 = JM ( I, Z_MAX )
      CALL CALCUL_VELOCITY( SIG_ZZ_2(JM1), - PRES_2(JM1), - PRES_3(JM1), MU(JM1), H*REL1_W(JM1)/DT, H*REL2_QZ(JM1)/DT, H*REL3_W(JM1)/DT, VELPF , VELPS, VELS )
      VEL  = VELPF

      ODA_ZH   (:,I) = 1._PP
      ODA2_ZH  (:,I) = 1._PP

      ODB_ZH   (:,I) = 1._PP
      ODB2_ZH  (:,I) = 1._PP

      ODC_ZH   (:,I) = 0._PP
      ODC2_ZH  (:,I) = 0._PP

      DO J = 1, TPML
        OMGPDT = OMG*VEL* REAL(J)     **POW
        AP     = 1.+SF*OMGPDT*0.075/0.91
        DENOM = EXP(-OMGPDT/AP-WC*DT)
        ODA_ZH  (J,I) = (DENOM+1.)/2./AP
        ODB_ZH  (J,I) =  DENOM
        ODC_ZH  (J,I) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

        OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
        AP     = 1.+SF*OMGPDT*0.075/0.91
        DENOM = EXP(-OMGPDT/AP-WC*DT)
        ODA2_ZH (J,I) = (DENOM+1.)/2./AP
        ODB2_ZH (J,I) =  DENOM
        ODC2_ZH (J,I) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

      END DO

      J = TPML +1
        OMGPDT = OMG*VEL*(REAL(J)-0.5)**POW
        AP     = 1.+SF*OMGPDT*0.075/0.91
        DENOM  = EXP(-OMGPDT/AP-WC*DT)
        ODA2_ZH (J,I) = (DENOM+1.)/2./AP
        ODB2_ZH (J,I) =  DENOM
        ODC2_ZH (J,I) =  OMGPDT / (OMGPDT+WC*AP*DT) * (DENOM-1.)

      ODD_ZH   (:,I) = (ODC_ZH  (:,I)+2.)/(ODB_ZH  (:,I)+1.)
      ODD2_ZH  (:,I) = (ODC2_ZH (:,I)+2.)/(ODB2_ZH (:,I)+1.)
    END DO

  END IF

END IF


END SUBROUTINE PPML
