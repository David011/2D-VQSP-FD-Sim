!=======================================================================

SUBROUTINE PRE_NR_COEF1

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MZ,  DT,                                     &
                          KTLE, KTRI,             KTBO
  USE AUXIL       , ONLY:                                              &
                           PMLL, PMLXH ,         PMLZH,                &
                                             X_MIN, X_MAX,             &
                                             Z_MIN, Z_MAX
  USE GRID_MEDIUM , ONLY:                                              &
                          JM,                                          &
                          MU,                                          &
                          SIG_XX_1, SIG_ZZ_2, PRES_1, PRES_2, PRES_3,  &
                          REL1_U,   REL2_QX,  REL3_U,                  &
                          REL1_W,   REL2_QZ,  REL3_W,                  &
                          H
  USE NONREF_BOUND, ONLY:                                              &
                          AXU,                AXW,                     &
                          ABW01,   ABW02,   ABW10,   ABW11,            &
                          ABW12,   ABW20,   ABW21,   ABW22,            &
                          ABH01,   ABH02,   ABH10,   ABH11,            &
                          ABH12,   ABH20,   ABH21,   ABH22
  
  USE INTERFACES  , ONLY: CALCUL_VELOCITY

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: L, JM1, JM2, JM3, JM4

  REAL (KIND=PP) :: VEL, VELS, VELPF, VELPS, KT, K1, K2

!-----------------------------------------------------------------------

  INTERFACE

    SUBROUTINE     NONREF (H, KEY, VEL, VELS, VELP,                    &
                           AX01,AX02,AX10,AX11,AX12,AX20,AX21,AX22)
      USE PRECISION, ONLY: PP
      INTEGER,        INTENT (IN)  :: KEY
      REAL (KIND=PP), INTENT (IN)  :: H, VEL, VELS, VELP
      REAL (KIND=PP), INTENT (OUT) ::                                  &
                           AX01,AX02,AX10,AX11,AX12,AX20,AX21,AX22
    END SUBROUTINE NONREF

  END INTERFACE

!-----------------------------------------------------------------------
!                                                   HETEROGENEOUS MEDIUM
  DO  L = MAX(Z_MIN,0), MIN(Z_MAX,MZ-1)

!______________________________________________________ LEFT TRANSPARENT

   IF ( X_MIN == PMLL ) THEN
     IF ( KTLE >=0 ) THEN

      JM1 = JM ( PMLL  , L )
      JM2 = JM ( PMLL+1, L )
      JM4 = JM ( PMLL+1, L+1 )
      
      CALL CALCUL_VELOCITY( SIG_XX_1(JM2), - PRES_1(JM2), - PRES_3(JM2), MU(JM2), H*REL1_U(JM2)/DT, H*REL2_QX(JM2)/DT, H*REL3_U(JM2)/DT, VELPF , VELPS, VELS ) 
      VEL  = VELPF
      CALL NONREF ( H, KTLE, VEL, VELS, VELPF,                          &
                                AXU(  L)%ALE01, AXU(  L)%ALE02,        &
                AXU(  L)%ALE10, AXU(  L)%ALE11, AXU(  L)%ALE12,        &
                AXU(  L)%ALE20, AXU(  L)%ALE21, AXU(  L)%ALE22   )

      IF ( L > 0 ) THEN
        JM3 = JM ( PMLL, L-1 )
      ELSE
        JM3 = JM1  
      END IF

      K1  = SIG_ZZ_2(JM1)
      K2  = SIG_ZZ_2(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      CALL CALCUL_VELOCITY( SIG_XX_1(JM1), - PRES_1(JM1), - PRES_3(JM1), MU(JM1), H*REL1_U(JM1)/DT, H*REL2_QX(JM1)/DT, H*REL3_U(JM1)/DT, VELPF , VELPS, VELS ) 
      VEL  = VELS
      CALL NONREF ( H, KTLE, VEL, VELS, VELPF,                          &
                                AXW(  L)%ALE01, AXW(  L)%ALE02,        &
                AXW(  L)%ALE10, AXW(  L)%ALE11, AXW(  L)%ALE12,        &
                AXW(  L)%ALE20, AXW(  L)%ALE21, AXW(  L)%ALE22   )

     ELSE

      AXU (  L)%ALE01 = -1.
      AXU (  L)%ALE02 =  0.
      AXU (  L)%ALE10 =  0.
      AXU (  L)%ALE11 =  0.
      AXU (  L)%ALE12 =  0.
      AXU (  L)%ALE20 =  0.
      AXU (  L)%ALE21 =  0.
      AXU (  L)%ALE22 =  0.

      AXW (  L)%ALE01 =  0.
      AXW (  L)%ALE02 =  1.
      AXW (  L)%ALE10 =  0.
      AXW (  L)%ALE11 =  0.
      AXW (  L)%ALE12 =  0.
      AXW (  L)%ALE20 =  0.
      AXW (  L)%ALE21 =  0.
      AXW (  L)%ALE22 =  0.

     END IF

    IF ( Z_MAX >= MZ ) THEN
        AXU(MZ:Z_MAX) = AXU(MZ-1) 
        AXW(MZ:Z_MAX) = AXW(MZ-1)
    END IF
   END IF

!_____________________________________________________ RIGHT TRANSPARENT

   IF ( X_MAX == PMLXH ) THEN

     IF ( KTRI >=0 ) THEN

      JM1 = JM(PMLXH  ,L)
      JM2 = JM(PMLXH-1,L)

      JM3 = JM2
      JM4 = JM( PMLXH, L+1 )

      CALL CALCUL_VELOCITY( SIG_XX_1(JM2), - PRES_1(JM2), - PRES_3(JM2), MU(JM2), H*REL1_U(JM2)/DT, H*REL2_QX(JM2)/DT, H*REL3_U(JM2)/DT, VELPF , VELPS, VELS )
      VEL  = VELPF
      CALL NONREF ( H, KTRI, VEL, VELS, VELPF,                          &
                                AXU(  L)%ARI01, AXU(  L)%ARI02,        &
                AXU(  L)%ARI10, AXU(  L)%ARI11, AXU(  L)%ARI12,        &
                AXU(  L)%ARI20, AXU(  L)%ARI21, AXU(  L)%ARI22   )

      JM3 = JM1

      IF ( L > 0 ) THEN
        JM3 = JM ( PMLXH , L-1 )
      ELSE
        JM3 = JM1
      END IF

      K1  = SIG_ZZ_2(JM1)
      K2  = SIG_ZZ_2(JM3)
      KT  = 2._PP*K1*K2/(K1+K2)

      CALL CALCUL_VELOCITY( SIG_XX_1(JM1), - PRES_1(JM1), - PRES_3(JM1), MU(JM1), H*REL1_U(JM1)/DT, H*REL2_QX(JM1)/DT, H*REL3_U(JM1)/DT, VELPF , VELPS, VELS ) 
      VEL  = VELS
      CALL NONREF ( H, KTRI, VEL, VELS, VELPF,                          &
                                AXW(  L)%ARI01, AXW(  L)%ARI02,        &
                AXW(  L)%ARI10, AXW(  L)%ARI11, AXW(  L)%ARI12,        &
                AXW(  L)%ARI20, AXW(  L)%ARI21, AXW(  L)%ARI22   )

     ELSE

      AXU (  L)%ARI01 = -1.
      AXU (  L)%ARI02 =  0.
      AXU (  L)%ARI10 =  0.
      AXU (  L)%ARI11 =  0.
      AXU (  L)%ARI12 =  0.
      AXU (  L)%ARI20 =  0.
      AXU (  L)%ARI21 =  0.
      AXU (  L)%ARI22 =  0.

      AXW (  L)%ARI01 =  0.
      AXW (  L)%ARI02 =  1.
      AXW (  L)%ARI10 =  0.
      AXW (  L)%ARI11 =  0.
      AXW (  L)%ARI12 =  0.
      AXW (  L)%ARI20 =  0.
      AXW (  L)%ARI21 =  0.
      AXW (  L)%ARI22 =  0.

     END IF

    IF ( Z_MAX >= MZ ) THEN
        AXU(  MZ:Z_MAX) = AXU(  MZ-1)
        AXW(  MZ:Z_MAX) = AXW(  MZ-1)
    END IF
   END IF
  END DO

!____________________________________________________ BOTTOM TRANSPARENT

  IF ( Z_MAX >= MZ ) THEN
    JM1 = JM ( X_MIN, MZ   )
    JM2 = JM ( X_MIN, MZ-1 )
    CALL CALCUL_VELOCITY( SIG_XX_1(JM1), - PRES_1(JM1), - PRES_3(JM1), MU(JM1), H*REL1_U(JM1)/DT, H*REL2_QX(JM1)/DT, H*REL3_U(JM1)/DT, VELPF , VELPS, VELS )
    VEL  = VELPF
    CALL NONREF ( H,KTBO, VEL, VELS, VELPF,                             &
                          ABW01, ABW02,                                &
                   ABW10, ABW11, ABW12,                                &
                   ABW20, ABW21, ABW22  )

    VEL  = VELS
    CALL NONREF ( H,KTBO, VEL, VELS, VELPF,                             &
                          ABH01, ABH02,                                &
                   ABH10, ABH11, ABH12,                                &
                   ABH20, ABH21, ABH22  )
  END IF


END SUBROUTINE PRE_NR_COEF1
