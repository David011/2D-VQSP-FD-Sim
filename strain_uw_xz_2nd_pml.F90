!=======================================================================

SUBROUTINE STRAIN_UW_XZ_2ND_PML ( L, LL, C, TPML,                      &
                          MXL, MXH, PL, MXT, PXH,                      &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          UTM, WTM, EXZT,                              &
                           !,XXZT,
                          PZXLT, PZXHT, PXZHT,                         &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH,          &
                          ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH           &
                                )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, LL, TPML,  &
                                           X_MIN, X_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                               C

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT

  REAL   (PP), DIMENSION (PL :  2,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT

  REAL   (PP), DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT

  REAL   (PP), DIMENSION (MXL:MXH,         0:TPML+1), INTENT(INOUT) :: &
                                                                 PXZHT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM

  !REAL   (PP), DIMENSION (N_FREQ,  MXL:MXH, MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT

  REAL   (PP), DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

  REAL   (PP), DIMENSION(0:TPML+1, MXL:MXH         ), INTENT(INOUT) :: &
                                     ODA2_ZH, ODB2_ZH, ODC2_ZH, ODD2_ZH

!-----------------------------------------------------------------------

  INTEGER       :: I, IO, LP1 , L1  , L2  , MXT1

  REAL(PP)      :: DZXA, DXZA

!-----------------------------------------------------------------------

LP1  = L  + 1
L1   = L  - 1
L2   = L  - 2

MXT1  = MXT - 1

!================================== EXZT ===============================
  DO  I = MAX(2-TPML,X_MIN), MIN(2,X_MAX)

    DZXA          = C*(   WTM (I  ,     L  ) - WTM (I-1,     L ) )
    DXZA          = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )
    
    !EXZT(I  )     = DZXA + DXZA

    EXZT(I  )     = ODA2_ZH(LL ,I  ) * ( ODD2_ZH(LL ,I  )*DXZA         &
                                       + PXZHT (I,  LL)            )   &
                  + ODA2_XL(3-I,  L) * ( ODD2_XL(3-I,  L)*DZXA         &
                                       + PZXLT (I,  L )            )

    PXZHT(I,   LL)= ODB2_ZH(LL ,I  ) * PXZHT(I,  LL)                   &      
                  + ODC2_ZH(LL ,I  ) * DXZA

    PZXLT(I  , L )= ODB2_XL(3-I,  L) * PZXLT(I,  L)                    &       
                  + ODC2_XL(3-I,  L) * DZXA

  END DO
  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)                                         

    DZXA          = C*(   WTM (I  ,     L  ) - WTM (I-1,     L ) )
    DXZA          = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )
    
    !EXZT(I  ) = DZXA + DXZA
  
    EXZT(I  )     = ODA2_ZH(LL ,I  ) * ( ODD2_ZH(LL ,I  )*DXZA         &      
                                       + PXZHT (I  ,LL)            )   &
                  +                      DZXA

    PXZHT(I,   LL)= ODB2_ZH(LL ,I  ) * PXZHT(I,  LL)                   &     
                  + ODC2_ZH(LL ,I  ) * DXZA

  END DO
  DO  I = MAX(MXT,X_MIN), MIN(MXT+TPML,X_MAX)

    DZXA          = C*(   WTM (I  ,     L  ) - WTM (I-1,     L ) )
    DXZA          = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )
    
    !EXZT(I  )      = DZXA + DXZA

    EXZT(I  )     = ODA2_ZH(LL    ,I  ) * ( ODD2_ZH(LL    ,I  )*DXZA   &
                                          + PXZHT (I,  LL)           ) &
                  + ODA2_XH(I-MXT1,  L) * ( ODD2_XH(I-MXT1,  L)*DZXA   &
                                          + PZXHT (I,  L )           )

    PXZHT(I,   LL)= ODB2_ZH(LL ,I  ) * PXZHT(I,  LL)                   &        
                  + ODC2_ZH(LL ,I  ) * DXZA

    PZXHT(I,   L )= ODB2_XH(I-MXT1,  L) * PZXHT(I,  L)                 &           
                  + ODC2_XH(I-MXT1,  L) * DZXA

  END DO

END SUBROUTINE  STRAIN_UW_XZ_2ND_PML
