SUBROUTINE STRAIN_UW_XZ ( L, MXL, MXH, PL, MXT, PXH,                   &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          A, B, C, TPML,                               &
                          UTM, WTM,  EXZT,                             &
                          !XXZT,                                        &
                          PZXLT, PZXHT,                                &
                          ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL,          &
                          ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH           & 
                          )

  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, TPML,      &
                                            X_MIN, X_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXZT

  REAL   (PP), DIMENSION (PL :  2,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXLT

  REAL   (PP), DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PZXHT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                              UTM, WTM

  !REAL   (PP), DIMENSION (N_FREQ, MXL:MXH,  MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXZT

  REAL   (PP), DIMENSION(0:TPML+1,          MZL:MZH), INTENT(INOUT) :: &
 ODA2_XL, ODB2_XL, ODC2_XL, ODD2_XL, ODA2_XH, ODB2_XH, ODC2_XH, ODD2_XH

!-----------------------------------------------------------------------

  INTEGER       :: I, IO, LP1 , L1  , L2  , MXT1

  REAL(PP)      :: DZXA

!-----------------------------------------------------------------------

LP1  = L  + 1
L1   = L  - 1
L2   = L  - 2

MXT1  = MXT - 1

!================================== EXZT ===============================
!----- FIRST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(2-TPML,X_MIN), MIN(2,X_MAX)                                    
                                                                             
    DZXA          = C*(   WTM (I  ,     L  ) - WTM (I-1,     L ) )
    
    !EXZT(I   )    = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )      &
                  !+ DZXA
    

    EXZT(I  )     = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )     &
                  + ODA2_XL(3-I,  L) * ( ODD2_XL(3-I,  L)*DZXA         &
                                       + PZXLT (I,  L)         )

    PZXLT(I,   L )= ODB2_XL(3-I,  L) * PZXLT(I,  L)                    &
                  + ODC2_XL(3-I,  L) * DZXA

  END DO
! INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(3,X_MIN), MIN(MXT1,X_MAX)                                         
                                                                                 
                                                                        
    EXZT(I  )     =  ( A*( UTM (I  ,     LP1 ) - UTM (I  ,     L2  )   &
                         + WTM (I+1,     L   ) - WTM (I-2,     L   ) ) &
                     + B*( UTM (I  ,     L   ) - UTM (I  ,     L1  )   &
                         + WTM (I  ,     L   ) - WTM (I-1,     L   ) ) )

  END DO
  
!----- LAST ROW OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(MXT,X_MIN), MIN(MXT+TPML,X_MAX)

    DZXA          = C*(   WTM (I  ,     L  ) - WTM (I-1,     L ) )
    
    !EXZT(I  )     = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )     &
                  !+ DZXA

    EXZT(I  )     = C*(   UTM (I  ,     L  ) - UTM (I  ,     L1) )     &
                  + ODA2_XH(I-MXT1,  L) * ( ODD2_XH(I-MXT1,  L)*DZXA   &
                                          + PZXHT(I,  L)             )

    PZXHT(I,   L )= ODB2_XH(I-MXT1,  L) * PZXHT(I,  L)                 &
                  + ODC2_XH(I-MXT1,  L) * DZXA

  END DO

END SUBROUTINE  STRAIN_UW_XZ
