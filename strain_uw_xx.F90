!=======================================================================

SUBROUTINE  STRAIN_UW_XX (L, MXL, MXH, PL, MXT, PXH,                   &
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

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH, PL, MXT, PXH, &
                                               MZL, MZH, L, TPML,      &
                                           X_MIN, X_MAX

  REAL   (PP),                                       INTENT (IN)    :: &
                                                               A, B, C

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  EXXT

  REAL   (PP), DIMENSION (PL :  1,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PXXLT

  REAL   (PP), DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 PXXHT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   UTM

  !REAL   (PP), DIMENSION (N_FREQ, MXL:MXH,  MZL:MZH), INTENT(INOUT) :: &
                                                                  !XXXT

  REAL   (PP), DIMENSION ( 0:TPML,          MZL:MZH), INTENT(INOUT) :: &
        ODA_XL, ODB_XL, ODC_XL, ODD_XL, ODA_XH, ODB_XH, ODC_XH, ODD_XH

!-----------------------------------------------------------------------

  INTEGER       :: I, IO, MXT1

  REAL (PP)     :: DXXA

!-----------------------------------------------------------------------

MXT1  = MXT - 1

!================================== EXXT ===============================

! FIRST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(2-TPML,X_MIN), MIN(1,X_MAX)                                    

    DXXA          = C*( UTM (I+1,     L  ) - UTM (I  ,     L ) )
    
    EXXT(I  )     = ODA_XL(2-I,  L) * ( ODD_XL(2-I,  L)*DXXA           &
                                      + PXXLT(I,  L)         )

    PXXLT(I,   L )= ODB_XL(2-I,  L) * PXXLT(I,  L)                     &
                  + ODC_XL(2-I,  L) * DXXA


  END DO

!----- FIRST COLUMN COMPUTED BY 2ND ORDER FD APPROX.
  IF ( ( X_MIN <= 2 ) .AND. ( 2 <= X_MAX ) ) THEN
    I = 2
    EXXT(I  )     =    C*( UTM (I+1,     L  ) - UTM (I  ,     L ) )
  END IF

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    EXXT(I  ) = A*( UTM (I+2,     L  ) - UTM (I-1,     L  ) )          &
              + B*( UTM (I+1,     L  ) - UTM (I  ,     L  ) )           
    

  END DO

!----- LAST COLUMN COMPUTED BY 2ND ORDER FD APPROX.
  IF ( ( X_MIN <= MXT1 ) .AND. ( MXT1 <= X_MAX) ) THEN
    I = MXT1
    EXXT(I  )     =    C*( UTM (I+1,     L  ) - UTM (I      , L ) )
  END IF

! LAST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)      

    DXXA          = C*( UTM (I+1,     L  ) - UTM (I  ,     L ) )


    EXXT(I  )     = ODA_XH(I-MXT1,  L) * ( ODD_XH(I-MXT1,  L)*DXXA     &
                                         + PXXHT(I,  L)            )

    PXXHT(I,   L )= ODB_XH(I-MXT1,  L) * PXXHT(I,  L)                  &
                  + ODC_XH(I-MXT1,  L) * DXXA


  END DO

END SUBROUTINE  STRAIN_UW_XX
