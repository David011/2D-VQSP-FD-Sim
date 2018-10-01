SUBROUTINE  STRAIN_Q_XX (L, MXL, MXH, PL, MXT, PXH,                    &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          A, B, C, TPML,                               &
                          QXT, QXXT,                                   &
                          OXXLT, OXXHT,                                &
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
                                                                  QXXT
  
  REAL   (PP), DIMENSION (PL :  1,          MZL:MZH), INTENT(INOUT) :: &
                                                                 OXXLT

  REAL   (PP), DIMENSION (MXT:PXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                 OXXHT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   QXT
  
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

    
    DXXA          = C*( QXT (I+1,     L  ) - QXT (I  ,     L ) )
    
    QXXT(I  )     = ODA_XL(2-I,  L) * ( ODD_XL(2-I,  L)*DXXA           &
                                      + OXXLT(I,  L)         )
    
    OXXLT(I,   L )=  ODB_XL(2-I,  L) * OXXLT(I,  L)                     &
                  + ODC_XL(2-I,  L) * DXXA

  END DO

!----- FIRST COLUMN COMPUTED BY 2ND ORDER FD APPROX.
  IF ( ( X_MIN <= 2 ) .AND. ( 2 <= X_MAX ) ) THEN
    I = 2
    QXXT(I  )     =    C*( QXT (I+1,     L  ) - QXT (I  ,     L ) )
  END IF

!----- INTERIOR POINTS COMPUTED BY 4TH ORDER FD APPROX.
  DO  I = MAX(3,X_MIN), MIN(MXT-2,X_MAX)

    QXXT(I  ) = A*( QXT (I+2,     L  ) - QXT (I-1,     L  ) )          &
              + B*( QXT (I+1,     L  ) - QXT (I  ,     L  ) )            
    

  END DO

!----- LAST COLUMN COMPUTED BY 2ND ORDER FD APPROX.
  IF ( ( X_MIN <= MXT1 ) .AND. ( MXT1 <= X_MAX) ) THEN
    I = MXT1
    QXXT(I  )     =    C*( QXT (I+1,     L  ) - QXT (I      , L ) )
  END IF

! LAST COLUMN OF PML COMPUTED BY 2ND ORDER FD APPROX.
  DO  I = MAX(MXT,X_MIN), MIN(MXT1+TPML,X_MAX)                  
 
    DXXA          = C*( QXT (I+1,     L  ) - QXT (I  ,     L ) )
    
    QXXT(I  )     = ODA_XH(I-MXT1,  L) * ( ODD_XH(I-MXT1,  L)*DXXA     &
                                         + OXXHT(I,  L)            )

    OXXHT(I,   L )= ODB_XH(I-MXT1,  L) * OXXHT(I,  L)                  &
                  + ODC_XH(I-MXT1,  L) * DXXA

  END DO



END SUBROUTINE  STRAIN_Q_XX