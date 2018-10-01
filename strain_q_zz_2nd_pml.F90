SUBROUTINE STRAIN_Q_ZZ_2ND_PML ( L, LL, C, TPML,                       &
                          MXL, MXH,     MXT,                           &
                          MZL, MZH,                                    &
                          X_MIN, X_MAX,                                &
                          QZT,  QZZT,                                  &
                          OZZHT,                                       &
                          ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH               &
                                )


  USE PRECISION   , ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE


  INTEGER,                                            INTENT(IN)    :: &
                                               MXL, MXH,     MXT,      &
                                               MZL, MZH, L, LL, TPML,  &
                                            X_MIN, X_MAX

  REAL   (PP),                                        INTENT(IN)    :: &
                                                                     C

  REAL   (PP), DIMENSION (MXL:MXH                  ), INTENT(INOUT) :: &
                                                                  QZZT

  REAL   (PP), DIMENSION (MXL:MXH,          MZL:MZH), INTENT(INOUT) :: &
                                                                   QZT
  
  REAL   (PP), DIMENSION (MXL:MXH,           0:TPML), INTENT(INOUT) :: &
                                                                 OZZHT

  REAL   (PP), DIMENSION ( 0:TPML,          MXL:MXH), INTENT(INOUT) :: &
                                        ODA_ZH, ODB_ZH, ODC_ZH, ODD_ZH


!-----------------------------------------------------------------------

  INTEGER       :: I, LP1, IO, MXT1

  REAL (PP)     :: DZZA

!-----------------------------------------------------------------------

  LP1    = L + 1

  MXT1  = MXT - 1

  DO  I = MAX(2-TPML,X_MIN), MIN(MXT1+TPML,X_MAX)

    

    DZZA          = C*(   QZT (I  ,     LP1 ) - QZT (I  ,     L ) )
    
    QZZT(I  )     = ODA_ZH(LL,I  ) * ( ODD_ZH(LL,I  )*DZZA           &
                                     + OZZHT(I,  LL)       )

    OZZHT(I,  LL) = ODB_ZH(LL,I  ) * OZZHT(I,  LL)                   &
                  + ODC_ZH(LL,I  ) * DZZA
    
    
  END DO


END SUBROUTINE  STRAIN_Q_ZZ_2ND_PML