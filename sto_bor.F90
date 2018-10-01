!=======================================================================

SUBROUTINE  STO_BOR ( MXL, MXH,           PXH,      I1,                &
                      X_MIN, X_MAX,               T,                   &
                      SX_T,       SX_TP,        TPML )                     

  USE PRECISION   , ONLY: PP
  USE NONREF_BOUND, ONLY: NR_DATA_X

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER                                      , INTENT (IN)    ::     &
                                     MXL, MXH,           I1,     TPML, &
                                     PXH,                              &
                                     X_MIN, X_MAX
  TYPE (NR_DATA_X),                              INTENT (INOUT) ::     &
                                                            SX_T, SX_TP
  REAL (KIND=PP)  , DIMENSION (MXL:MXH        ), INTENT (INOUT) :: T


  INTEGER :: PXH1, PXH2, IT1, IT2, IT3

!-----------------------------------------------------------------------
  
  !PXH = MX + TPML
  PXH1 = PXH - 1
  PXH2 = PXH - 2
  IT1  = I1  - TPML
  IT2  = IT1 + 1
  IT3  = IT1 + 2

  SX_TP = SX_T

  IF ( ( X_MIN <= IT1 ) .AND. ( IT1 <= X_MAX ) ) THEN
    SX_T%T1KL   = T ( IT1 )
    SX_T%T2KL   = T ( IT2 )
    SX_T%T3KL   = T ( IT3 )
  END IF

  IF ( X_MAX == PXH ) THEN
    SX_T%TMX0KL = T ( PXH  )
    SX_T%TMX1KL = T ( PXH1 )
    SX_T%TMX2KL = T ( PXH2 )
  END IF

END SUBROUTINE STO_BOR
