!-----------------------------------------------------------------------

SUBROUTINE  HIGCE_FD  (WB, ALFA, BETA, DT, H,                          &
                       A00, A01, A02, A10, A11, A12, A20, A21, A22)
  USE PRECISION

  IMPLICIT NONE

  REAL (PP), INTENT (IN)  :: WB, ALFA, BETA, DT, H

  REAL (DP), INTENT (OUT) :: A00, A01, A02, A10, A11, A12, A20, A21, A22

  REAL (DP)               :: H1X, H1T, HXT, C1T, CXT, GP, GS, GS1, HC


  GP = ALFA*DT/H
  GS = BETA*DT/H

  GS1 = 1._DP + GS

  H1X = ( GS    - WB*GS1 ) / ( GS1*(1-WB) )
  H1T = ( 1._DP - WB*GS1 ) / ( GS1*(1-WB) )
  HXT = WB / ( 1._DP - WB )

  C1T = 1._DP - GP
  CXT = GP

  HC = H1T * C1T


  A00 =     0._DP

  A01 =     H1X
  A02 =     0._DP
  A10 =         C1T + H1T
  A11 =         CXT + HXT    - H1X*C1T
  A12 =   - H1X*CXT
  A20 =   - H1T*C1T
  A21 =   - HXT*C1T - H1T*CXT
  A22 =   - HXT*CXT

END SUBROUTINE  HIGCE_FD
