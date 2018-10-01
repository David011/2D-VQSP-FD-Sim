!-----------------------------------------------------------------------

SUBROUTINE  NONREF (H, KEY, VEL, VELS, VELP,                           &
                    AX01, AX02, AX10, AX11, AX12, AX20, AX21, AX22)
  USE PRECISION   , ONLY:                                              &
                          PP, DP

  USE CONTROL_DATA, ONLY:                                              &
                          DT, OMG, WB,                                 &
                          THPPLE, THPPRI,                 THPPBO,      &
                          THSSLE, THSSRI,                 THSSBO

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: KEY

  REAL (KIND=PP), INTENT (IN)  :: H, VEL, VELS, VELP

  REAL (KIND=PP), INTENT (OUT) ::                                      &
                         AX01, AX02, AX10, AX11, AX12, AX20, AX21, AX22

  REAL (KIND=DP)               ::  gplus, gminus, dtaux, haux,         &
                         aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22

  REAL (KIND=PP)               :: TMP, TMP1, TMP2

  INTERFACE
    SUBROUTINE  HIGDON_FD (WB, ALFA, BETA, DT, H,                      &
                           A00, A01, A02, A10, A11, A12, A20, A21, A22)
      USE PRECISION
      IMPLICIT NONE
      REAL (KIND=PP), INTENT (IN)  ::                                  &
                        WB, ALFA, BETA, DT, H
      REAL (KIND=DP), INTENT (OUT) ::                                  &
                        A00, A01, A02, A10, A11, A12, A20, A21, A22
    END SUBROUTINE HIGDON_FD

    SUBROUTINE  OptAbsBc(gplus, gminus, dt, dx,                        &
                     a00, a01, a02, a10, a11, a12, a20, a21, a22)
      USE PRECISION, ONLY: DP
        real(DP), INTENT(IN)  :: gplus, gminus, dt, dx
        real(DP), INTENT(OUT) :: a00,a01,a02,a10,a11,a12,a20,a21,a22
    END SUBROUTINE OptAbsBc

    SUBROUTINE  HIGCE_FD  (WB, ALFA, BETA, DT, H,                      &
                           A00, A01, A02, A10, A11, A12, A20, A21, A22)
      USE PRECISION
      IMPLICIT NONE
      REAL (PP), INTENT (IN)  :: WB, ALFA, BETA, DT, H
      REAL (DP), INTENT (OUT) :: A00,A01,A02,A10,A11,A12,A20,A21,A22
    END SUBROUTINE HIGCE_FD

    SUBROUTINE  HIGCE2_FD  (WB, VELH, VELCE, DT, H,                    &
                            A00, A01, A02, A10, A11, A12, A20, A21, A22)
      USE PRECISION
      IMPLICIT NONE
      REAL (KIND=PP), INTENT (IN)  ::                                  &
                        WB, VELH, VELCE, DT, H
      REAL (KIND=DP), INTENT (OUT) ::                                  &
                        A00, A01, A02, A10, A11, A12, A20, A21, A22
    END SUBROUTINE HIGCE2_FD

  END INTERFACE

!-----------------------------------------------------------------------

    SELECT CASE ( KEY )
!                             Solid Boundary
    CASE ( 0 )
      AX01 =  0._PP
      AX02 =  0._PP
      AX10 =  0._PP
      AX11 =  0._PP
      AX12 =  0._PP
      AX20 =  0._PP
      AX21 =  0._PP
      AX22 =  0._PP

!                             Higdon
    CASE ( 1 )
      CALL HIGDON_FD ( WB, VELP, VELS, DT, H,                          &
                       aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22)
      AX01 = REAL (aa01,PP)
      AX02 = REAL (aa02,PP)
      AX10 = REAL (aa10,PP)
      AX11 = REAL (aa11,PP)
      AX12 = REAL (aa12,PP)
      AX20 = REAL (aa20,PP)
      AX21 = REAL (aa21,PP)
      AX22 = REAL (aa22,PP)

!                             Reynolds
    CASE ( 2 )
      TMP   = -DT*VEL/H
      AX01 =  0._PP
      AX02 =  0._PP
      AX10 =  1._PP + TMP
      AX11 =  1._PP - TMP
      AX12 =  0._PP
      AX20 =  0._PP
      AX21 = -1._PP - TMP
      AX22 =          TMP

!                               Peng & Toksoz
    CASE ( 3 )
      gplus = DBLE (THPPLE)
      gminus= DBLE (THSSLE)
      dtaux = DBLE (OMG*DT)
       haux = DBLE (OMG/VEL*H)
       CALL OptAbsBc(gplus, gminus, dtaux, haux,                        &
                    aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22)
      AX01 = REAL (aa01)
      AX02 = REAL (aa02)
      AX10 = REAL (aa10)
      AX11 = REAL (aa11)
      AX12 = REAL (aa12)
      AX20 = REAL (aa20)
      AX21 = REAL (aa21)
      AX22 = REAL (aa22)

!                               Emerman & Stephen
    CASE ( 4 )
      TMP1 = ( DT - H/VEL )/( DT + H/VEL )
      TMP2 =  2._PP*H/VEL  /( DT + H/VEL )
      AX01 =  TMP1
      AX02 =  0._PP
      AX10 =  TMP2
      AX11 =  TMP2
      AX12 =  0._PP
      AX20 =  TMP1
      AX21 = -1._PP
      AX22 =  0._PP

!                               Clayton & Engquist A1
    CASE ( 5 )
      TMP1 =  VEL*DT/H
      AX01 = - ( 1._PP - TMP1 )/( 1._PP + TMP1 )
      AX02 =  0._PP
      AX10 =   ( 1._PP - TMP1 )/( 1._PP + TMP1 )
      AX11 =  1._PP
      AX12 =  0._PP
      AX20 =  0._PP
      AX21 =  0._PP
      AX22 =  0._PP

!              Higdon (S) + Clayton & Engquist A1Simple(P)
    CASE ( 6 )
      CALL HIGCE_FD ( WB, VELP, VELS, DT, H,                           &
                      aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22)
      AX01 = REAL (aa01,PP)
      AX02 = REAL (aa02,PP)
      AX10 = REAL (aa10,PP)
      AX11 = REAL (aa11,PP)
      AX12 = REAL (aa12,PP)
      AX20 = REAL (aa20,PP)
      AX21 = REAL (aa21,PP)
      AX22 = REAL (aa22,PP)

!                  Higdon  + Clayton & Engquist A1
    CASE ( 7 )
      CALL HIGCE2_FD( WB, VEL , VEL , DT, H,                           &
                      aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22)
      AX01 = REAL (aa01,PP)
      AX02 = REAL (aa02,PP)
      AX10 = REAL (aa10,PP)
      AX11 = REAL (aa11,PP)
      AX12 = REAL (aa12,PP)
      AX20 = REAL (aa20,PP)
      AX21 = REAL (aa21,PP)
      AX22 = REAL (aa22,PP)

!                  Higdon (S)  + Clayton & Engquist A1 (P)
    CASE ( 8 )
      CALL HIGCE2_FD( WB, VELS, VELP, DT, H,                           &
                      aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22)
      AX01 = REAL (aa01,PP)
      AX02 = REAL (aa02,PP)
      AX10 = REAL (aa10,PP)
      AX11 = REAL (aa11,PP)
      AX12 = REAL (aa12,PP)
      AX20 = REAL (aa20,PP)
      AX21 = REAL (aa21,PP)
      AX22 = REAL (aa22,PP)

!                               Peng & Toksoz
    CASE ( 9 )
      gplus = DBLE (THPPLE)
      gminus= DBLE (THSSLE)
      dtaux = DBLE (OMG*DT)
       haux = DBLE (OMG/VELS*H)
       CALL OptAbsBc(gplus, gminus, dtaux, haux,                        &
                    aa00,aa01,aa02,aa10,aa11,aa12,aa20,aa21,aa22)
      AX01 = REAL (aa01)
      AX02 = REAL (aa02)
      AX10 = REAL (aa10)
      AX11 = REAL (aa11)
      AX12 = REAL (aa12)
      AX20 = REAL (aa20)
      AX21 = REAL (aa21)
      AX22 = REAL (aa22)

    END SELECT

END SUBROUTINE NONREF
