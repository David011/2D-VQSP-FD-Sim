!-----------------------------------------------------------------------

SUBROUTINE  OptAbsBc(gplus, gminus, dt, dx,                            &
                     a00, a01, a02, a10, a11, a12, a20, a21, a22)
!======================================================================
! ====================================================================
!      Geophysical Center For Parallel Computer Application
!              Earth  Resources Laboratory
!          Massachusetts Institute Of Technology
!                        OPTABSBC
!
!                    Copyright @1992 of
!           Massachusetts Institute Of Technology
!
!                   written by CHENGBIN PENG
! ====================================================================
! Design an optimal absorbing condition
! inputs:
! gplus : theta_plus,  incidence angle in radiant with zero reflection
! gminus: theta_minus, incidence angle in radiant with zero reflection
!     dt: dimensionless time step
!     dx: dimensionless grid size
! outputs:
! a_00 -> a_22: nine coeffficients of the optimal absorbing boundary
!
USE PRECISION, ONLY: DP
  INTEGER :: i, j

  real(DP), INTENT(IN)  :: gplus, gminus, dt, dx
  real(DP), INTENT(OUT) :: a00, a01, a02, a10, a11, a12, a20, a21, a22
  real(DP) :: tplus, tminus
  real(DP) :: Rz, Iz, Rx, Ix, Ry, Iy
  real(DP) :: ratio, temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
  real(DP) :: phi, psi, delta
  real(DP) :: r0(0:2,0:2), r1(0:2,0:2)

  tplus =  dx * dcos( gplus )
  tminus = dx * dcos( gminus)
  Rz = dcos( dt )
  Iz = dsin( dt )
  Rx = dcos( tplus + tminus )
  Ix = dsin( tplus + tminus )
  Ry = -( dcos( tplus ) + dcos( tminus ) )
  Iy = -( dsin( tplus ) + dsin( tminus ) )


  temp1 = Rx * Rz - Ix * Iz
  temp2 = Rx * Iz + Rz * Ix
  temp3 = Ry * Rz - Iy * Iz
  temp4 = Ry * Iz + Rz * Iy

  delta = ( Ry + Iy + 2.0* Rz) * ( temp1 + temp2 -Rz + Iz )          &
        + ( 1.0 - Rx - Ix    ) * ( 2.0 + temp3 + temp4    )
  delta = 1.0/delta

  r1(1,0) = (1.0 - Rz) *( temp1 + temp2 -Rz + Iz ) * delta
  r1(0,2) = (1.0 - Rz) *( 1.0 - Rx - Ix    ) * delta
  r1(0,1) = -0.50 - r1(1,0)
  r1(1,1) = 1.0 - 2.0 * r1(0,2)
  r1(0,0) = 0.0
  r1(2,2) = 0.0
  r1(1,2) = r1(1,0)
  r1(2,0) = r1(0,2)
  r1(2,1) = r1(0,1)

  temp5 = Rx * Rz + Ix * Iz
  temp6 = Rx * Iz - Rz * Ix
  temp7 = Ry * Rz + Iy * Iz
  temp8 = Ry * Iz - Rz * Iy
  phi =- ( Rz + Iz - temp5 + temp6  ) * ( Ry + Iy + 2.0 * Rz )       &
       + ( 1.0 + Rz + temp7 - temp8 ) * ( 1.0 - Rx - Ix )
  r0(0,2) = phi * delta

  psi = ( 1.0 + Rz + temp7 - temp8 ) * ( temp1 + temp2 -Rz + Iz )    &
      + ( Rz + Iz - temp5 + temp6  ) * ( 2.0 + temp3 + temp4    )
  r0(1,0) = psi * delta

  r0(0,0) = 0.
  r0(0,1) = 0.50 -     r0(1,0)
  r0(1,1) = 1.0  - 2.0*r0(0,2)
  r0(2,2) = -1.0
  r0(1,2) = r0(1,0)
  r0(2,0) = r0(0,2)
  r0(2,1) = r0(0,1)

  temp1 = 0.0
  temp2 = 0.0
  do i=0, 2
    do j=0, 2
      temp1 = temp1 + r0(i,j)*r1(i,j)
      temp2 = temp2 + r1(i,j)*r1(i,j)
    end do
  end do
  ratio = -temp1/temp2

  a00 = 0.0
  a01 = r0(0,1) + ratio * r1(0,1)
  a02 = r0(0,2) + ratio * r1(0,2)
  a10 = r0(1,0) + ratio * r1(1,0)
  a11 = r0(1,1) + ratio * r1(1,1)
  a12 = r0(1,2) + ratio * r1(1,2)
  a20 = r0(2,0) + ratio * r1(2,0)
  a21 = r0(2,1) + ratio * r1(2,1)
  a22 = r0(2,2) + ratio * r1(2,2)

END SUBROUTINE OptAbsBc
