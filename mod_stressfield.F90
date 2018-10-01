!_______________________________________________________________________

MODULE STRESSFIELD
  USE PRECISION
  IMPLICIT NONE

  REAL    (KIND=PP), DIMENSION (:,:), ALLOCATABLE ::                   &
                                        TXX , TXZ , TZZ, PRES
!$SGI DISTRIBUTE TXZ  (*,BLOCK)
!$SGI DISTRIBUTE TXX  (*,BLOCK)
!$SGI DISTRIBUTE TZZ  (*,BLOCK)

END MODULE STRESSFIELD
