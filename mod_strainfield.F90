!_______________________________________________________________________

MODULE STRAINFIELD
  USE PRECISION
  IMPLICIT NONE

  REAL    (KIND=PP), DIMENSION (:), ALLOCATABLE ::                   &
                                        EXX , EXZ , EZZ,             &
                                        QXX ,       QZZ
!$SGI DISTRIBUTE EXZ  (BLOCK)
!$SGI DISTRIBUTE EXX  (BLOCK)
!$SGI DISTRIBUTE EZZ  (BLOCK)

END MODULE STRAINFIELD
