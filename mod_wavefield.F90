!_______________________________________________________________________

MODULE WAVEFIELD
  USE PRECISION, ONLY: PP
  IMPLICIT NONE

  REAL    (KIND=PP), DIMENSION (:)    , ALLOCATABLE ::                 &
                    SEISU, SEISW, DSEISU, DSEISW, SSU, SSW, SSUF, SSWF
  
  REAL    (KIND=PP), DIMENSION (:)    , ALLOCATABLE ::                 &
                    SEISQX, SEISQZ, DSEISQX, DSEISQZ, SSQX, SSQZ

  REAL    (KIND=PP), DIMENSION (:)    , ALLOCATABLE ::                 &
                    SEISP, SEIST

  REAL    (KIND=PP), DIMENSION (:,:)  , ALLOCATABLE ::                 &
                    SNAPU, SNAPW, SNAPQX, SNAPQZ

  REAL    (KIND=PP), DIMENSION (:,:), ALLOCATABLE ::                   &
                     UM, WM, QX, QZ, QXP, QZP                             

END MODULE WAVEFIELD
