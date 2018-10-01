!_______________________________________________________________________

MODULE GRID_MEDIUM
  USE PRECISION, ONLY: PP, IK
  IMPLICIT NONE

  CHARACTER(LEN=40) :: TEXT

  INTEGER           :: JMNUM                                            

  INTEGER(IK), DIMENSION (:,:), ALLOCATABLE :: JM                       
!$SGI DISTRIBUTE JM(*,BLOCK,*)

  REAL (KIND=PP) :: H

  REAL (KIND=PP), DIMENSION (:), ALLOCATABLE ::                        &
                           MU,                                         &
                           SIG_XX_1,       SIG_XX_2,       SIG_XX_3,   &
                           SIG_ZZ_1,       SIG_ZZ_2,       SIG_ZZ_3,   &
                           PRES_1,         PRES_2,         PRES_3,     &
                           REL1_U,         REL2_U,         REL3_U,     &
                           REL1_W,         REL2_W,         REL3_W,     &
                           REL1_QX,        REL2_QX,        REL3_QX,    &
                           REL1_QZ,        REL2_QZ,        REL3_QZ,    &
                           AUXIL_W,        AUXIL_QZ

END MODULE GRID_MEDIUM
