!_______________________________________________________________________

MODULE CONTROL_DATA
  USE PRECISION, ONLY: PP

  IMPLICIT NONE

  CHARACTER(LEN=20) :: SR_FILE_NAME,                                   &
                       MO_FILE_NAME, JMH_FILE_NAME

  CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE :: REC_NAME

  LOGICAL       ::  KEY_TLV, KEY_SNV, KEY_TLQ, KEY_TLP, STRESS_IMAGING

  INTEGER       ::  MX,     MZ,                                        &
                    MT1, MT2,        IPAS2, MR,                        &
                    TPML, LIN_X, LIN_Z                                 

  INTEGER       ::  KTLE, KTRI, KTBO, KEY_SOUR

  REAL (KIND=PP)::  DT, TAU_EPS, DT_VS, QR
                                                                          
  REAL (KIND=PP)::  R, TAU, POW, SF, WC, PHI, SOURSP, SOURFP              

  REAL (KIND=PP)::  OMG, WB,                                           &
                    THPPLE, THPPRI,                 THPPBO,            & 
                    THSSLE, THSSRI,                 THSSBO               

  INTEGER, DIMENSION (:) , ALLOCATABLE ::                       &
                    IRECX, LRECX, IRECZ, LRECZ                                            


END MODULE CONTROL_DATA
