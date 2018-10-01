!_______________________________________________________________________

MODULE AUXIL
  USE PRECISION, ONLY: PP, DP
  USE NONREF_BOUND, ONLY: NR_DATA_X

!----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------

  INTEGER       , PARAMETER :: ORDER = 2 ! FOR VALUE 4: STORE_SEISMOGRAMS HAS TO BE CHANGED !!!

#IF(SECOND)
  REAL (KIND=PP), PARAMETER :: PI=3.141592653589793,                   &
                               A = 0., B = 1.
#ELSE
  REAL (KIND=PP), PARAMETER :: PI=3.141592653589793,                   &
                               A = -1./24., B = 9./8.
#ENDIF

  INTEGER       ::  PMLL, PMLXH,        PMLZH, PMLZH1, PMLZH2,         &
                           X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &
                           Z_MIN_D , Z_MAX_D , Z_MIN , Z_MAX ,         &
                    MX1,   source_id0A , dest_id0A ,                   &
                           source_id0B , dest_id0B ,                   &
                           source_id1A , dest_id1A ,                   &
                           source_id1B , dest_id1B

  REAL (KIND=PP)::  AH, BH, CH

  REAL (KIND=PP), DIMENSION (:), ALLOCATABLE ::                        &
                    BUFS, BUFQ

  REAL (KIND=PP), DIMENSION (4) ::                                     &
                    OTA, OTB

  REAL (KIND=PP), DIMENSION (:), ALLOCATABLE ::   BUFXS_P, BUFXR_P
                                                           
  REAL (KIND=PP), DIMENSION (:,:), ALLOCATABLE :: BUFXS ,       BUFZS ,&      
                                                  BUFXR ,       BUFZR ,&         
                                                  SOURTF_A

END MODULE AUXIL
