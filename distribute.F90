!=======================================================================

SUBROUTINE  DISTRIBUTE

  USE GRID_MEDIUM , ONLY:                                              &
                          H

  USE CONTROL_DATA, ONLY:                                              &
                          KEY_TLV, KEY_SNV,                            &
                          MT1, MT2,        IPAS2,                      &
                          DT,       MR, MZ, TPML,                      &
                          IRECX, LRECX, IRECZ, LRECZ, REC_NAME,        &
                          TAU, R, POW, SF, WC,                         &
                          KTLE, KTRI,             KTBO,                &
                          OMG, WB,                                     &
                    THPPLE, THPPRI,                 THPPBO,            &
                    THSSLE, THSSRI,                 THSSBO

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED

!-----------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=6) :: REC_NAME_T
  INTEGER          :: J

!-----------------------------------------------------------------------

  CALL MPI_Bcast(H,         1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(KEY_SNV  , 1, MPI_LOGICAL, root, comm, errcode)
  CALL MPI_Bcast(MT1,       1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MT2,       1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(IPAS2,     1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(DT,        1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(TAU,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(R,         1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(POW,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(SF,        1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(WC,        1, MPI_REAL   , root, comm, errcode)

  CALL MPI_Bcast(KTLE,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KTRI,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KTBO,      1, MPI_INTEGER, root, comm, errcode)

  CALL MPI_Bcast(OMG,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(WB ,       1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPLE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPRI,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THPPBO,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSLE,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSRI,    1, MPI_REAL   , root, comm, errcode)
  CALL MPI_Bcast(THSSBO,    1, MPI_REAL   , root, comm, errcode)

  IF (KEY_TLV) THEN
      
    DO J = 1, MR
      REC_NAME_T = REC_NAME (J)
      CALL MPI_Bcast(REC_NAME_T, 6, MPI_CHARACTER, root, comm, errcode)                                                                  
      REC_NAME(J) = REC_NAME_T                                         
    END DO
    CALL MPI_Bcast(IRECX,   MR, MPI_INTEGER, root, comm, errcode)
    CALL MPI_Bcast(LRECX,   MR, MPI_INTEGER, root, comm, errcode)
    CALL MPI_Bcast(IRECZ,   MR, MPI_INTEGER, root, comm, errcode)
    CALL MPI_Bcast(LRECZ,   MR, MPI_INTEGER, root, comm, errcode)
  END IF

#endif

END SUBROUTINE  DISTRIBUTE
