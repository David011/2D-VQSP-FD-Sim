SUBROUTINE  STORE_SEISMOGRAMS_PRESSURE (T)
    
    
    
    
    
  USE PRECISION   , ONLY: PP

  USE GRID_MEDIUM , ONLY: H

  USE CONTROL_DATA, ONLY: MR, DT, REC_NAME, IRECZ, LRECX
  
  USE AUXIL       , ONLY: X_MIN ,X_MAX , Z_MIN ,Z_MAX , ORDER
  
  USE STRESSFIELD , ONLY: PRES, TXX

  USE WAVEFIELD   , ONLY: SEISP, SEIST

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif 
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: T

  CHARACTER (LEN = 11) :: TLBUF, STRING

  INTEGER :: J, IR, LR

  REAL(PP), DIMENSION (1:MR) ::  SSP, SST
  
  CHARACTER(LEN=20) :: STR, FMT

!-----------------------------------------------------------------------

SSP = 0.
SST = 0.


DO J = 1, MR
    
    IR = IRECZ(J)
    LR = LRECX(J)

    IF ( ( X_MIN <= IR ) .AND. ( IR <= X_MAX) .AND.                  &
           ( Z_MIN <= LR ) .AND. ( LR <= Z_MAX)       ) THEN
        
           SSP(J) =   PRES(IR, LR)
           SST(J) = - TXX (IR, LR)
        
    END IF
   
END DO
#if(USE_MPI)
    CALL MPI_REDUCE(SSP, SEISP, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(SST, SEIST, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)

    IF ( myid == root ) THEN
         WRITE ( STRING, '(I1)' ) MR
         FMT = '('//TRIM(STRING)//'(1X,E13.6))'
          DO J = 1, MR
             WRITE ( TLBUF, '(I4.4)' ) J
             OPEN ( 27, FILE='REC_STRESS_'//TRIM(STR(J))//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
             WRITE (27, '(3(1X,E13.6))')  DT*(T), SEISP (J), SEIST (J)
      
             CLOSE (27)      
          END DO
    END IF
#else  
    DO J = 1, MR
       WRITE ( TLBUF, '(I4.4)' ) J
       OPEN ( 27, FILE='REC_STRESS_'//TRIM(STR(J))//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
       WRITE (27, '(3(1X,E13.6))')  DT*(T), SSP(J), SST(J)
      
       CLOSE (27)   
     END DO
#endif

END SUBROUTINE STORE_SEISMOGRAMS_PRESSURE