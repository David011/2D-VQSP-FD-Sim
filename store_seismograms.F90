SUBROUTINE  STORE_SEISMOGRAMS (T)
                                                    



  USE PRECISION   , ONLY: PP

  USE GRID_MEDIUM , ONLY: H

  USE CONTROL_DATA, ONLY: MR, DT, REC_NAME, KEY_TLQ, IRECX, LRECX, IRECZ, LRECZ

  USE AUXIL       , ONLY: X_MIN ,X_MAX , Z_MIN ,Z_MAX , ORDER

  USE WAVEFIELD   , ONLY: UM     , WM     , QX      , QZ      ,                  &
                          SEISU  , SEISW  , DSEISU  , DSEISW  , SSU  , SSW ,     &      
                          SEISQX , SEISQZ , DSEISQX , DSEISQZ , SSQX , SSQZ,     &    
                          SSUF, SSWF  
  
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

  REAL(PP) :: SP1, SM1, S, ZP1, ZM1, Z, INTERPOL

  REAL(PP), DIMENSION (1:MR) :: SU, SW, SPU, SPW, ZU, ZW, ZPU, ZPW
  
  CHARACTER(LEN=20) :: STR, FMT


  SU = 0.
  SW = 0.
  SPU = 0.
  SPW = 0.

  DO J = 1, MR

      IR = IRECX(J)
      LR = LRECX(J)
      
      IF ( ( X_MIN <= IR ) .AND. ( IR <= X_MAX) .AND.                  &
           ( Z_MIN <= LR ) .AND. ( LR <= Z_MAX)       ) THEN
        S     = UM (IR, LR)                                                 
                                                                           
        SP1   = UM (IR, LR)                                                
 
        SM1   = UM (IR, LR)                                              
        IF (KEY_TLQ) THEN
            Z     = QX (IR, LR)
            ZP1   = QX (IR, LR)
            ZM1   = QX (IR, LR)
        END IF
        
        SPU(J) = (S - SSU(J))/DT
        SSU(J) =  S
        SU (J) = (SP1-SM1)/H/2.
        
        IF (KEY_TLQ) THEN
            ZPU(J)  = (Z - SSQX(J))/DT
            SSQX(J) =  Z
            ZU (J)  = (ZP1-ZM1)/H/2.
        END IF
        
      END IF
                                                                                
                                                                            
      IR = IRECZ(J)        
      LR = LRECZ(J)
      
      IF ( ( X_MIN <= IR ) .AND. ( IR <= X_MAX) .AND.                  &    
           ( Z_MIN <= LR ) .AND. ( LR <= Z_MAX)       ) THEN
        S     = -WM (IR, LR)
        SP1   = -WM (IR, LR)
        SM1   = -WM (IR, LR)
        IF (KEY_TLQ) THEN
            Z     = -QZ (IR, LR)
            ZP1   = -QZ (IR, LR)
            ZM1   = -QZ (IR, LR)
        END IF
        
        SPW(J) = (S - SSW(J))/DT
        SSW(J) =  S
        SW (J) = (SP1-SM1)/H/2.
        
        IF (KEY_TLQ) THEN
            ZPW(J) = (Z - SSQZ(J))/DT
            SSQZ(J)=  Z
            ZW (J) = (ZP1-ZM1)/H/2.
        END IF
        
      END IF

  END DO

#if(USE_MPI)
    CALL MPI_REDUCE(SU, DSEISU, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(SSU, SEISU, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(SW ,DSEISW, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(SSW, SEISW, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    
IF (KEY_TLQ) THEN   
    CALL MPI_REDUCE(ZU, DSEISQX, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(SSQX, SEISQX, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(ZW ,DSEISQZ, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
    CALL MPI_REDUCE(SSQZ, SEISQZ, MR, MPI_REAL, MPI_SUM, root, comm_cart, errcode)
END IF
    
    IF ( myid == root ) THEN
         WRITE ( STRING, '(I1)' ) MR
         FMT = '('//TRIM(STRING)//'(1X,E13.6))'
          DO J = 1, MR
             WRITE ( TLBUF, '(I4.4)' ) J
             OPEN ( 21, FILE='REC_UW_'//TRIM(STR(J))//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
             WRITE (21, '(5(1X,E13.6))')  DT*(T), SEISU (J), SEISW (J), DSEISU (J),DSEISW (J)
      
             CLOSE (21)
             IF (KEY_TLQ) THEN   
                 OPEN ( 22, FILE='REC_Q_'//TRIM(STR(J))//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
                 WRITE (22, '(5(1X,E13.6))')  DT*(T), SEISQX(J), SEISQZ(J), DSEISQX(J), DSEISQZ(J)
                 CLOSE (22)
             END IF      
          END DO
    END IF
#else
  WRITE ( STRING, '(I1)' ) MR
  FMT = '('//TRIM(STRING)//'(1X,E13.6))'
  DO J = 1, MR
     WRITE ( TLBUF, '(I4.4)' ) J
     OPEN ( 21, FILE='REC_UW_'//TRIM(STR(J))//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
     WRITE (21, '(5(1X,E13.6))')  DT*(T), SSU(J), SSW(J), SPU (J), SPW (J)
      
     CLOSE (21)
     IF (KEY_TLQ) THEN   
         OPEN ( 22, FILE='REC_Q_'//TRIM(STR(J))//'.asc', STATUS='UNKNOWN',           &
                                                   POSITION='APPEND')
         WRITE (22, '(5(1X,E13.6))')  DT*(T), SSQX(J), SSQZ(J), ZPU (J), ZPW (J)
      
         CLOSE (22)
     END IF        
  END DO
#endif

END SUBROUTINE  STORE_SEISMOGRAMS
