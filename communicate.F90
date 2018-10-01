SUBROUTINE COMMUNICATE (VARIABLE)

  USE PRECISION   , ONLY : PP
  USE AUXIL       , ONLY : X_MIN_D, X_MAX_D, X_MIN, X_MAX,             &
                           Z_MIN_D, Z_MAX_D, Z_MIN, Z_MAX,             &
                           BUFXS,        BUFZS,                        &
                           BUFXR,        BUFZR,                        &
                           source_id0A, dest_id0A,                     &
                           source_id0B, dest_id0B,                     &
                           source_id1A, dest_id1A,                     &
                           source_id1B, dest_id1B

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  REAL(PP), INTENT(INOUT),                                             &
    DIMENSION(X_MIN_D:X_MAX_D, Z_MIN_D:Z_MAX_D) ::  VARIABLE


  INTEGER       :: stag, rtag, SIZE

!-----------------------------------------------------------------------

#if(USE_MPI)
!!  CALL MPI_BARRIER (comm_cart,errcode)

! --- communicate in x direction

  SIZE = OVERLAP*(Z_MAX-Z_MIN+1)   !OVERLAP=2

  IF ( SIZE > 0 )  THEN

    stag = 0
    rtag = 0                                                                
    IF ( dest_id0A /= MPI_PROC_NULL ) THEN                                 
      BUFXS(1:OVERLAP,1:Z_MAX-Z_MIN+1) =                               &
        VARIABLE (X_MAX+1-OVERLAP:X_MAX,Z_MIN:Z_MAX)
      IF ( source_id0A == MPI_PROC_NULL ) THEN                              
        CALL MPI_SEND            (BUFXS, SIZE, MPI_REAL,               &    
                                   dest_id0A, stag,                    &    
                                   comm_cart, errcode)                      
                                                                            
                                                                            
                                                                            
        
      ELSE                                                                  
        CALL MPI_ISEND           (BUFXS, SIZE, MPI_REAL,               &
                                   dest_id0A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &
                                                    source_id0A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_D:X_MIN_D+OVERLAP-1,            Z_MIN:Z_MAX) = &
          BUFXR(1:OVERLAP,                1:Z_MAX-Z_MIN+1)
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0A /= MPI_PROC_NULL ) THEN                              
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &    
                                                    source_id0A, rtag, &    
                                   comm_cart, status, errcode)
        VARIABLE (X_MIN_D:X_MIN_D+OVERLAP-1,            Z_MIN:Z_MAX) = &
          BUFXR(1:OVERLAP,                1:Z_MAX-Z_MIN+1)
      END IF
    END IF

    stag = 1
    rtag = 1
    IF ( dest_id0B /= MPI_PROC_NULL ) THEN                                
      BUFXS = VARIABLE (X_MIN:X_MIN+OVERLAP-1,            Z_MIN:Z_MAX)    
      IF ( source_id0B == MPI_PROC_NULL ) THEN                            
        CALL MPI_SEND            (BUFXS, SIZE, MPI_REAL,               &  
                                   dest_id0B, stag,                    &  
                                   comm_cart, errcode)                    
        
        
      ELSE                                                                
        CALL MPI_ISEND           (BUFXS, SIZE, MPI_REAL,               &
                                   dest_id0B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &
                                                    source_id0B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_D+1-OVERLAP:X_MAX_D,            Z_MIN:Z_MAX) = &  
          BUFXR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id0B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFXR, SIZE, MPI_REAL,               &   
                                                    source_id0B, rtag, &   
                                   comm_cart, status, errcode)
        VARIABLE (X_MAX_D+1-OVERLAP:X_MAX_D,            Z_MIN:Z_MAX) = &
          BUFXR
      END IF
    END IF

  END IF


! --- communicate in z direction

  SIZE = (X_MAX_D-X_MIN_D+1)*OVERLAP

  IF ( Z_MIN <= Z_MAX ) THEN

    stag = 4
    rtag = 4
    IF ( dest_id1A /= MPI_PROC_NULL ) THEN
      BUFZS = VARIABLE                                                 &
                (X_MIN_D:X_MAX_D,                Z_MAX+1-OVERLAP:Z_MAX)
      IF ( source_id1A == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id1A, stag,                    &
                                   comm_cart, errcode)
      ELSE
        CALL MPI_ISEND           (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id1A, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
         (X_MIN_D:X_MAX_D,                Z_MIN_D:Z_MIN_D+OVERLAP-1) = &
         BUFZR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1A /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id1A, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
         (X_MIN_D:X_MAX_D,                Z_MIN_D:Z_MIN_D+OVERLAP-1) = &
         BUFZR
      END IF
    END IF


    stag = 5
    rtag = 5
    IF ( dest_id1B /= MPI_PROC_NULL ) THEN
      BUFZS = VARIABLE                                                 &
                (X_MIN_D:X_MAX_D,                Z_MIN:Z_MIN+OVERLAP-1)
      IF ( source_id1B == MPI_PROC_NULL ) THEN
        CALL MPI_SEND            (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id1B, stag,                    &
                                   comm_cart, errcode)
      ELSE
        CALL MPI_ISEND           (BUFZS, SIZE, MPI_REAL,               &
                                   dest_id1B, stag,                    &
                                   comm_cart, request, errcode)
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
           (X_MIN_D:X_MAX_D,                Z_MAX_D+1-OVERLAP:Z_MAX_D) &
           = BUFZR
        CALL MPI_WAIT (request, status, errcode)
      END IF
    ELSE
      IF ( source_id1B /= MPI_PROC_NULL ) THEN
        CALL MPI_RECV            (BUFZR, SIZE, MPI_REAL,               &
                                                    source_id1B, rtag, &
                                   comm_cart, status, errcode)
        VARIABLE                                                       &
           (X_MIN_D:X_MAX_D,                Z_MAX_D+1-OVERLAP:Z_MAX_D) &
           = BUFZR
      END IF
    END IF

  END IF
#endif


END SUBROUTINE COMMUNICATE
