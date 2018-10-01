!=======================================================================

SUBROUTINE  STORE_SNAPSHOTS(ITIST)

  USE PRECISION   , ONLY: PP

  USE CONTROL_DATA, ONLY: MZ, STRESS_IMAGING, KEY_TLQ

  USE AUXIL       , ONLY: X_MIN, X_MAX,               Z_MIN, Z_MAX,    &
                          BUFS, BUFQ

  USE WAVEFIELD   , ONLY: UM, WM, QX, QZ, SNAPU, SNAPW, SNAPQX, SNAPQZ

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ITIST

  CHARACTER (LEN = 14) :: TLBUF

  INTEGER :: J, JJ, I, L, SIZE

!-----------------------------------------------------------------------

#if(USE_MPI)
  IF ( myid == root ) THEN
      
     DO J = 0, numprocs-1
         
        CALL MPI_CART_COORDS (comm_cart, J, ndims, coord, errcode)

        SIZE = ( X_MAX_A(coord(1)) - X_MIN_A(coord(1)) + 1 )             &
           * ( Z_MAX_A(coord(2)) - Z_MIN_A(coord(2)) + 1 )

      IF ( J /= root ) THEN

        tag  = 9

        CALL MPI_RECV (BUFS , SIZE, MPI_REAL, J, tag, comm_cart,       &
                                                 status, errcode)
        JJ=1
        DO L = Z_MIN_A (coord(2)), Z_MAX_A (coord(2))
          DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
            SNAPU  ( I, L )=BUFS (JJ)
            JJ = JJ + 1
          END DO
        END DO

        CALL MPI_RECV (BUFS , SIZE, MPI_REAL, J, tag, comm_cart,       &
                                                 status, errcode)
        JJ=1
        DO L = Z_MIN_A (coord(2)), Z_MAX_A (coord(2))
          DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
            SNAPW  ( I, L )=BUFS (JJ)
            JJ = JJ + 1
          END DO
        END DO
        
        IF(KEY_TLQ) THEN
           tagq  = 10
        
           CALL MPI_RECV (BUFQ , SIZE, MPI_REAL, J, tagq, comm_cart,       &
                                                 status, errcode)
           JJ=1
           DO L = Z_MIN_A (coord(2)), Z_MAX_A (coord(2))
              DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
                 SNAPQX  ( I, L )=BUFQ (JJ)
                 JJ = JJ + 1
              END DO
           END DO

           CALL MPI_RECV (BUFQ , SIZE, MPI_REAL, J, tagq, comm_cart,       &
                                                 status, errcode)
           JJ=1
           DO L = Z_MIN_A (coord(2)), Z_MAX_A (coord(2))
              DO I = X_MIN_A (coord(1)), X_MAX_A (coord(1))
                 SNAPQZ  ( I, L )=BUFQ (JJ)
                 JJ = JJ + 1
              END DO
           END DO
        END IF

      ELSE
        SNAPU( X_MIN : X_MAX, Z_MIN : Z_MAX ) =                        &
          UM ( X_MIN : X_MAX, Z_MIN : Z_MAX    )
        SNAPW( X_MIN : X_MAX, Z_MIN : Z_MAX ) =                        &
          WM ( X_MIN : X_MAX, Z_MIN : Z_MAX    )
        
        IF(KEY_TLQ) THEN
           SNAPQX( X_MIN : X_MAX, Z_MIN : Z_MAX ) =                        &
               QX ( X_MIN : X_MAX, Z_MIN : Z_MAX    )
           SNAPQZ( X_MIN : X_MAX, Z_MIN : Z_MAX ) =                        &
               QZ ( X_MIN : X_MAX, Z_MIN : Z_MAX    )
        END IF
      END IF

    END DO


    WRITE ( TLBUF,'(A6, I6.6)' )  'SNAP.V', ITIST
    OPEN  ( 18, FILE   = TLBUF(1:12),  FORM   ='UNFORMATTED',          &
            STATUS = 'REPLACE')

    WRITE ( 18 ) SNAPU, SNAPW
    CLOSE ( 18 )
    
    IF(KEY_TLQ) THEN
       WRITE ( TLBUF,'(A6, I6.6)' )  'SNAP.Q', ITIST
       OPEN  ( 19, FILE   = TLBUF(1:12),  FORM   ='UNFORMATTED',          &
                   STATUS = 'REPLACE')

       WRITE ( 19 ) SNAPQX, SNAPQZ
       CLOSE ( 19 )
    END IF
    
  ELSE

    SIZE = ( X_MAX  - X_MIN  + 1 ) * ( Z_MAX  - Z_MIN  + 1 )

    tag=9

    JJ=1
    DO L = Z_MIN , Z_MAX
      DO I = X_MIN , X_MAX
        BUFS (JJ) = UM  ( I, L )
        JJ = JJ + 1
      END DO
    END DO
    CALL MPI_SEND(BUFS ,SIZE,MPI_REAL,root,tag,comm_cart,errcode)

    JJ=1
    DO L = Z_MIN , Z_MAX
      DO I = X_MIN , X_MAX
        BUFS (JJ) = WM  ( I, L )
        JJ = JJ + 1
      END DO
    END DO
    CALL MPI_SEND(BUFS ,SIZE,MPI_REAL,root,tag,comm_cart,errcode)
    
    IF(KEY_TLQ) THEN
        
       tagq=10

       JJ=1
       DO L = Z_MIN , Z_MAX
          DO I = X_MIN , X_MAX
             BUFQ (JJ) = QX  ( I, L )
             JJ = JJ + 1
          END DO
       END DO
       CALL MPI_SEND(BUFQ ,SIZE,MPI_REAL,root,tagq,comm_cart,errcode)

       JJ=1
       DO L = Z_MIN , Z_MAX
          DO I = X_MIN , X_MAX
             BUFQ (JJ) = QZ  ( I, L )
             JJ = JJ + 1
          END DO
       END DO
       CALL MPI_SEND(BUFQ ,SIZE,MPI_REAL,root,tagq,comm_cart,errcode)
    END IF
    
  END IF

#else
  WRITE ( TLBUF,'(A6, I6.6)' )  'SNAP.V', ITIST

  OPEN  ( 18, FILE   = TLBUF(1:12),  FORM   ='UNFORMATTED',            &
              STATUS = 'REPLACE')

  WRITE ( 18 ) UM (: ,: ), WM (: ,: )

  CLOSE ( 18 )
  
  IF(KEY_TLQ) THEN
     WRITE ( TLBUF,'(A6, I6.6)' )  'SNAP.Q', ITIST

     OPEN  ( 19, FILE   = TLBUF(1:12),  FORM   ='UNFORMATTED',            &
                 STATUS = 'REPLACE')

     WRITE ( 19 ) QX (: ,: ), QZ (: ,: )

     CLOSE ( 19 )
  
  END IF
#endif

END SUBROUTINE  STORE_SNAPSHOTS
