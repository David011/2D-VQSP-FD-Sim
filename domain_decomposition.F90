!=======================================================================

SUBROUTINE  DOMAIN_DECOMPOSITION

  USE PRECISION   , ONLY:                                              & 
                           PP                                            
  USE AUXIL       , ONLY:                                              &
                           PMLL , PMLXH  ,          PMLZH ,            &
                           X_MIN_D , X_MAX_D , X_MIN , X_MAX ,         &   
                           Z_MIN_D , Z_MAX_D , Z_MIN , Z_MAX,          &
                           BUFXS, BUFZS , BUFXR , BUFZR,               &
                           BUFXS_P, BUFXR_P, BUFS, BUFQ,               &
                           source_id0A , dest_id0A ,                   &
                           source_id0B , dest_id0B ,                   &
                           source_id1A , dest_id1A ,                   &
                           source_id1B , dest_id1B

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

#if(USE_MPI)
  LOGICAL :: REORDER

  INTEGER :: PROCX, PROCZ, J, root1, NDX, NDZ
  REAL(PP):: DX, DZ
#endif

!-----------------------------------------------------------------------

#if(USE_MPI)

  OVERLAP = 2
  
  periods = .FALSE.
  REORDER = .TRUE.
  
  IF ( myid == 0 ) THEN
      iamroot = .TRUE.
  ELSE
    iamroot = .FALSE.
  END IF
  
  CALL MPI_DIMS_CREATE (numprocs,ndims,dims, errcode)
  PROCX = dims(1)
  PROCZ = dims(2)
  
  CALL MPI_CART_CREATE                                                 &
               (comm, ndims, dims, periods, REORDER, comm_cart, errcode)
  
  CALL MPI_CART_GET                                                    &
          (comm_cart, ndims, dims, periods, my_coord, errcode)
  CALL MPI_COMM_RANK (comm_cart, myid, errcode)
  
  IF ( iamroot) root1 = myid
  CALL MPI_Bcast(root1,      1, MPI_INTEGER, root, comm, errcode)
  
  root = root1
  
   ALLOCATE ( X_MIN_A (0:PROCX-1), X_MAX_A (0:PROCX-1),                 &
             Z_MIN_A (0:PROCZ-1), Z_MAX_A (0:PROCZ-1) )
   
  X_MIN_A(0        ) = PMLL
  X_MAX_A(dims(1)-1) = PMLXH
  DX = REAL(PMLXH-PMLL+1)/REAL(PROCX)
  NDX = INT(DX)
  DO J = 0, (DX-NDX)*PROCX-1
    X_MAX_A(J  ) = (J+1)*(NDX+1) + PMLL
    X_MIN_A(J+1) = (J+1)*(NDX+1) + PMLL + 1
  END DO
  DO J = (DX-NDX)*PROCX, dims(1)-2
    X_MAX_A(J  ) = (J+1)*NDX + PMLL
    X_MIN_A(J+1) = (J+1)*NDX + PMLL + 1
  END DO
  
  Z_MIN_A(0        ) = 0
  Z_MAX_A(dims(2)-1) = PMLZH
  DZ = REAL(PMLZH+1)/REAL(PROCZ)
  NDZ = INT(DZ)
  DO J = 0, (DZ-NDZ)*PROCZ-1
    Z_MAX_A(J  ) = (J+1)*(NDZ+1)
    Z_MIN_A(J+1) = (J+1)*(NDZ+1) + 1
  END DO
  DO J = (DZ-NDZ)*PROCZ, dims(2)-2
    Z_MAX_A(J  ) = (J+1)*NDZ
    Z_MIN_A(J+1) = (J+1)*NDZ + 1
  END DO
  
  X_MIN = X_MIN_A ( my_coord(1) )
  X_MAX = X_MAX_A ( my_coord(1) )
  Z_MIN = Z_MIN_A ( my_coord(2) )
  Z_MAX = Z_MAX_A ( my_coord(2) )
  
  ALLOCATE ( X_MIN_DA (0:PROCX-1), X_MAX_DA (0:PROCX-1),               &
             Z_MIN_DA (0:PROCZ-1), Z_MAX_DA (0:PROCZ-1) )

  X_MIN_DA(0        ) = PMLL
  X_MAX_DA(dims(1)-1) = PMLXH
  DO J = 0, dims(1)-2
    X_MAX_DA(J  ) = X_MAX_A(J  ) + OVERLAP
    X_MIN_DA(J+1) = X_MIN_A(J+1) - OVERLAP
  END DO

  Z_MIN_DA(0        ) = 0
  Z_MAX_DA(dims(2)-1) = PMLZH
  DO J = 0, dims(2)-2
    Z_MAX_DA(J  ) = Z_MAX_A(J  ) + OVERLAP
    Z_MIN_DA(J+1) = Z_MIN_A(J+1) - OVERLAP
  END DO
  
  X_MIN_D = X_MIN_DA ( my_coord(1) )
  X_MAX_D = X_MAX_DA ( my_coord(1) )
  Z_MIN_D = Z_MIN_DA ( my_coord(2) )
  Z_MAX_D = Z_MAX_DA ( my_coord(2) )
  
  ALLOCATE   ( BUFXS(      1:OVERLAP, 1:Z_MAX-Z_MIN+1),                &
               BUFZS(X_MIN_D:X_MAX_D, 1:OVERLAP) )
  ALLOCATE   ( BUFXR(      1:OVERLAP, 1:Z_MAX-Z_MIN+1),                &
               BUFZR(X_MIN_D:X_MAX_D, 1:OVERLAP) )

  ALLOCATE ( BUFXS_P (1:OVERLAP                ) )
  ALLOCATE ( BUFXR_P (1:OVERLAP                ) )

  ALLOCATE ( BUFS ( MAXVAL(X_MAX_DA-X_MIN_DA+1)                        &
                  * MAXVAL(Z_MAX_DA-Z_MIN_DA+1) ),                     &
             BUFQ ( MAXVAL(X_MAX_DA-X_MIN_DA+1)                        &
                  * MAXVAL(Z_MAX_DA-Z_MIN_DA+1) )                      &
  )
  
  CALL MPI_CART_SHIFT (comm_cart, 0, 1, source_id0A, dest_id0A, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 0,-1, source_id0B, dest_id0B, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 1, 1, source_id1A, dest_id1A, errcode)
  CALL MPI_CART_SHIFT (comm_cart, 1,-1, source_id1B, dest_id1B, errcode)
  
  print *,myid, 'X',X_MIN_D , X_MIN , X_MAX , X_MAX_D
  print *,myid, 'Z',Z_MIN_D , Z_MIN , Z_MAX , Z_MAX_D
#else 

  X_MIN_D = PMLL           
  X_MAX_D = PMLXH          
  Z_MIN_D = 0              
  Z_MAX_D = PMLZH          

  X_MIN = PMLL
  X_MAX = PMLXH
  Z_MIN = 0
  Z_MAX = PMLZH

#endif

END SUBROUTINE  DOMAIN_DECOMPOSITION
