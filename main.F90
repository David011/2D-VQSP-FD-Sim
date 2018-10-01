!=======================================================================

PROGRAM FD2D_VS_PORO

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                          H
  USE INTERFACES  , ONLY:                                              &
                          ALLOC_NR, ALLOC_T_SEIS, AUX1,                &
                          OPEN_SR_FILE,                                &
                          PPML,        PRE_NR_COEF1,                   &
                          READ_MODEL, READ_MODEL_PAR, READ_OPEN,       &
                          TIME_LOOP
  USE TIME_FD     , ONLY:                                              &
                          time_begin, time_end, DATE_TIME, BIG_BEN

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=20) :: JOB_NAME

!-----------------------------------------------------------------------

#if(USE_MPI)
  CALL MPI_INIT (errcode)
  comm = MPI_COMM_WORLD
  CALL MPI_COMM_RANK (comm, myid,     errcode)
  CALL MPI_Get_ProCessor_name   (nam_proc, len, errcode)
  CALL MPI_COMM_SIZE (comm, numprocs, errcode)
  WRITE(0,*) "Process ", myid, "of", numprocs, "is alive : ",nam_proc
  CALL MPI_Barrier   (comm,           errcode)
  root = 0
#endif

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    CALL CPU_TIME ( time_begin )
    CALL DATE_AND_TIME ( BIG_BEN(1), BIG_BEN(2), BIG_BEN(3), DATE_TIME )
    WRITE (*,'(''Current time: '',I2,'':'',I2,'':'',F6.3)')            &
      DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)+REAL(DATE_TIME(8)/1000.)

    OPEN  (10, FILE = 'HF_2DFD_VS_PORO', STATUS='OLD')
    READ  (10,'(A)') JOB_NAME
    CLOSE (10)
    
#if(USE_MPI)
  END IF
#endif

  CALL  READ_OPEN ( JOB_NAME (1:LEN_TRIM(JOB_NAME)) )
                         print *,  'read_open'
  CALL  DISTRIBUTE
                         print *,  'distribute'
  CALL  AUX1
                         print *,  'aux1'
  CALL  DOMAIN_DECOMPOSITION
                         print *,  'domain_decomposition'
  CALL  ALLOC_NR
                         print *,  'alloc_nr'
  CALL  ALLOC_T_SEIS
                         print *,  'alloc_t_seis'
  CALL  READ_MODEL
                         print *,  'read_model'
                         
  CALL  READ_MODEL_PAR
  
                         print *,  'read_model_par'
  CALL  OPEN_SR_FILE
                         print *,  'open_sr_file'
  CALL  PRE_NR_COEF1
                         print *,  'pre_nr_coef1'
  CALL  PPML
                         print *,  'ppml'
  CALL  TIME_LOOP
                         print *,  'time_loop'

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    CALL CPU_TIME ( time_end )
    CALL DATE_AND_TIME ( BIG_BEN(1), BIG_BEN(2), BIG_BEN(3), DATE_TIME )
    PRINT     *, 'CPU time: ', time_end - time_begin, ' seconds'
    WRITE (*,'(''Current time: '',I2,'':'',I2,'':'',F6.3)')            &
      DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)+REAL(DATE_TIME(8)/1000.)
    WRITE (11,'(F10.2,1X,I4,I3,I3,I5,3I3,I4)')                         &
             time_end - time_begin, DATE_TIME
#if(USE_MPI)
  END IF
#endif

END PROGRAM FD2D_VS_PORO
