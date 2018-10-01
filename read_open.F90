!=======================================================================

SUBROUTINE  READ_OPEN  ( WORK_NAME )

  USE PRECISION   , ONLY:                                              &
                          PP

  USE GRID_MEDIUM , ONLY:                                              &
                          TEXT, JM, H

  USE CONTROL_DATA, ONLY:                                              &
                          SR_FILE_NAME, REC_NAME,                      &     
                          KEY_TLV, KEY_SNV, KEY_TLQ, KEY_TLP,          &          
                          MX,     MZ,                                  &          
                          MT1, MT2,        IPAS2, MR,                  &    
                          KTLE, KTRI,      KTBO,  KEY_SOUR,            &      
                          DT,                TAU_EPS,                  &
                          OMG, WB, QR,                                 &
                          THPPLE, THPPRI,                 THPPBO,      &
                          THSSLE, THSSRI,                 THSSBO,      &
                          IRECX, LRECX, IRECZ, LRECZ,                  &      
                          TPML, R, TAU, POW, SF, WC,                   &
                          MO_FILE_NAME, JMH_FILE_NAME,                 &
                          STRESS_IMAGING,                              &   
                          LIN_X, LIN_Z, DT_VS

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT (IN) :: WORK_NAME

  CHARACTER (LEN=25) LOG_FILE_NAME,  IN_FILE_NAME, REC_FILE_NAME   

  LOGICAL :: EX

  INTEGER :: IOS, I, ALLOSTAT, PROCX, PROCZ                    
                                                                    
  REAL    :: XBMIN, X_SOURCE, Z_SOURCE, MT                  
                                                                    
  NAMELIST /NAMES/        MO_FILE_NAME,  JMH_FILE_NAME

  NAMELIST /NAMES_SR_REC/ SR_FILE_NAME, REC_FILE_NAME

  NAMELIST /KEYS/         KEY_TLV, KEY_SNV, KEY_TLQ, KEY_TLP

  NAMELIST /CONTROLDATA1/ MT, R, TAU , POW, SF, WC, TAU_EPS, IPAS2

  NAMELIST /CONTROLDATA2/ MX, MZ,  H, DT_VS, XBMIN

  NAMELIST /CONTROLDATA/  TPML, PROCX, PROCZ, STRESS_IMAGING

  NAMELIST /SOURCE/       X_SOURCE, Z_SOURCE, KEY_SOUR, QR

  NAMELIST /NONREF/       OMG, WB,                                     &
                            KTLE,   KTRI,                   KTBO,      &
                          THPPLE, THPPRI,                 THPPBO,      &
                          THSSLE, THSSRI,                 THSSBO

  NAMELIST /TXT/          TEXT

  NAMELIST /REC/          MR
  
#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
!______________________________________________ LOG FILE __________ 11 _

    LOG_FILE_NAME = WORK_NAME//'.LOG'                                    
                                                                     
    OPEN  (11, FILE = LOG_FILE_NAME, STATUS='UNKNOWN' )
    WRITE (11,*) 'LOG_FILE: JOB_NAME =', WORK_NAME

!______________________________________________ CONTROL DATA FILE _ 12 _

    IN_FILE_NAME = WORK_NAME//'.IN'                                        
                                                                           
    INQUIRE ( FILE = IN_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (IN_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) IN_FILE_NAME,' NOT FOUND'
        STOP
      END IF
    OPEN  (12, FILE = IN_FILE_NAME, STATUS='OLD' )

!______________________________________________ READ CONTROL DATA ______

    OPEN ( 22, FILE = WORK_NAME//'.INM', STATUS = 'OLD' )        
                                                           
     MO_FILE_NAME = WORK_NAME//'.MO'                          
    JMH_FILE_NAME = WORK_NAME//'.JMH '
     SR_FILE_NAME = 'STF.DAT'                 
    REC_FILE_NAME = 'REC.DAT'               

    TPML         = 0
    PROCX        = 0
    PROCZ        = 0
    STRESS_IMAGING = .TRUE.

    READ  (12, NML=CONTROLDATA)
    WRITE (11, NML=CONTROLDATA)

    READ  (22, NML=NAMES)
    WRITE (11, NML=NAMES)

    QR           = 0
    R            = 0.001
    TAU          = 1.5
    POW          = 2
    SF           = 0.
    WC           = 0.
    MT1          = 1
    XBMIN        = 0.
    IPAS2        = 1000
    READ  (12, NML=CONTROLDATA1)
    WRITE (11, NML=CONTROLDATA1)

    READ  (12, NML=NAMES_SR_REC)
    WRITE (11, NML=NAMES_SR_REC)

    KEY_TLV = .TRUE.
    KEY_SNV = .FALSE.
    KEY_TLP = .FALSE.
    READ  (12, NML=KEYS)
    WRITE (11, NML=KEYS)

    IF  ( .NOT. ( KEY_TLV .OR. KEY_SNV ) ) THEN
      WRITE (11,*) ' NONE OF TL_ AND SN_ OUTPUT FILES IS REQUIRED'
      STOP
    END IF

    READ ( 22, NML = CONTROLDATA2 )
    WRITE( 11, NML = CONTROLDATA2 )
    
    DT = DT_VS

    IF (STRESS_IMAGING) THEN
      MZ = MZ+2                             
    END IF                               

    MT2 = INT(MT/DT)+1
    KEY_SOUR = 0
    
#if(USE_MPI)
    dims(1) = PROCX
    dims(2) = PROCZ
#endif

    READ  (12, NML=SOURCE)
    WRITE (11, NML=SOURCE)
    X_SOURCE = ABS(X_SOURCE)
    Z_SOURCE = ABS(Z_SOURCE)
    IF (STRESS_IMAGING) THEN
      Z_SOURCE = Z_SOURCE + 2.*H         
    END IF
    LIN_X = NINT((X_SOURCE-XBMIN)/H)
    LIN_Z = NINT(Z_SOURCE/H)
    
    TAU_EPS = 0.                          
    CLOSE (22)

    WB     = 0.4_PP                 
    THPPLE = 0.19757_PP             
    THPPRI = 0.19757_PP
    THPPBO = 0.19757_PP
    THSSLE = 0.70511_PP
    THSSRI = 0.70511_PP
    THSSBO = 0.70511_PP
    
    KTLE = 0
    KTRI = 0
    KTBO = 0
    
    READ  (12, NML=NONREF)
    WRITE (11, NML=NONREF)
    OMG = OMG/6.283185307179586476925286766559  

    READ  (12, NML=TXT)
    WRITE (11, NML=TXT)
#if(USE_MPI)
  END IF

  CALL MPI_Bcast( JMH_FILE_NAME, 20, MPI_CHARACTER, root, comm, errcode)
  CALL MPI_Bcast(  MO_FILE_NAME, 20, MPI_CHARACTER, root, comm, errcode)

  CALL MPI_Bcast(TPML,      1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MX,        1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(MZ,        1, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(KEY_TLV  , 1, MPI_LOGICAL, root, comm, errcode)
  CALL MPI_Bcast(dims,  ndims, MPI_INTEGER, root, comm, errcode)
  CALL MPI_Bcast(STRESS_IMAGING, 1, MPI_LOGICAL, root, comm, errcode)

  CALL MPI_Bcast(LIN_X      , 1, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast(LIN_Z      , 1, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast(KEY_SOUR   , 1, MPI_INTEGER, root, comm_cart, errcode)
  CALL MPI_Bcast(KEY_TLQ    , 1, MPI_LOGICAL, root, comm_cart, errcode)
  CALL MPI_Bcast(KEY_TLP    , 1, MPI_LOGICAL, root, comm_cart, errcode)
  CALL MPI_Bcast(QR         , 1, MPI_DOUBLE_PRECISION, root, comm_cart, errcode)
#endif

  IF  ( KEY_TLV )  THEN

#if(USE_MPI)
    IF ( myid == root ) THEN
#endif
      READ  (12, NML=REC)
      WRITE (11, NML=REC)
#if(USE_MPI)
    END IF
    
    CALL MPI_Bcast(MR, 1, MPI_INTEGER, root, comm, errcode)
#endif
    ALLOCATE ( IRECX (1:MR), LRECX (1:MR), IRECZ (1:MR), LRECZ (1:MR), REC_NAME(1:MR), STAT=ALLOSTAT )
      IF  ( ALLOSTAT > 0 )  THEN
        WRITE (11,*) '  ALLOCATION ERROR R_O_02, STAT=', ALLOSTAT
        STOP
      END IF
      
#if(USE_MPI)
    IF ( myid == root ) THEN
#endif
      OPEN ( 22, FILE = REC_FILE_NAME, STATUS = 'OLD' )
      READ (22,*)  ( REC_NAME(I), IRECX(I), LRECX(I), IRECZ(I), LRECZ(I), I=1,MR )
      CLOSE( 22 )
#if(USE_MPI)
    END IF
#endif
  END IF

#if(USE_MPI)
  IF ( myid == root ) THEN
#endif
    CLOSE (12)

#if(USE_MPI)
#else

!___________________________________________ MODEL FILES _ 13,14,15 - vysledok model preparation

    INQUIRE ( FILE = MO_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (MO_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) MO_FILE_NAME,' NOT FOUND'
        STOP
      END IF

    INQUIRE ( FILE = JMH_FILE_NAME, EXIST=EX, IOSTAT=IOS )
      IF  ( IOS > 0 )  THEN
        WRITE (11,*) 'ERROR HAS OCCURED DURING THE EXECUTION           &
                     &OF THE INQUIRE STATEMENT.',                      &
                     'IOSTAT (JMH_FILE_NAME)=',IOS
        STOP
      END IF
      IF  ( .NOT.EX )  THEN
        WRITE (11,*) JMH_FILE_NAME,' NOT FOUND'
        STOP
      END IF

    OPEN ( 13, FILE   = JMH_FILE_NAME, FORM = 'UNFORMATTED',           &   
                                                          STATUS='OLD' )

    OPEN ( 14, FILE = MO_FILE_NAME, FORM = 'UNFORMATTED', STATUS='OLD' )
#endif

#if(USE_MPI)
  END IF
#endif

END SUBROUTINE  READ_OPEN
