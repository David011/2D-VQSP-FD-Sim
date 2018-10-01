!=======================================================================

SUBROUTINE ALLOC_T_SEIS

  USE PRECISION   , ONLY:                                                &
                          PP

  USE CONTROL_DATA, ONLY:                                                       &
                          MR, MZ, SOURSP, SOURFP, KEY_TLQ

  USE AUXIL       , ONLY:                                                       &
                          PMLL , PMLXH , PMLZH,                                 & 
                          X_MIN_D , X_MAX_D                                       

  USE STRAINFIELD , ONLY:                                                       &
                          EXX , EXZ , EZZ, QXX, QZZ

  USE WAVEFIELD   , ONLY:                                                       &
                          SEISU ,  SEISW , DSEISU  , DSEISW  , SSU  , SSW  ,    &
                          SEISQX , SEISQZ, DSEISQX , DSEISQZ , SSQX , SSQZ ,    &
                          SEISP , SEIST,                                        &
                          SSUF  , SSWF  ,                                       &
                          SNAPU , SNAPW, SNAPQX, SNAPQZ

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: ALLOSTAT

!-----------------------------------------------------------------------

  ALLOCATE ( EXX (X_MIN_D:X_MAX_D                 ),                            &
             EZZ (X_MIN_D:X_MAX_D                 ),                            &
             EXZ (X_MIN_D:X_MAX_D                 ),                            &
             QXX (X_MIN_D:X_MAX_D                 ),                            &
             QZZ (X_MIN_D:X_MAX_D                 ),                            &
                                                        STAT=ALLOSTAT )
    IF  ( ALLOSTAT > 0 )  THEN
      WRITE (11,*) '  ALLOCATION ERROR ALLOC_T_SEIS_2, STAT= ', ALLOSTAT
      STOP
    END IF

  ALLOCATE (  SEISU  (1:MR), SEISW  (1:MR), SEISQX  (1:MR), SEISQZ  (1:MR),     &
              SEISP  (1:MR), SEIST  (1:MR),                                     &
              SSU    (1:MR), SSW    (1:MR), SSQX    (1:MR), SSQZ    (1:MR),     &
              SSUF    (1:MR), SSWF    (1:MR),                                   &
              DSEISU (1:MR), DSEISW (1:MR), DSEISQX (1:MR), DSEISQZ (1:MR),     &             
              STAT=ALLOSTAT )
      IF  ( ALLOSTAT > 0 )  THEN
        WRITE (11,*) 'ALLOCATION ERROR ALLOC_T_SEIS_5, STAT= ', ALLOSTAT
        STOP
      END IF

   SSU  = 0.
   SSW  = 0.
   SSQX = 0.
   SSQZ = 0.
   SSUF = 0.
   SSWF = 0.
   
   SOURSP = 0.
   SOURFP = 0.
   
#if(USE_MPI)
  IF ( myid == root ) THEN
    ALLOCATE ( SNAPU (PMLL:PMLXH,0: PMLZH   ),                                     &
               SNAPW (PMLL:PMLXH,0: PMLZH   )                                      &
             )
        
    IF(KEY_TLQ) THEN    
       ALLOCATE ( SNAPQX(PMLL:PMLXH,0: PMLZH   ),                                  &
                  SNAPQZ(PMLL:PMLXH,0: PMLZH   )                                   &
                )
    END IF
  END IF
#endif

END SUBROUTINE ALLOC_T_SEIS
