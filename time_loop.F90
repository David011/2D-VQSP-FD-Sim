!=======================================================================

SUBROUTINE  TIME_LOOP

  USE CONTROL_DATA, ONLY:                                              &
                          KEY_SNV, KEY_TLV, KEY_TLP,                   &
                          MZ, MT1, MT2, DT,                            &
                                 IPAS2,        TPML,                   &
                          LIN_X, LIN_Z, STRESS_IMAGING
  USE GRID_MEDIUM , ONLY:                                              &
                          H
  USE T_LOOP      , ONLY:                                              &
                          ITILE
  USE AUXIL       , ONLY:                                              &
                          PMLZH, PMLZH1, PMLZH2,                       &
                          Z_MIN, Z_MAX, SOURTF_A
  USE WAVEFIELD   , ONLY:                                              &
                          UM, WM, QX, QZ
  USE STRESSFIELD , ONLY:                                              &
                          TXX , TZZ , TXZ, PRES
  
  USE TIME_FD     , ONLY:                                              &
                          time_begin, time_end, DATE_TIME, BIG_BEN
  USE NONREF_BOUND, ONLY:                                              &
                          SZ_UP, SZ_WP, SZ_U, SZ_W,                    &
                          SZ_QXP, SZ_QZP, SZ_QX, SZ_QZ
  USE INTERFACES  , ONLY:                                              &
                    HPLANE_0_STRESS       , HPLANE_0_VELOCITY ,        &
                    HPLANE_1_STRESS       , HPLANE_1_VELOCITY ,        &
                    HPLANE_IC_STRESS      , HPLANE_IC_VELOCITY,        &
                    HPLANE_LIN_STRESS     ,                            &
                    HPLANE_MZ1_STRESS     ,                            &
                    HPLANE_PML_STRESS     , HPLANE_PML_VELOCITY,       &
                    HPLANE_PMLZH0_STRESS  , HPLANE_PMLZH0_VELOCITY,    &
                    STORE_SEISMOGRAMS     , STORE_SEISMOGRAMS_PRESSURE,&
                    STORE_SNAPSHOTS
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL :: EX

  INTEGER :: L

!-----------------------------------------------------------------------

  INTERFACE

    SUBROUTINE     TL_ALLOC
    END SUBROUTINE TL_ALLOC

  END INTERFACE

!-----------------------------------------------------------------------

  CALL TL_ALLOC

  DO ITILE = MT1, MT2

#if(USE_MPI)
    IF ( myid == root ) THEN
#endif
      CALL CPU_TIME ( time_end )
      CALL DATE_AND_TIME ( BIG_BEN(1),BIG_BEN(2),BIG_BEN(3), DATE_TIME )
      PRINT     *, 'CPU time: ', time_end - time_begin, ' seconds'
      WRITE (*,'(''Current time: '',I2,'':'',I2,'':'',F6.3)')            &
      DATE_TIME(5), DATE_TIME(6), DATE_TIME(7)+REAL(DATE_TIME(8)/1000.)
      PRINT     *, 'time level:', ITILE
      WRITE (11,'(I6,1X,F10.2,1X,I4,I3,I3,I5,3I3,I4)')                   &
                ITILE, time_end - time_begin, DATE_TIME

      INQUIRE ( FILE = 'STOP', EXIST = EX )
      IF ( EX ) THEN
        MT2 = ITILE-1
        OPEN  ( 35, FILE='STOP', STATUS='UNKNOWN')
        WRITE ( 35,* ) ITILE
        CLOSE ( 35 )
        STOP
      END IF

#if(USE_MPI)
    END IF
#endif
!-- UPDATE STRESS FIELD

    DO L = MAX(Z_MIN, 2), MIN (Z_MAX, LIN_Z-1)
      CALL  HPLANE_IC_STRESS ( L )
    END DO
    
    
    IF ( ( Z_MIN <= LIN_Z   ) .AND. ( LIN_Z   <= Z_MAX ) ) THEN
       CALL  HPLANE_LIN_STRESS ( SOURTF_LIN = SOURTF_A(1,ITILE), &            
                                 SOURTF_LOU = SOURTF_A(2,ITILE) )
    END IF
    
    
    DO L = MAX(Z_MIN, LIN_Z+1), MIN (Z_MAX, MZ-2)
      CALL  HPLANE_IC_STRESS ( L )                                             
    END DO

    IF ( Z_MIN == 0 ) THEN
      CALL  HPLANE_0_STRESS                                                    
      CALL  HPLANE_1_STRESS                                                   
    END IF

    IF ( ( Z_MIN <= MZ-1 ) .AND. ( MZ-1 <= Z_MAX ) ) THEN
      CALL  HPLANE_MZ1_STRESS
    END IF

    DO L = MAX(Z_MIN, MZ), MIN(Z_MAX, PMLZH1) 
        
      CALL  HPLANE_PML_STRESS ( L, L-MZ+1 ) 
    END DO

    IF ( Z_MAX == PMLZH ) THEN
      CALL  HPLANE_PMLZH0_STRESS
    END IF
!--- COMMUNICATE STRESSFIELD

#IF(USE_MPI)
    CALL COMMUNICATE_STRESSFIELD
#ENDIF

!--- STORE BOTTOM BOUNDARY  

    IF ( Z_MAX == PMLZH ) THEN       
      SZ_UP (:)         = SZ_U  (:)       
      SZ_WP (:)         = SZ_W  (:)       
      SZ_QXP(:)         = SZ_QX (:)
      SZ_QZP(:)         = SZ_QZ (:)

      SZ_U  (:)%TIKMZ2  = UM( :, PMLZH2 )
      SZ_W  (:)%TIKMZ2  = WM( :, PMLZH2 )
      SZ_QX (:)%TIKMZ2  = QX( :, PMLZH2 )
      SZ_QZ (:)%TIKMZ2  = QZ( :, PMLZH2 )

      SZ_U  (:)%TIKMZ1  = UM( :, PMLZH1 )
      SZ_W  (:)%TIKMZ1  = WM( :, PMLZH1 )
      SZ_QX (:)%TIKMZ1  = QX( :, PMLZH1)
      SZ_QZ (:)%TIKMZ1  = QZ( :, PMLZH1)
      

      SZ_U  (:)%TIKMZ0  = UM( :, PMLZH  )
      SZ_W  (:)%TIKMZ0  = WM( :, PMLZH  )
      SZ_QX (:)%TIKMZ0  = QX( :, PMLZH  )
      SZ_QZ (:)%TIKMZ0  = QZ( :, PMLZH  )
      
    END IF

!--------------------------

    DO L = MAX(Z_MIN, 2), MIN (Z_MAX, MZ-2)
      CALL  HPLANE_IC_VELOCITY ( L )
    END DO

    IF ( Z_MIN == 0 ) THEN
      IF (STRESS_IMAGING) THEN
        CALL  HPLANE_1_VELOCITY
      ELSE
        CALL  HPLANE_0_VELOCITY
        CALL  HPLANE_1_VELOCITY
      END IF
    END IF

    DO L = MAX(Z_MIN, MZ-1), MIN(Z_MAX, PMLZH1)
        
      CALL  HPLANE_PML_VELOCITY ( L, L-MZ+1 )
    END DO

    IF ( Z_MAX == PMLZH ) THEN
      CALL  HPLANE_PMLZH0_VELOCITY            
    END IF
!--- COMMUNICATE WAVEFIELD

#if(USE_MPI)
    CALL COMMUNICATE_WAVEFIELD
#endif

!--- SYNCHRONIZATION

#if(USE_MPI)
    CALL MPI_BARRIER (comm_cart,errcode)
#endif

!--- STORE RESULTS

    IF ( KEY_TLV ) THEN
      CALL STORE_SEISMOGRAMS(ITILE)
    END IF

    IF ( KEY_TLP ) THEN
      CALL STORE_SEISMOGRAMS_PRESSURE(ITILE)
    END IF

    IF ( KEY_SNV .AND. (MOD(ITILE,IPAS2)==0) ) THEN
      CALL STORE_SNAPSHOTS(ITILE)
    END IF

  END DO


END  SUBROUTINE  TIME_LOOP
