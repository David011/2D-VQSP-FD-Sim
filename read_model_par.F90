!=======================================================================

SUBROUTINE  READ_MODEL_PAR

  USE PRECISION   , ONLY:                                              &
                          PP
  USE GRID_MEDIUM , ONLY:                                              &
                           JMNUM,                                      & 
                           H,                                          &
                           MU,                                         &
                           SIG_XX_1,       SIG_XX_2,       SIG_XX_3,   &
                           SIG_ZZ_1,       SIG_ZZ_2,       SIG_ZZ_3,   &
                           PRES_1,         PRES_2,         PRES_3,     &
                           REL1_U,         REL2_U,         REL3_U,     &
                           REL1_W,         REL2_W,         REL3_W,     &
                           REL1_QX,        REL2_QX,        REL3_QX,    &
                           REL1_QZ,        REL2_QZ,        REL3_QZ,    &
                           AUXIL_W,        AUXIL_QZ
  USE CONTROL_DATA, ONLY:                                              &
                          DT

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: J, JM1, ALLOSTAT, JMNUMR, I

!-----------------------------------------------------------------------

    READ ( 14 ) JMNUMR

    JMNUM = ABS(JMNUMR)

  ALLOCATE ( MU(JMNUM),    SIG_XX_1(JMNUM),    SIG_XX_2(JMNUM),    SIG_XX_3(JMNUM),                          &
                           SIG_ZZ_1(JMNUM),    SIG_ZZ_2(JMNUM),    SIG_ZZ_3(JMNUM),                          &
                           PRES_1  (JMNUM),    PRES_2  (JMNUM),    PRES_3  (JMNUM),                          &
                           REL1_U  (JMNUM),    REL1_W  (JMNUM),    REL2_U  (JMNUM),                          &
                           REL2_W  (JMNUM),    REL3_U  (JMNUM),    REL3_W  (JMNUM),                          &
                           REL1_QX (JMNUM),    REL1_QZ (JMNUM),    REL2_QX (JMNUM),                          &
                           REL2_QZ (JMNUM),    REL3_QX (JMNUM),    REL3_QZ (JMNUM),                          &
                           AUXIL_W (JMNUM),    AUXIL_QZ(JMNUM)                                               )                                                               

    

   READ  ( 14 )                                                                                              &
           ( MU(JM1),      SIG_XX_1(JM1)  ,    SIG_XX_2(JM1)  ,    SIG_XX_3(JM1)  ,                          &
                           SIG_ZZ_1(JM1)  ,    SIG_ZZ_2(JM1)  ,    SIG_ZZ_3(JM1)  ,                          &
                           PRES_1  (JM1)  ,    PRES_2  (JM1)  ,    PRES_3  (JM1)  ,                          &
                           REL1_U  (JM1)  ,    REL1_W  (JM1)  ,    REL2_U  (JM1)  ,                          &
                           REL2_W  (JM1)  ,    REL3_U  (JM1)  ,    REL3_W  (JM1)  ,                          &
                           REL1_QX (JM1)  ,    REL1_QZ (JM1)  ,    REL2_QX (JM1)  ,                          &
                           REL2_QZ (JM1)  ,    REL3_QX (JM1)  ,    REL3_QZ (JM1)  ,                          &
                           AUXIL_W (JM1)  ,    AUXIL_QZ(JM1)  ,    JM1 = 1, JMNUM                            )
        
   
  REL1_U  =  DT*REL1_U/H
  REL2_U  =  DT*REL2_U
  REL3_U  =  DT*REL3_U/H
  REL1_W  =  DT*REL1_W/H
  REL2_W  =  DT*REL2_W
  REL3_W  =  DT*REL3_W/H
  REL1_QX =  DT*REL1_QX/H
  REL2_QX =  DT*REL2_QX/H
  REL3_QX =  DT*REL3_QX
  REL1_QZ =  DT*REL1_QZ/H
  REL2_QZ =  DT*REL2_QZ/H
  REL3_QZ =  DT*REL3_QZ

  

      !DO JM1 = 1, JMNUM
         !IF (CAP_M(JM1) == 0) THEN
             !AUXIL_W (JM1) = LAMC(JM1)/L2MC(JM1)
             !AUXIL_QZ(JM1) = 0._PP
         !ELSE
             !AUXIL_W(JM1)  = (LAMC(JM1) -  (ALPHA_M(JM1) * ALPHA_M(JM1) / CAP_M(JM1)))/(L2MC(JM1) - (ALPHA_M(JM1) * ALPHA_M(JM1) / CAP_M(JM1)))
             !AUXIL_QZ(JM1) = ((ALPHA_M(JM1) / CAP_M(JM1)) * 2. * MU(JM1)) / (L2MC(JM1) - (ALPHA_M(JM1) * ALPHA_M(JM1) / CAP_M(JM1)))
             
         !END IF
      !END DO
  
  CLOSE(14)

END SUBROUTINE  READ_MODEL_PAR
