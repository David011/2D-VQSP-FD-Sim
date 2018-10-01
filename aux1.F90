!=======================================================================

SUBROUTINE  AUX1
  USE PRECISION   , ONLY :                                             &
                           PP

  USE GRID_MEDIUM , ONLY :                                             &
                           H

  USE CONTROL_DATA, ONLY :                                             &
                           MX, MZ, DT, TPML

  USE AUXIL       , ONLY :                                             &
                           PI,                                         &
                           A, B, AH, BH, CH,                           & 
                           PMLL, PMLXH ,         PMLZH, PMLZH1, PMLZH2,& 
                           MX1                                           

#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER        :: J

  REAL (KIND=PP) :: OH

!-----------------------------------------------------------------------

  MX1 = MX - 1

  PMLL   =   1 - TPML               
  PMLXH  = MX  + TPML
  PMLZH  = MZ  + TPML  
  PMLZH1 = MZ  + TPML - 1
  PMLZH2 = MZ  + TPML - 2

  AH  = A    /H
  BH  = B    /H
  CH  = 1._PP/H

END SUBROUTINE  AUX1
