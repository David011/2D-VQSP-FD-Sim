!=======================================================================

SUBROUTINE ALLOC_NR

  USE PRECISION   , ONLY:                                              &
                          PP
  USE CONTROL_DATA, ONLY:                                              &
                          MX, MZ

  USE NONREF_BOUND, ONLY:                                              &
                          NR_DATA_X,            NR_DATA_Z,             &    
                                   AXU, AXW,                           &   
                                   SX_U , SX_UP , SX_QX , SX_QXP ,     &   
                                   SX_W , SX_WP , SX_QZ , SX_QZP ,     &    
                                   SZ_U , SZ_UP , SZ_QX , SZ_QXP ,     &   
                                   SZ_W , SZ_WP , SZ_QZ , SZ_QZP            

  USE AUXIL       , ONLY:                                              &
                          X_MIN_D, X_MAX_D, Z_MIN_D, Z_MAX_D

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  ALLOCATE (  AXU   ( Z_MIN_D:Z_MAX_D ),                               &
              AXW   ( Z_MIN_D:Z_MAX_D ) )
  ALLOCATE (  SX_U  ( Z_MIN_D:Z_MAX_D ),                               &
              SX_W  ( Z_MIN_D:Z_MAX_D ),                               &
              SX_QX ( Z_MIN_D:Z_MAX_D ),                               &
              SX_QZ ( Z_MIN_D:Z_MAX_D ),                               &
              SZ_U  ( X_MIN_D:X_MAX_D ),                               &
              SZ_W  ( X_MIN_D:X_MAX_D ),                               &
              SZ_UP ( X_MIN_D:X_MAX_D ),                               &
              SZ_WP ( X_MIN_D:X_MAX_D ),                               &
              SZ_QX ( X_MIN_D:X_MAX_D ),                               &
              SZ_QZ ( X_MIN_D:X_MAX_D ),                               &
              SZ_QXP ( X_MIN_D:X_MAX_D ),                              &
              SZ_QZP ( X_MIN_D:X_MAX_D ) )
  
  SX_U   = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )    
  SX_W   = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )
  SX_QX  = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )    
  SX_QZ  = NR_DATA_X ( 0._PP, 0._PP, 0._PP, 0._PP, 0._PP, 0._PP )

  SZ_U   = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_W   = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_UP  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_WP  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_QX  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_QZ  = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_QXP = NR_DATA_Z ( 0._PP, 0._PP, 0._PP )
  SZ_QZP = NR_DATA_Z ( 0._PP, 0._PP, 0._PP ) 

END SUBROUTINE ALLOC_NR
