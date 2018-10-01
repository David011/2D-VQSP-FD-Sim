!_______________________________________________________________________

MODULE NONREF_BOUND

  USE PRECISION, ONLY: PP

!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------

  TYPE NR_PAR_A_X                                                          
    REAL (KIND=PP) ::                                                  &
                    ALE01,   ALE02,   ALE10,   ALE11,                  &
                    ALE12,   ALE20,   ALE21,   ALE22,                  &   
                    ARI01,   ARI02,   ARI10,   ARI11,                  &
                    ARI12,   ARI20,   ARI21,   ARI22
  END TYPE NR_PAR_A_X

  TYPE NR_DATA_Z
    REAL (KIND=PP) ::                                                  &
                    TIKMZ0  , TIKMZ1  , TIKMZ2                          
                                                                        
                                                                        
    
  END TYPE NR_DATA_Z

  TYPE NR_DATA_X
    REAL (KIND=PP) ::                                                  &
                    T1KL    , T2KL    , T3KL    ,                      &
                    TMX0KL  , TMX1KL  , TMX2KL
  END TYPE NR_DATA_X

!-----------------------------------------------------------------------

  TYPE (NR_PAR_A_X), DIMENSION (:), ALLOCATABLE ::   AXU, AXW                

  TYPE (NR_DATA_X), DIMENSION (:),  ALLOCATABLE ::  SX_U, SX_W, SX_QX, SX_QZ

  TYPE (NR_DATA_Z), DIMENSION (:),  ALLOCATABLE ::  SZ_U,  SZ_W,       &   
                                                    SZ_UP, SZ_WP,      &   
                                                    SZ_QX,  SZ_QZ,     &
                                                    SZ_QXP, SZ_QZP
  
  TYPE (NR_DATA_X) ::  SX_UP, SX_WP, SX_QXP, SX_QZP

  REAL (KIND=PP)   ::                                                  &
                    ABW01,   ABW02,   ABW10,   ABW11,                  &
                    ABW12,   ABW20,   ABW21,   ABW22,                  &
                    ABH01,   ABH02,   ABH10,   ABH11,                  &
                    ABH12,   ABH20,   ABH21,   ABH22

END MODULE NONREF_BOUND
