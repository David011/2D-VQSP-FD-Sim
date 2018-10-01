SUBROUTINE CALCUL_VELOCITY ( SIG_1, P,    AINV, MU,  &
                             REL1,  REL2, REL3,      &
                             VPF ,   VPS, VS         )
    USE PRECISION   , ONLY:                                              &
                          PP
    
    IMPLICIT NONE
    
    REAL (KIND=PP), INTENT (IN) :: SIG_1, P, AINV, MU
    REAL (KIND=PP), INTENT (IN) :: REL1, REL2, REL3
    REAL (KIND=PP), INTENT (OUT) :: VPF, VPS, VS
    REAL (KIND=PP), DIMENSION(2,2) :: AINVER, B, TT
    REAL (KIND=PP) :: DELTA
    REAL (KIND=PP), DIMENSION(1:2) :: D
    
    AINVER(1,1) =   REL1
    AINVER(2,2) =   REL2
    AINVER(1,2) = - REL3
    AINVER(2,1) = - REL3
    
    B(1,1) = SIG_1
    B(2,2) = AINV
    B(1,2) = P
    B(2,1) = P
    
    TT    = MATMUL(AINVER, B)
    DELTA = (TT(1,1) + TT(2,2))*(TT(1,1) + TT(2,2)) + 4._PP*(TT(1,2)*TT(2,1) - TT(1,1)*TT(2,2))
    D(1)  = ((TT(1,1) + TT(2,2)) + SQRT(DELTA))/2._PP
    D(2)  = ((TT(1,1) + TT(2,2)) - SQRT(DELTA))/2._PP
    
    VPF = SQRT(D(1))
    VPS = SQRT(D(2))
    VS  = SQRT(REL1*MU)
    IF (MU == 0.) THEN
        VS  = VPF
    END IF
    
END SUBROUTINE CALCUL_VELOCITY