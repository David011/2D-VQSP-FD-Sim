!=======================================================================

SUBROUTINE  OPEN_SR_FILE

  USE PRECISION   , ONLY: PP, DP

  USE CONTROL_DATA, ONLY: SR_FILE_NAME, KEY_SOUR,                      &
                          MT1, MT2, DT, LIN_X, LIN_Z, STRESS_IMAGING, QR

  USE GRID_MEDIUM , ONLY:                                              &
                          JM,                                          &
                          MU,                                          &
                          SIG_XX_1, SIG_ZZ_2, PRES_1, PRES_2, PRES_3,  &
                          REL1_U,   REL2_QX,  REL3_U,                  &
                          REL1_W,   REL2_QZ,  REL3_W,                  &
                          H

  USE AUXIL       , ONLY: PI, SOURTF_A, X_MIN, X_MAX, Z_MIN, Z_MAX
  
  USE INTERFACES  , ONLY: CALCUL_VELOCITY
  
#if(USE_MPI)
#if(MPI2)
  USE MPI
#endif
  USE MPI_SHARED
#endif

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER       ::  I, NSTFS, IOS, JM1

  REAL (KIND=PP) :: TMP, VELX, VEL, VELS, VELPF, VELPS, T

  REAL (KIND=PP), DIMENSION(:), ALLOCATABLE ::  TIME, STF, STF2

!-----------------------------------------------------------------------

  ALLOCATE( SOURTF_A(2,MT2) )
  SOURTF_A = 0.

  IF ( ( Z_MIN <=LIN_Z ) .AND. ( LIN_Z <= Z_MAX ) .AND. ( X_MIN <=LIN_X ) .AND. ( LIN_X <= X_MAX ) ) THEN
         JM1 = JM (LIN_X, LIN_Z)
         CALL CALCUL_VELOCITY( SIG_ZZ_2(JM1), - PRES_2(JM1), - PRES_3(JM1), MU(JM1), H*REL1_W(JM1)/DT, H*REL2_QZ(JM1)/DT, H*REL3_W(JM1)/DT, VELPF , VELPS, VELS )
         VELX = VELPF
  ELSE
         VELX = 0.
  END IF

#if(USE_MPI)
  CALL MPI_REDUCE                                                      &
             (VELX ,VEL, 1, MPI_REAL, MPI_MAX, root, comm_cart, errcode)
  IF ( myid == root ) THEN
#else
  VEL = VELX
#endif

  
    OPEN  ( 10, FILE = TRIM(SR_FILE_NAME), STATUS='OLD' )

    I = 0
    DO
      READ(10,*,IOSTAT=IOS) TMP, TMP
      IF ( IOS /= 0 ) EXIT
      I = I+1
    END DO
    REWIND (10)
    NSTFS = I+1
    ALLOCATE ( TIME(NSTFS), STF(NSTFS), STF2(NSTFS) )

    DO I = 1, NSTFS-1
      READ(10,*) TIME(I), STF(I)
    END DO
    CLOSE(10)
   
    TIME(NSTFS) = MT2*DT
    STF (NSTFS) = STF(NSTFS-1)

    CALL spline( TIME, STF, 0., 0., STF2)

     DO I = MT1, MT2
        
      T = (I-1)*DT
      IF (T > TIME(NSTFS-1)) EXIT                !zero padding
      SOURTF_A (1,I) = splint(TIME,STF,STF2,T)
      
      SELECT CASE ( KEY_SOUR )
      !                       bulk source
      CASE ( 0 )
          
      SOURTF_A (2,I) = QR*SOURTF_A (1,I)/(H**2)
      SOURTF_A (1,I) = (1-QR)*SOURTF_A (1,I)/(H**2)
      !                       solid source
      CASE ( 1 )
      SOURTF_A (2,I) = 0.
      SOURTF_A (1,I) = SOURTF_A (1,I)/(H**2)
      !                       fluid source
      CASE ( 2 )
      SOURTF_A (2,I) = SOURTF_A (1,I)/(H**2)
      SOURTF_A (1,I) = QR*SOURTF_A (1,I)/(H**2)
      !                       unit fluid source
      CASE ( 3 )
      SOURTF_A (2,I) = SOURTF_A (1,I)/(H**2)
      SOURTF_A (1,I) = 0.
            
      END SELECT
      
    END DO
     
#if(USE_MPI)
  END IF

  CALL MPI_Bcast  (SOURTF_A, MT2*2, MPI_REAL, root, comm_cart, errcode)
  CALL MPI_Barrier(comm_cart,errcode)
#endif

CONTAINS

  SUBROUTINE spline(x,y,yp1,ypn,y2)
    USE PRECISION, ONLY: PP
    IMPLICIT NONE
    REAL(PP), DIMENSION(:), INTENT(IN)  :: x,y
    REAL(PP),               INTENT(IN)  :: yp1,ypn
    REAL(PP), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER :: n
    REAL(PP), DIMENSION(size(x)) :: a,b,c,r
    n=size(x)
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0*(c(2:n-1)+a(2:n-1))
    b(1)=1.0
    b(n)=1.0
    if (yp1 > 0.99e30) then
      r(1)=0.0
      c(1)=0.0
    else
      r(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      c(1)=0.5
    end if

  if (ypn > 0.99e30) then
      r(n)=0.0
      a(n)=0.0
    else
      r(n)=(-3.0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
      a(n)=0.5
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))

  END SUBROUTINE spline

  SUBROUTINE tridag(a,b,c,r,u)
    USE PRECISION, ONLY: PP
    IMPLICIT NONE
    REAL(PP), DIMENSION(:), INTENT(IN)  :: a,b,c,r
    REAL(PP), DIMENSION(:), INTENT(OUT) :: u
    REAL(PP), DIMENSION(size(b))        :: gam
    INTEGER :: n,j
    REAL(PP) :: bet
    n=size(b)
    bet=b(1)
    if (bet == 0.0) STOP 'tridag_ser: Error at code stage 1'
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j-1)*gam(j)
      if (bet == 0.0) STOP 'tridag_ser: Error at code stage 2'
      u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag

  FUNCTION splint(xa,ya,y2a,x)
    USE PRECISION, ONLY: PP
    IMPLICIT NONE
    REAL(PP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(PP), INTENT(IN) :: x
    REAL(PP) :: splint
    INTEGER :: khi,klo,n
    REAL(PP) :: a,b,h
    n=size(xa)
    klo=max(min(locate(xa,x),n-1),1)
    khi=klo+1
    h=xa(khi)-xa(klo)
    if (h == 0.0) STOP 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0

  END FUNCTION splint

  FUNCTION locate(xx,x)
    USE PRECISION, ONLY: PP
    IMPLICIT NONE
    REAL(PP), DIMENSION(:), INTENT(IN) :: xx
    REAL(PP)              , INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
      if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd .eqv. (x >= xx(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do
    if (x == xx(1)) then
      locate=1
    else if (x == xx(n)) then
      locate=n-1
    else
      locate=jl
    end if
  END FUNCTION locate


END SUBROUTINE  OPEN_SR_FILE
