!=======================================================================
! 2-D UPWINDING
! VAN LEER FOR THE THIN LAYER GAS FLOW
! WITH AUSMD AND RNG K-E TURBULENCE MODEL
! PREDICITE PARTICULATE TRAJECTORY BY LAGRANGIAN METHOD
!=======================================================================

SUBROUTINE solver(Qin,alpha,beta,imethod,Qout)

  IMPLICIT NONE

  REAL*8, intent(in) :: Qin(6,-1:1201,-1:251)
  REAL*8, intent(in) :: alpha,beta,imethod

  REAL*8, intent(out) :: Qout(6,-1:1201,-1:251)

!-----------------------------------------------------------------------

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

!-----------------------------------------------------------------------

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N

!-----------------------------------------------------------------------

  CALL QTERM(Qin)

  CALL STATETFUN(Qin)

  CALL sour1(Qin)

  CALL EXPLICIT(Qin)

  DO K=1,KXM1
    
    CALL ETERM1(K)

    CALL MUSCL(JMAX)

    CALL THINCEM(beta,JMAX)

    if (imethod == 1) then

      CALL ORI(alpha,JMAX)

    elseif (imethod == 2) then

      CALL ATM(alpha,JMAX)

    elseif (imethod == 3) then

      CALL JCP(JMAX)

    elseif (imethod == 4) then

      CALL PMUS(JMAX)

    elseif (imethod == 5) then

      CALL PTHI(JMAX)

    endif

    CALL ETERM2(K)

  ENDDO

  DO J=1,JXM1

    CALL FTERM1(J)

    CALL MUSCL(KMAX)

    CALL THINCEM(beta,KMAX)

    if (imethod == 1) then

      CALL ORI(alpha,KMAX)

    elseif (imethod == 2) then

      CALL ATM(alpha,KMAX)

    elseif (imethod == 3) then

      CALL JCP(KMAX)

    elseif (imethod == 4) then

      CALL PMUS(KMAX)

    elseif (imethod == 5) then

      CALL PTHI(KMAX)

    endif

    CALL FTERM2(J)

  ENDDO

  CALL Qoo1(Qin)

  CALL EXPLICIT(Qin)

  DO K=1,KXM1
    
    CALL ETERM1(K)

    CALL MUSCL(JMAX)

    CALL THINCEM(beta,JMAX)

    if (imethod == 1) then

      CALL ORI(alpha,JMAX)

    elseif (imethod == 2) then

      CALL ATM(alpha,JMAX)

    elseif (imethod == 3) then

      CALL JCP(JMAX)

    elseif (imethod == 4) then

      CALL PMUS(JMAX)

    elseif (imethod == 5) then

      CALL PTHI(JMAX)

    endif

    CALL ETERM2(K)

  ENDDO

  DO J=1,JXM1

    CALL FTERM1(J)

    CALL MUSCL(KMAX)

    CALL THINCEM(beta,KMAX)

    if (imethod == 1) then

      CALL ORI(alpha,KMAX)

    elseif (imethod == 2) then

      CALL ATM(alpha,KMAX)

    elseif (imethod == 3) then

      CALL JCP(KMAX)

    elseif (imethod == 4) then

      CALL PMUS(KMAX)

    elseif (imethod == 5) then

      CALL PTHI(KMAX)

    endif

    CALL FTERM2(J)

  ENDDO

  CALL Qoo2(Qin)

  CALL sour2(Qin)

  CALL LOWERWALL(Qin)

  CALL UPWALL(Qin)

  CALL ghostpoint(Qin)

  CALL controlz(Qin)

   DO N=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Qout(N,J,K) = Qin(N,J,K)
      ENDDO
    ENDDO
  ENDDO

RETURN

END

!=======================================================================
! 2-D UPWINDING
! VAN LEER FOR THE THIN LAYER GAS FLOW
! WITH AUSMD AND RNG K-E TURBULENCE MODEL
! PREDICITE PARTICULATE TRAJECTORY BY LAGRANGIAN METHOD
!=======================================================================

SUBROUTINE plot_bestresult(alpha,beta,ngrd,imethod,ISAVE,ITMAX)

  IMPLICIT NONE

  REAL*8, intent(in) :: alpha(ITMAX)
  REAL*8, intent(in) :: beta(ITMAX)

  INTEGER, intent(in) :: ITMAX,ISAVE,ngrd,imethod

!-----------------------------------------------------------------------

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

!-----------------------------------------------------------------------

  REAL*8 :: Q(6,-1:1201,-1:251)

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,IT

!-----------------------------------------------------------------------

  CALL inic(Q,ngrd)

  CALL QTERM(Q)

  DO IT=1,ITMAX

  CALL STATETFUN(Q)

  CALL sour1(Q)

  CALL EXPLICIT(Q)

  DO K=1,KXM1
    
      CALL ETERM1(K)

      CALL MUSCL(JMAX)

      CALL THINCEM(beta(IT),JMAX)

      if (imethod == 1) then

        CALL ORI(alpha(IT),JMAX)

      elseif (imethod == 2) then

        CALL ATM(alpha(IT),JMAX)

      elseif (imethod == 3) then

        CALL JCP(JMAX)

      elseif (imethod == 4) then

        CALL PMUS(JMAX)

      elseif (imethod == 5) then

        CALL PTHI(JMAX)

      endif

      CALL ETERM2(K)

    ENDDO

    DO J=1,JXM1

      CALL FTERM1(J)

      CALL MUSCL(KMAX)

      CALL THINCEM(beta(IT),KMAX)

      if (imethod == 1) then

        CALL ORI(alpha(IT),KMAX)

      elseif (imethod == 2) then

        CALL ATM(alpha(IT),KMAX)

      elseif (imethod == 3) then

        CALL JCP(KMAX)

      elseif (imethod == 4) then

        CALL PMUS(KMAX)

      elseif (imethod == 5) then

        CALL PTHI(KMAX)

      endif

      CALL FTERM2(J)

    ENDDO

    CALL Qoo1(Q)

    CALL EXPLICIT(Q)

    DO K=1,KXM1
      
      CALL ETERM1(K)

      CALL MUSCL(JMAX)

      CALL THINCEM(beta(IT),JMAX)

      if (imethod == 1) then

        CALL ORI(alpha(IT),JMAX)

      elseif (imethod == 2) then

        CALL ATM(alpha(IT),JMAX)

      elseif (imethod == 3) then

        CALL JCP(JMAX)

      elseif (imethod == 4) then

        CALL PMUS(JMAX)

      elseif (imethod == 5) then

        CALL PTHI(JMAX)

      endif

      CALL ETERM2(K)

    ENDDO

    DO J=1,JXM1

      CALL FTERM1(J)

      CALL MUSCL(KMAX)

      CALL THINCEM(beta(IT),KMAX)

      if (imethod == 1) then

        CALL ORI(alpha(IT),KMAX)

      elseif (imethod == 2) then

        CALL ATM(alpha(IT),KMAX)

      elseif (imethod == 3) then

        CALL JCP(KMAX)

      elseif (imethod == 4) then

        CALL PMUS(KMAX)

      elseif (imethod == 5) then

        CALL PTHI(KMAX)

      endif

      CALL FTERM2(J)

    ENDDO

    CALL Qoo2(Q)

    CALL sour2(Q)

    CALL LOWERWALL(Q)

    CALL UPWALL(Q)

    CALL ghostpoint(Q)

    CALL controlz(Q)

    CALL SavePlt(Q,IT,ISAVE)

  enddo

RETURN

END

!=======================================================================
! Subroutine GRD
! readgrdreadgrd
!=======================================================================

SUBROUTINE readgrd(ngrd)

  IMPLICIT NONE

  REAL*8 :: X,Y
  REAL*8 :: xini,xend,yini,yend,dx,dy
  
  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: mx,i
  INTEGER :: my,j
  INTEGER :: ngrd

  COMMON/XY/X(0:1201,0:251),Y(0:1201,0:251)

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  xini = 0.0
  xend = 0.025

  yini = 0.0
  yend = 0.005

  mx = 400*ngrd
  my = 80*ngrd

  dx = (xend-xini) / mx
  dy = (yend-yini) / my

  do j = 1, my+1
    do i = 1, mx+1
      x(i,j) = xini + (i-1)*dx
      y(i,j) = yini + (j-1)*dy
    enddo
  enddo

  JMAX=mx+1
  KMAX=my+1
  JMP1=JMAX+1
  JXM1=JMAX-1
  KMP1=KMAX+1
  KXM1=KMAX-1

RETURN

END

!=======================================================================
! Subroutine Initial Condition
! Chan
!=======================================================================

SUBROUTINE inic(Q,ngrd)

  IMPLICIT NONE

  INTEGER, intent(in) :: ngrd

  REAL*8, intent(out) :: Q(6,-1:1201,-1:251)

  REAL*8 :: X,Y
  REAL*8 :: akp1,aky,ajp1,ajx
  REAL*8 :: Q1Z,U,VINF,PZ,Z
  REAL*8 :: Q1Z2,U2,PZ2,Z2
  REAL*8 :: q0,r

  INTEGER :: JMAX,KMAX,JMP1,JXM1,KMP1,KXM1
  INTEGER :: J,K

  COMMON/XY/X(0:1201,0:251),Y(0:1201,0:251)

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  CALL readgrd(ngrd)

  q0=0.5196d10
  r=1.4

  akp1=-2.0

  DO  K=-1,kMP1 

    akp1=akp1+1.0
    aky=dabs( akp1*y(1,2) -0.0025)

    if( aky .GE. 0.001 )then
      aky=0.004
    elseif(aky .LT. 0.001)then
      aky=0.005-aky
    endif

    ajp1=-2.0

    DO  J=-1,JMP1 
      ajp1=ajp1+1.0
      ajx=ajp1*y(1,2)

      if(ajx .LE. aky)then
        q1z=1.945d-3
        U=8.162d4
        VINF=0.0
        pz=6.27d6
        z=0.0

        Q(1,J,K)=q1z
        Q(2,J,K)=q1z*U
        Q(3,J,K)=q1z*VINF
        Q(4,J,K)=pz/(r-1)+0.5*q1z*(U**2+VINF**2)+q0*q1z*z
        Q(5,J,K)=q1z*z

      elseif(ajx .GT. aky)then
        q1z2=1.201d-3
        U2=0.0 
        VINF=0.0
        pz2=8.321d5
        z2=1.0

        Q(1,J,K)=q1z2
        Q(2,J,K)=q1z2*U2
        Q(3,J,K)=q1z2*VINF
        Q(4,J,K)=pz2/(r-1)+0.5*q1z2*(U2**2+VINF**2)+q0*q1z2*z2
        Q(5,J,K)=q1z2*z2
      endif

    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine SET UP INDEX
! QTERM
!=======================================================================

SUBROUTINE QTERM(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)  

  REAL*8 :: X,Y,RJ,VOL
  REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS
  REAL*8 :: XIX,XIY,ETX,ETY

  INTEGER :: JMAX,KMAX,JMP1,JXM1,KMP1,KXM1
  INTEGER :: N,K,J,JP1,KP1

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/JACO/RJ(-1:1201,-1:251)

  COMMON/XY/X(0:1201,0:251),Y(0:1201,0:251)

  COMMON/Qold/Qold1ds(-1:1201,-1:251,5),Qold1s(-1:1201,-1:251,5),&
              Qold1d(-1:1201,-1:251,5),Qold1(-1:1201,-1:251,5),&
              Qold2d(-1:1201,-1:251,5),Qold2ss(-1:1201,-1:251,5),&
              Qolddss(-1:1201,-1:251,5)

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

!-----------------------------------------------------------------------

  do n=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Qold1ds(J,K,n) = Q(n,J,K)
        Qold1s(J,K,n)  = Q(n,J,K)
        Qold1d(J,K,n)  = Q(n,J,K)
        Qold1(J,K,n)   = Q(n,J,K)
        Qold2d(J,K,n)  = Q(n,J,K)
        Qold2ss(J,K,n) = Q(n,J,K)
        Qolddss(J,K,n) = Q(n,J,K)
      enddo
    enddo
  enddo

!-----------------ETX,ETY-----------------------------------------------

  DO K=1,KMAX
    DO J=1,JXM1
      JP1=J+1
      ETX(J,K)=-Y(JP1,K)+Y(J,K)
      ETY(J,K)=X(JP1,K)-X(J,K)
      ETX(0,K)=ETX(1,K)
      ETY(0,K)=ETY(1,K)
      ETX(JMAX,K)=ETX(JXM1,K)
      ETY(JMAX,K)=ETY(JXM1,K)
    enddo
  enddo

!-----------------XIX,XIY-----------------------------------------------     

  DO K=1,KXM1
    DO J=1,JMAX
      KP1=K+1
      XIX(J,K)=Y(J,KP1)-Y(J,K)
      XIY(J,K)=-X(J,KP1)+X(J,K)
    enddo
  enddo

  DO J=1,JMAX
    XIX(J,0)=XIX(J,1)
    XIY(J,0)=XIY(J,1)
    XIX(J,KMAX)=XIX(J,KXM1)
    XIY(J,KMAX)=XIY(J,KXM1)
  enddo

  DO K=1,KMAX
    XIX(0,K)=XIX(1,K)
    XIY(0,K)=XIY(1,K)
    XIX(JMAX,K)=XIX(JXM1,K)
    XIY(JMAX,K)=XIY(JXM1,K)
  enddo

!-----------------RJ=1./VOL---------------------------------------------        

  DO K=1,KXM1
    DO J=1,JXM1
      KP1=K+1
      JP1=J+1
      VOL=((X(J,K)-X(JP1,KP1))*(Y(JP1,K)-Y(J,KP1))-(X(JP1,K)-X(J,KP1))*(Y(J,K)-Y(JP1,KP1)))/2.
      RJ(J,K)=1./VOL
    enddo
  enddo

  DO J=1,JXM1
    RJ(J,0)=RJ(J,1)
    RJ(J,-1)=RJ(J,0)
    RJ(J,KMAX)=RJ(J,KXM1)
    RJ(J,KMP1)=RJ(J,KMAX)
  enddo

  DO K=1,KXM1
    RJ(0,K)=RJ(1,K)
    RJ(-1,K)=RJ(0,K)
    RJ(JMAX,K)=RJ(JXM1,K)
    RJ(JMP1,K)=RJ(JMAX,K)
  enddo

RETURN

END

!=======================================================================
! Subroutine STATEMENT FUNCTION
! STATETFUN
!=======================================================================

SUBROUTINE STATETFUN(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  real*8 :: DT,RJ
  REAL*8 :: XIXC,XIYC,ETXC,ETYC
  REAL*8 :: XIX,XIY,ETX,ETY
  REAL*8 :: UF,VF,PF,SONF
  real*8 :: TOT,DTMIN
  REAL*8 :: CFL,GM,q0

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: J,K

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Dtt/DT(0:1201,0:251)

  COMMON/JACO/RJ(-1:1201,-1:251)

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/ST/UF(0:1201,0:251),VF(0:1201,0:251),PF(0:1201,0:251),SONF(0:1201,0:251)

  COMMON/TIME/TOT,DTMIN

  CFL=0.1
  q0=0.5196d10
  GM=1.4

  DTMIN=100.
  tot =0.0

  DO K=1,KXM1
    DO J=1,JXM1
      XIXC=.5*(XIX(J,K)+XIX(J+1,K))
      XIYC=.5*(XIY(J,K)+XIY(J+1,K))
      ETXC=.5*(ETX(J,K)+ETX(J,K+1))
      ETYC=.5*(ETY(J,K)+ETY(J,K+1))

      UF(J,K)=Q(2,J,K)/Q(1,J,K)
      VF(J,K)=Q(3,J,K)/Q(1,J,K)
      PF(J,K)=(GM-1.)*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))

      SONF(J,K)=DSQRT(DABS(GM*PF(J,K)/Q(1,J,K)))

      DT(J,K)=CFL/((DABS(XIXC*UF(J,K)+XIYC*VF(J,K))+&
             DABS(ETXC*UF(J,K)+ETYC*VF(J,K))+SONF(J,K)*&
             DSQRT(DABS(XIXC**2+XIYC**2+ETXC**2+ETYC**2)))*RJ(J,K))
      DTMIN=DMIN1(DTMIN,DT(J,K))
    enddo
  enddo

  DO  K=1,KXM1
    DO  J=1,JXM1
      DT(j,k)=DTMIN
    enddo
  enddo

  TOT=TOT+DTMIN

RETURN

END

!=======================================================================
! Subroutine EXPLICIT
! CTOP
!=======================================================================

SUBROUTINE EXPLICIT(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: W
  REAL*8 :: q0,GM

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qw/w(6,-1:1201,-1:251)

  q0=0.5196d10
  GM=1.4

  do k=-1,kmp1
    do j=-1,jmp1
      w(1,j,k)=Q(1,j,k)
      w(2,j,k)=Q(2,J,K)/Q(1,J,K)
      w(3,j,k)=Q(3,J,K)/Q(1,J,K)
      w(4,j,k)=(GM-1.)*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
      w(5,j,k)=Q(5,J,K)/Q(1,J,K)
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine control z
! control z between 0 and 1
!=======================================================================

SUBROUTINE controlz(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: ZFF,DFF,UFF,VFF,PFF
  REAL*8 :: q0,r,GM

  INTEGER :: JMAX,KMAX,JMP1,JXM1,KMP1,KXM1
  INTEGER :: K,J

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  q0=0.5196d10
  r=1.4
  GM=1.4

  DO  K=-1,KMP1  
    DO  J=-1,JMP1
      ZFF=Q(5,J,K)/Q(1,J,K)
      ZFF=min(1.0,ZFF)  
      ZFF=max(0.0,ZFF)
      DFF=Q(1,J,K)
      UFF=Q(2,J,K)/Q(1,J,K)
      VFF=Q(3,J,K)/Q(1,J,K)
      PFF=(GM-1.)*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
      Q(4,J,K)=PFF/(r-1)+0.5*DFF*(UFF**2+VFF**2)+q0*DFF*ZFF
      Q(5,J,K)=DFF*ZFF
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine LOWERWALL
! LOWER WALL
!=======================================================================

SUBROUTINE LOWERWALL(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  real*8 :: TERM
  REAL*8 :: XIX,XIY,ETX,ETY
  REAL*8 :: BKX,BKY,RKX,RKY,RSXY,DB,DR,PB,PR,UB,VB,VR,ZF,UR
  REAL*8 :: UF,VF,PF,SONF
  REAL*8 :: q0,gm

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: J,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/ST/UF(0:1201,0:251),VF(0:1201,0:251),PF(0:1201,0:251),SONF(0:1201,0:251)

  q0=0.5196d10
  GM=1.4

  DO J=1,jxm1  
    RKX=ETX(J,2)
    RKY=ETY(J,2)
    RSXY=1./DSQRT(DABS(RKX*RKX+RKY*RKY))
    BKX=RKX*RSXY
    BKY=RKY*RSXY

    DR=Q(1,J,2)
    PR=PF(J,2)
    UR=UF(J,2)
    VR=VF(J,2)
    ZF=Q(5,J,2)/DR

    TERM=BKX*UR+BKY*VR

    PB=PR  
    DB=DR 
    UB=UR-BKX*TERM
    VB=VR-BKY*TERM

    Q(1,J,1)=DB
    Q(2,J,1)=DB*UB
    Q(3,J,1)=DB*VB
    Q(4,J,1)=PB*(1./(GM-1.))+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
    Q(5,J,1)=DB*ZF
  enddo

  DO  N=1,5 
    DO  J=1,jxm1
      Q(N,J,0)=Q(N,J,1)
      Q(N,J,-1)=Q(N,J,0)
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine UPWALL
! UP WALL
!=======================================================================

SUBROUTINE UPWALL(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: TERM
  REAL*8 :: XIX,XIY,ETX,ETY
  REAL*8 :: BKX,BKY,RKX,RKY,RSXY,DB,DR,PB,PR,UB,VB,VR,ZF,UR
  REAL*8 :: UF,VF,PF,SONF
  REAL*8 :: q0,GM

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: J,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/ST/UF(0:1201,0:251),VF(0:1201,0:251),PF(0:1201,0:251),SONF(0:1201,0:251)

  q0=0.5196d10
  GM=1.4

  DO  J=1,jxm1     
    RKX=ETX(J,KXM1-1)
    RKY=ETY(J,KXM1-1)
    RSXY=1./DSQRT(RKX*RKX+RKY*RKY)
    BKX=RKX*RSXY
    BKY=RKY*RSXY

    DR=Q(1,J,KXM1-1)
    PR=PF(J,KXM1-1)
    UR=UF(J,KXM1-1)
    VR=VF(J,KXM1-1)
    ZF=Q(5,J,KXM1-1)/DR
    TERM=BKX*UR+BKY*VR

    PB=PR
    DB=DR
    UB=UR-BKX*TERM
    VB=VR-BKY*TERM
    Q(1,J,kxm1)=DB
    Q(2,J,kxm1)=DB*UB
    Q(3,J,kxm1)=DB*VB
    Q(4,J,kxm1)=PB*(1./(GM-1.))+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
    Q(5,J,kxm1)=DB*ZF
  enddo

  DO  N=1,5 
    DO  J=1,jxm1
      Q(N,J,kmax)=Q(N,J,kxm1)
      Q(N,J,kmp1)=Q(N,J,kmax)
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine 
! MUSCL
!=======================================================================

SUBROUTINE MUSCL(imax)

  IMPLICIT NONE

  REAL*8 :: CK,HI,EPS,DTP,DTM,SS
  REAL*8 :: uuu1
  REAL*8 :: vpr,vpl

  INTEGER :: IMAX,K,J,JP1

  COMMON/U/uuu1(-1:2001,5)

  COMMON/MUS/vpr(-1:2001,5),vpl(-1:2001,5)

  CK=1./3.
  hi=1.
  eps=1e-14

  DO k=1,5
    DO J=1,imax-1
      JP1=J+1

      DTP=uuu1(JP1,k)-uuu1(J,k)
      DTM=uuu1(J,k)-uuu1(J-1,k)
      SS=(2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)

      vpl(j,k)=uuu1(J,k)+hi*.25*SS*((1.-CK*SS)*DTM+(1.+CK*SS)*DTP)

      DTM=DTP
      DTP=uuu1(JP1+1,k)-uuu1(JP1,k)
      SS=(2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)

      vpr(j,k)=uuu1(JP1,k)-hi*.25*SS*((1.-CK*SS)*DTP+(1.+CK*SS)*DTM)

    enddo
  enddo  

RETURN

END

!=======================================================================
! Subroutine 
! THINCEM
!=======================================================================

SUBROUTINE THINCEM(beta,imax)

  IMPLICIT NONE

  REAL*8 :: epsvof,CFL,cbvof,tbvof
  REAL*8 :: avofm,avof,svof
  REAL*8 :: BB,CC,tbvofp,tvofp,tvofm
  REAL*8 :: uuu1
  REAL*8 :: vr,vl
  REAL*8 :: beta

  INTEGER :: IMAX,K,J,JP,JM

  COMMON/U/uuu1(-1:2001,5)

  COMMON/THI/vr(-1:2001,5) ,vl(-1:2001,5)

  epsvof = 1e-12
  CFL=1e-6

  cbvof = dcosh(beta)
  tbvof = dtanh(beta)

  DO k=1,5
    DO j=1,imax-1

      jp=j+1
      jm=j-1

      jp = min(jp,imax)
      jm = max(jm,1)

      avofm = min(uuu1(jp,k),uuu1(jm,k))
      avof =  max(uuu1(jp,k),uuu1(jm,k))-avofm
      svof = -1.0

      if(avofm.eq.uuu1(jm,k)) svof = 1.0
        BB = exp(svof*beta*(2.0*(uuu1(j,k)-avofm+epsvof)/(avof+epsvof)-1.0))

        CC = (BB/cbvof - 1.0)/tbvof

        tbvofp = (CC+tbvof)/(1.0+CC*tbvof)
        tvofp = log( cosh(beta*cfl)-sinh(beta*cfl)*tbvofp)/(beta*(CFL+1e-12))
        tvofm = log( cosh(beta*cfl)+sinh(beta*cfl)*CC)/(beta*(CFL+1e-12))

        vl(j,k) = avofm + 0.5*avof*(1. - svof*tvofp)
        vr(j-1,k) = avofm + 0.5*avof*(1. + svof*tvofm)

      if(sign( 1d0,(uuu1(jp,k)-uuu1(j,k))*(uuu1(j,k)-uuu1(jm,k))).eq.-1d0) then
        vr(j-1,k) = uuu1(j,k)
        vl(j,k) = uuu1(j,k)
      endif

    ENDDO
  ENDDO

RETURN

END

!=======================================================================
! Subroutine 
! ORI
!=======================================================================

SUBROUTINE ORI(alpha,IMAX)

  IMPLICIT NONE

  REAL*8 :: EPS
  REAL*8 :: DTP,DTM,DIP0,DDP
  REAL*8 :: vpr,vpl
  REAL*8 :: vr,vl
  REAL*8 :: ur,ul
  REAL*8 :: uuu1
  REAL*8 :: alpha

  INTEGER :: JM1,JP1,I,K,IMAX

  COMMON/U/uuu1(-1:2001,5)

  COMMON/MUS/vpr(-1:2001,5),vpl(-1:2001,5)

  COMMON/THI/vr(-1:2001,5) ,vl(-1:2001,5)

  COMMON/ATMOUT/ur(-1:2001,5),ul(-1:2001,5)

  eps =1e-12

  DO  I=1,IMAX-1

    Jm1  = i-1
    Jp1  = i+1
    DTP  = ( uuu1(Jm1,4) - 2*uuu1(i,4) + uuu1(Jp1,4) )
    DTM  = ( uuu1(Jm1,4) + 2*uuu1(i,4) + uuu1(Jp1,4) )
    dip0 = alpha * abs((eps+dtp)/(dtm+eps))
    ddp  = max((1d0 - dip0),0.)

    DO  k=1,5
     ur(i,k) = ddp * vpr(i,k)  + ( 1.- ddp ) * vr(i,k)
     ul(i,k) = ddp * vpl(i,k)  + ( 1.- ddp ) * vl(i,k) 
    enddo

  enddo

RETURN

END

!=======================================================================
! Subroutine 
! ATM
!=======================================================================

SUBROUTINE ATM(alpha,IMAX)

  IMPLICIT NONE

  REAL*8 :: EPS
  REAL*8 :: DTP,DTM,DIP0,DDP
  REAL*8 :: vpr,vpl
  REAL*8 :: vr,vl
  REAL*8 :: ur,ul
  REAL*8 :: uuu1
  REAL*8 :: alpha

  INTEGER :: JM1,JP1,I,K,IMAX

  COMMON/U/uuu1(-1:2001,5)

  COMMON/MUS/vpr(-1:2001,5),vpl(-1:2001,5)

  COMMON/THI/vr(-1:2001,5) ,vl(-1:2001,5)

  COMMON/ATMOUT/ur(-1:2001,5),ul(-1:2001,5)

  eps =1e-12

  DO  k=1,5
    DO  I=1,IMAX-1

      Jm1  = i-1
      Jp1  = i+1
      DTP  = ( uuu1(Jm1,k) - 2*uuu1(i,k) + uuu1(Jp1,k) )
      DTM  = ( uuu1(Jm1,k) + 2*uuu1(i,k) + uuu1(Jp1,k) )
      dip0 = alpha * abs((eps+dtp)/(dtm+eps))
      ddp  = max((1d0 - dip0),0.)

      ur(i,k) = ddp * vpr(i,k)  + ( 1.- ddp ) * vr(i,k)
      ul(i,k) = ddp * vpl(i,k)  + ( 1.- ddp ) * vl(i,k) 

    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine 
! JCP
!=======================================================================

SUBROUTINE JCP(IMAX)

  IMPLICIT NONE

  REAL*8 :: epsvof
  REAL*8 :: vpr,vpl
  REAL*8 :: vr,vl
  REAL*8 :: ur,ul
  REAL*8 :: uuu1
  REAL*8 :: r,l,qmp1,qmp1e,dvnp,qmm1,qmm1e,dvnm,zeta

  INTEGER :: JM1,JP1,I,K,IMAX

  COMMON/U/uuu1(-1:2001,5)

  COMMON/MUS/vpr(-1:2001,5),vpl(-1:2001,5)

  COMMON/THI/vr(-1:2001,5) ,vl(-1:2001,5)

  COMMON/ATMOUT/ur(-1:2001,5),ul(-1:2001,5)

  epsvof = 1e-10

  do k=1,5
    do i=1,IMAX-1
      Jm1=i-1
      Jp1=i+1

      r = vpl(i,k)-vpr(Jm1,k)+epsvof
      qmp1 = uuu1(Jp1,k)-uuu1(i,k)
      qmp1e = uuu1(Jp1,k)-uuu1(i,k)+epsvof
      dvnp = r/qmp1e
      
      l = vpl(i,k)-vpr(Jm1,k)+epsvof
      qmm1 = uuu1(i,k)-uuu1(Jm1,k)
      qmm1e = uuu1(i,k)-uuu1(Jm1,k)+epsvof
      dvnm = l/qmm1e

      zeta = 1.0-min(dvnp,dvnm)

      if ((qmp1*qmm1).GT.0.0) then
        ur(i,k) = ( 1.- zeta ) * vpr(i,k) + zeta * vr(i,k)
        ul(i,k) = ( 1.- zeta ) * vpl(i,k) + zeta * vl(i,k)       
      else
        ur(i,k) = vpr(i,k)
        ul(i,k) = vpl(i,k)
      endif
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine 
! MUS
!=======================================================================

SUBROUTINE PMUS(IMAX)

  IMPLICIT NONE

  REAL*8 :: vpr,vpl
  REAL*8 :: ur,ul

  INTEGER :: I,K,IMAX

  COMMON/MUS/vpr(-1:2001,5),vpl(-1:2001,5)

  COMMON/ATMOUT/ur(-1:2001,5),ul(-1:2001,5)

  do k=1,5
    do i=1,IMAX-1

      ur(i,k) = vpr(i,k)
      ul(i,k) = vpl(i,k)     

    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine 
! THI
!=======================================================================

SUBROUTINE PTHI(IMAX)

  IMPLICIT NONE

  REAL*8 :: vr,vl
  REAL*8 :: ur,ul

  INTEGER :: I,K,IMAX

  COMMON/THI/vr(-1:2001,5) ,vl(-1:2001,5)

  COMMON/ATMOUT/ur(-1:2001,5),ul(-1:2001,5)

  do k=1,5
    do i=1,IMAX-1

      ur(i,k) = vr(i,k)
      ul(i,k) = vl(i,k)       

    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine sour1
! source ONE-S(k/2)
!=======================================================================

SUBROUTINE sour1(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: P,HH,temperature,source,source2
  REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS
  REAL*8 :: q0,temperature0,ddd,GM
  
  REAL*8 :: TOT,DTMIN

  INTEGER :: JMAX,KMAX,JMP1,JXM1,KMP1,KXM1
  INTEGER :: J,K,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/P/p(-1:1201,-1:251),hh(-1:1201,-1:251),&
           temperature(-1:1201,-1:251),&
           source(-1:1201,-1:251,5),source2(-1:1201,-1:251,5)

  COMMON/Qold/Qold1ds(-1:1201,-1:251,5),Qold1s(-1:1201,-1:251,5),&
              Qold1d(-1:1201,-1:251,5),Qold1(-1:1201,-1:251,5),&
              Qold2d(-1:1201,-1:251,5),Qold2ss(-1:1201,-1:251,5),&
              Qolddss(-1:1201,-1:251,5)

  COMMON/TIME/TOT,DTMIN

  q0=0.5196d10
  temperature0=0.1155d10
  ddd=1.0/0.5825d10
  GM=1.4

  DO  K=-1,KMP1
    DO  J=-1,JMP1
    p(j,k)=(GM-1.)*(Q(4,J,K)-0.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
    enddo
  enddo

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      temperature(j,k)=p(j,k)/Q(1,J,K)
    enddo
  enddo

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      if ((temperature(j,k)-temperature0)>0.0)then
        hh(j,k)=1.0
      else
        hh(j,k)=0.0  
      endif
    enddo
  enddo

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      source(j,k,1)=0.0
      source(j,k,2)=0.0
      source(j,k,3)=0.0
      source(j,k,4)=0.0
      source(j,k,5)=-(10d-1/ddd)*Q(5,J,K)*hh(j,k)
    enddo
  enddo

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      source2(j,k,1)=0.0
      source2(j,k,2)=0.0
      source2(j,k,3)=0.0
      source2(j,k,4)=0.0
      source2(j,k,5)=-10d-1*hh(j,k)
    enddo
  enddo

  do  n=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Qold1ds(j,k,n)=5d-1*DTMIN*source(j,k,n)/(10d-1-5d-1*5d-1*DTMIN*source2(j,k,n))
        Qold1s(j,k,n)=Q(n,J,K)+Qold1ds(j,k,n)
      enddo
    enddo
  enddo

  do  n=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Q(n,J,k)=Qold1s(j,k,n)
      enddo
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine sour2
! source TWO-S(k/2)
!=======================================================================

SUBROUTINE sour2(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: P,HH,temperature,source,source2
  REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS
  real*8 :: TOT,DTMIN
  REAL*8 :: q0,temperature0,ddd,gm

  INTEGER :: JMAX,KMAX,JMP1,JXM1,KMP1,KXM1
  INTEGER :: J,K,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/P/p(-1:1201,-1:251),hh(-1:1201,-1:251),&
           temperature(-1:1201,-1:251),&
           source(-1:1201,-1:251,5),source2(-1:1201,-1:251,5)

  COMMON/Qold/Qold1ds(-1:1201,-1:251,5),Qold1s(-1:1201,-1:251,5),&
              Qold1d(-1:1201,-1:251,5),Qold1(-1:1201,-1:251,5),&
              Qold2d(-1:1201,-1:251,5),Qold2ss(-1:1201,-1:251,5),&
              Qolddss(-1:1201,-1:251,5)


  COMMON/TIME/TOT,DTMIN

  q0=0.5196d10
  temperature0=0.1155d10
  ddd=1.0/0.5825d10
  GM=1.4

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      p(J,K)=(GM-1.)*(Qold2ss(J,K,4)-5d-1*(Qold2ss(J,K,2)**2+Qold2ss(J,K,3)**2)/Qold2ss(J,K,1)-q0*Qold2ss(J,K,5))
    enddo
  ENDDO

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      temperature(J,K)=p(J,K)/Qold2ss(J,K,1)
    enddo
  ENDDO

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      if ((temperature(J,K)-temperature0)>0.0)then
        hh(J,K)=10d-1
      else
        hh(J,K)=0.0  
      endif
    enddo
  ENDDO

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      source(J,K,1)=0.0
      source(J,K,2)=0.0
      source(J,K,3)=0.0
      source(J,K,4)=0.0
      source(J,K,5)=-(10d-1/ddd)*Qold2ss(J,K,5)*hh(J,K)
    enddo
  ENDDO

  DO  K=-1,KMP1
    DO  J=-1,JMP1
      source2(J,K,1)=0.0
      source2(J,K,2)=0.0
      source2(J,K,3)=0.0
      source2(J,K,4)=0.0
      source2(J,K,5)=-10d-1*hh(J,K)
    enddo
  ENDDO

  do N=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Qolddss(J,K,N)=(5d-1*DTMIN*source(J,K,N))/(10d-1-5d-1*5d-1*DTMIN*source2(J,K,N))
        Q(N,J,K)=Qold2ss(J,K,N)+Qolddss(J,K,N)
      enddo
    enddo
  ENDDO

RETURN

END

!=======================================================================
! Subroutine 
! Qoo1
!=======================================================================

SUBROUTINE Qoo1(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  real*8 :: S
  REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qold/Qold1ds(-1:1201,-1:251,5),Qold1s(-1:1201,-1:251,5),&
              Qold1d(-1:1201,-1:251,5),Qold1(-1:1201,-1:251,5),&
              Qold2d(-1:1201,-1:251,5),Qold2ss(-1:1201,-1:251,5),&
              Qolddss(-1:1201,-1:251,5)

  COMMON/SD/S(6,0:1201,0:251)

  do n=1,5
    DO  K=0,KMP1
      DO  J=0,JMP1
        Qold1d(j,k,n)=S(N,J,K)
        Qold1(j,k,n)=Qold1s(j,k,n)+Qold1d(j,k,n)
      enddo
    enddo
  enddo

  do  n=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Q(n,J,k)=Qold1(j,k,n)
      enddo
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine 
! Qoo2
!=======================================================================

SUBROUTINE Qoo2(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: S
  REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qold/Qold1ds(-1:1201,-1:251,5),Qold1s(-1:1201,-1:251,5),&
              Qold1d(-1:1201,-1:251,5),Qold1(-1:1201,-1:251,5),&
              Qold2d(-1:1201,-1:251,5),Qold2ss(-1:1201,-1:251,5),&
              Qolddss(-1:1201,-1:251,5)

  COMMON/SD/S(6,0:1201,0:251)

  DO N=1,5
    DO K=0,KMP1
      DO J=0,JMP1
        Qold2d(J,K,N)=S(N,J,K)
        Qold2ss(J,K,N)=Qold1s(J,K,N)+0.5*(Qold1d(J,K,N)+Qold2d(J,K,N))
      enddo
    enddo
  ENDDO

  do  n=1,5
    DO  K=-1,KMP1
      DO  J=-1,JMP1
        Q(n,J,k)=Qold2ss(j,k,n)
      enddo
    enddo
  enddo

RETURN

END

!=======================================================================
! Subroutine ghostpoint
! ghost point REFERENCE
!=======================================================================

SUBROUTINE ghostpoint(Q)

  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  DO  N=1,5    
    DO  K=2,KXM1-1
      Q(N,2,K)=Q(N,3,K)
      Q(N,1,K)=Q(N,2,K)
      Q(N,0,K)=Q(N,1,K)
      Q(N,-1,K)=Q(N,0,K)
      Q(N,JXM1-1,K)=Q(N,JXM1-1-1,K)  
      Q(N,JXM1,K)=Q(N,JXM1-1,K)  
      Q(N,JMAX,K)=Q(N,JXM1,K)  
      Q(N,Jmp1,K)=Q(N,jmax,K)
    enddo
  enddo

RETURN

END


!=======================================================================
! Subroutine E TERM-1
! ETERM1
!=======================================================================

SUBROUTINE ETERM1(K)
  
  IMPLICIT NONE

  real*8 :: W,S,DT,RJ,F
  real*8 :: RKX,RKY
  real*8 :: XIX,XIY,ETX,ETY
  real*8 :: uuu1
  real*8 :: QBF(6,2)

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qw/w(6,-1:1201,-1:251)

  COMMON/SD/S(6,0:1201,0:251)

  COMMON/Dtt/DT(0:1201,0:251)

  COMMON/JACO/RJ(-1:1201,-1:251)

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/AUSMF/F(6,2000)

  COMMON/U/uuu1(-1:2001,5)

  RKX=XIX(1,K)
  RKY=XIY(1,K)
  DO N=1,5
    QBF(N,1)=w(N,0,K)
    QBF(N,2)=w(N,1,K)
  enddo

  CALL FVSW(1,RKX,RKY,QBF)

  DO J=1,JXM1+1
    DO N=1,5
      uuu1(J,N)=W(N,J,K)
    ENDDO
  ENDDO

RETURN

END

!=======================================================================
! Subroutine E TERM-2
! ETERM2
!=======================================================================

SUBROUTINE ETERM2(K)
  
  IMPLICIT NONE

  real*8 :: W,S,DT,RJ,F
  real*8 :: RKX,RKY,HE1DT
  real*8 :: XIX,XIY,ETX,ETY
  real*8 :: ur,ul
  real*8 :: QBF(6,2)

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N,JP1

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qw/w(6,-1:1201,-1:251)

  COMMON/SD/S(6,0:1201,0:251)

  COMMON/Dtt/DT(0:1201,0:251)

  COMMON/JACO/RJ(-1:1201,-1:251)

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/AUSMF/F(6,2000)

  COMMON/ATMOUT/ur(-1:2001,5),ul(-1:2001,5)

  DO J=1,JXM1-1
    JP1=J+1

    RKX=XIX(JP1,K)
    RKY=XIY(JP1,K)

    DO N=1,5
      QBF(N,1)=ul(J,N)
      QBF(N,2)=ur(J,N)
    ENDDO

    CALL FVSW(JP1,RKX,RKY,QBF)

  enddo

  RKX=XIX(JMAX,K)
  RKY=XIY(JMAX,K)

  DO N=1,5
    QBF(N,1)=W(N,JXM1,K)
    QBF(N,2)=W(N,JMAX,K)
  enddo

  CALL FVSW(JMAX,RKX,RKY,QBF)

!-----------------RHS=-HE1*DT*(E(J+.5)-E(J-.5))+HE2*DQ(N-1)-------------

  DO J=1,JXM1
    HE1DT=DT(J,K)*RJ(J,K)
    DO N=1,5
      S(N,J,K)=-HE1DT*(F(N,J+1)-F(N,J))
    ENDDO
  ENDDO

RETURN

END

!=======================================================================
! Subroutine F TERM-1
! FTERM1
!=======================================================================

SUBROUTINE FTERM1(J)
  
  IMPLICIT NONE

  real*8 :: W,S,DT,RJ,F
  real*8 :: RKX,RKY
  real*8 :: XIX,XIY,ETX,ETY
  real*8 :: uuu1
  real*8 :: QBF(6,2)

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qw/w(6,-1:1201,-1:251)

  COMMON/SD/S(6,0:1201,0:251)

  COMMON/Dtt/DT(0:1201,0:251)

  COMMON/JACO/RJ(-1:1201,-1:251)

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/AUSMF/F(6,2000)

  COMMON/U/uuu1(-1:2001,5)

  RKX=ETX(J,1)
  RKY=ETY(J,1)
  DO N=1,5
    QBF(N,1)=W(N,J,0)
    QBF(N,2)=W(N,J,1)
  ENDDO

  CALL FVSW(1,RKX,RKY,QBF)

  DO  K=1,KXM1+1
    DO N=1,5
      uuu1(K,N)=W(N,J,K)
    ENDDO
  ENDDO

RETURN

END

!=======================================================================
! Subroutine F TERM-2
! FTERM2
!=======================================================================

SUBROUTINE FTERM2(J)
  
  IMPLICIT NONE

  real*8 :: W,S,DT,RJ,F
  real*8 :: RKX,RKY,HE1DT
  real*8 :: XIX,XIY,ETX,ETY
  real*8 :: ur,ul
  real*8 :: QBF(6,2)

  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
  INTEGER :: K,J,N,KP1

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/Qw/w(6,-1:1201,-1:251)

  COMMON/SD/S(6,0:1201,0:251)

  COMMON/Dtt/DT(0:1201,0:251)

  COMMON/JACO/RJ(-1:1201,-1:251)

  COMMON/XEXY/XIX(0:1201,0:251),XIY(0:1201,0:251),ETX(0:1201,0:251),ETY(0:1201,0:251)

  COMMON/AUSMF/F(6,2000)

  COMMON/ATMOUT/ur(-1:2001,5) ,ul(-1:2001,5)

  DO K=1,KXM1-1
    KP1=K+1

    RKX=ETX(J,KP1)
    RKY=ETY(J,KP1)

    DO N=1,5
      QBF(N,1)=ul(K,N)
      QBF(N,2)=ur(K,N)
    ENDDO
    
    CALL FVSW(KP1,RKX,RKY,QBF)

  ENDDO

  RKX=ETX(J,KMAX)
  RKY=ETY(J,KMAX)

  DO N=1,5
    QBF(N,1)=W(N,J,KXM1)
    QBF(N,2)=W(N,J,KMAX)
  ENDDO

  CALL FVSW(KMAX,RKX,RKY,QBF)

!-----------------RHS=RHS-HE1*DT*((F(K+.5)-F(K-.5))---------------------

  DO K=1,KXM1
    HE1DT=DT(J,K)*RJ(J,K)
    DO N=1,5
      KP1=K+1
      S(N,J,K)=S(N,J,K)-HE1DT*(F(N,KP1)-F(N,K))
    ENDDO
  ENDDO

RETURN

END

!=======================================================================
! Subroutine AFVS
! AUSMDV FLUX DIFFERENCE SPLITTING
!=======================================================================

SUBROUTINE FVSW(I,xiw,eta,QA)

  IMPLICIT NONE

  REAL*8 :: Q0,A1,A2,A12,F
  REAL*8 :: TXY,XIW,ETA
  REAL*8 :: D1,D2,QA,U1,U2,V1,V2,P1,P2,Z1,Z2,C1,C2
  REAL*8 :: PD1,PD2,UU1,UU2,CU1,CU2
  REAL*8 :: ALPHAL,ALPHAR
  REAL*8 :: CML,CMR,FLR
  REAL*8 :: AUSMDU,AUSMDV,AUSMVU,AUSMVV,DP,COEF,ST

  INTEGER :: N,I,K

  COMMON/AUSMF/F(6,2000)

  DIMENSION QA(6,2),fLR(2,6)

  q0=0.5196d10
  a1=0.
  a2=0.

  DO N=1,5
    F(N,I)=0.
  enddo

  txy=dsqrt(xiw*xiw+eta*eta)

  d1=qa(1,1)
  u1=qa(2,1)
  v1=qa(3,1)
  p1=qa(4,1)
  z1=qa(5,1)

  c1=dsqrt(dabs(p1*1.4/qa(1,1)))
  pd1=p1/qa(1,1)

  uu1=(xiw*u1+eta*v1)

!----------------2------------------------------------------------------

  d2=qa(1,2)
  u2=qa(2,2)
  v2=qa(3,2)
  p2=qa(4,2)
  z2=qa(5,2)

  c2=dsqrt(dabs(p2*1.4/qa(1,2)))
  pd2=p2/qa(1,2)

  uu2=(xiw*u2+eta*v2)

!----------------interface----------------------------------------------

  a12=(c1+c2)/2.
  cu1=uu1/txy/a12
  cu2=uu2/txy/a12
  alphaL=2*pd1/(pd1+pd2)
  alphaR=2*pd2/(pd1+pd2)

!----------------for particle based on density--------------------------

  IF (ABS(cu1).Le.1.0) THEN
    CML=alphaL*((cu1+1.)**2/4-(cu1+abs(cu1))/2.)+(cu1+abs(cu1))/2.
    a1=a1+CML*d1
    a2=a2+p1*(cu1+1.)**2*(2-cu1)/4.
  ELSE
    CML=(cu1+abs(cu1))/2.
    a1=a1+CML*d1
    a2=a2+p1*(cu1+abs(cu1))/(2*cu1)
  endif

!----------------2------------------------------------------------------

  IF (ABS(cu2).Le.1.0) THEN
    CMR=alphaR*(-(cu2-1.)**2/4-(cu2-abs(cu2))/2.)+(cu2-abs(cu2))/2.
    a1=a1+CMR*d2
    a2=a2+p2*(cu2-1.)**2*(2+cu2)/4.

  ELSE
    CMR=(cu2-abs(cu2))/2.
    a1=a1+CMR*d2
    a2=a2+p2*(cu2-abs(cu2))/(2*cu2)
   endif


  fLR(1,1)=1.
  fLR(1,2)=u1
  fLR(1,3)=v1
  fLR(1,4)=(p1*3.5+.5*d1*(u1*u1+v1*v1) +q0*z1*d1)/d1
  fLR(1,5)=z1

  fLR(2,1)=1.
  fLR(2,2)=u2
  fLR(2,3)=v2
  fLR(2,4)=(p2*3.5+.5*d2*(u2*u2+v2*v2)+q0*z2*d2)/d2
  fLR(2,5)=z2

  do  k=1,1
    f(k,I)=.5*txy*a1*a12*(fLR(1,k)+fLR(2,k))-abs(.5*txy*a1)*a12*(fLR(2,k)-fLR(1,k))
  enddo

  do  k=5,5
    f(k,I)=.5*txy*a1*a12*(fLR(1,k)+fLR(2,k))-abs(.5*txy*a1)*a12*(fLR(2,k)-fLR(1,k))
  enddo

  ausmdU=.5*txy*a1*a12*(fLR(1,2)+fLR(2,2))-abs(.5*txy*a1)*a12*(fLR(2,2)-fLR(1,2))

  ausmdV=.5*txy*a1*a12*(fLR(1,3)+fLR(2,3))-abs(.5*txy*a1)*a12*(fLR(2,3)-fLR(1,3))

  ausmvU=(CML*d1*u1+CMR*d2*u2)*a12*txy

  ausmvV=(CML*d1*v1+CMR*d2*v2)*a12*txy

  dp = p1-p2
  coef = 10.*abs(dp)/(p1+p2)
  st = .5*min(1.d0,coef)

  f(2,i)=(.5+st)*ausmvU+(.5-st)*ausmdU
  f(3,i)=(.5+st)*ausmvv+(.5-st)*ausmdv

!-----------------------------------------------------------------------

  do  k=4,4
    f(k,I)=.5*txy*a1*a12*(fLR(1,k)+fLR(2,k))-abs(.5*txy*a1)*a12*(fLR(2,k)-fLR(1,k))
  enddo

  f(2,I)=f(2,I)+xiw*a2
  f(3,I)=f(3,I)+eta*a2

RETURN

END

!=======================================================================
! Subroutine PRINT/PLOT FINAL DATA
! SavePlt
!=======================================================================

SUBROUTINE SavePlt(Q,IT,ISAVE)
  
  IMPLICIT NONE

  REAL*8 :: Q(6,-1:1201,-1:251)

  REAL*8 :: S,X,Y
  REAL*8 :: RESID,DOM
  REAL*8 :: UF,VF,PF,SONF
  real*8 :: TOT,DTMIN

  INTEGER :: J,K,IT,ISAVE,itt
  INTEGER :: JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  CHARACTER*80 :: filename

  COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

  COMMON/SD/S(6,0:1201,0:251)

  COMMON/XY/X(0:1201,0:251),Y(0:1201,0:251)

  COMMON/ST/UF(0:1201,0:251),VF(0:1201,0:251),PF(0:1201,0:251),SONF(0:1201,0:251)

  COMMON/TIME/TOT,DTMIN

  RESID=0.

  DO K=1,KXM1
    DO J=1,JXM1
      RESID=RESID+S(1,J,K)*S(1,J,K)
    enddo
  enddo

  RESID=DSQRT(DABS(RESID/(JXM1*KXM1)))

!=======================================================================
!  IF (MOD(IT,ISAVE).EQ.0.or.it.eq.1) THEN
!    WRITE (*,*) IT,RESID,tot,dtmin
!  ENDIF
!=======================================================================

  IF (MOD(IT,ISAVE).EQ.0.or.it.eq.1) THEN
    itt = it/ISAVE + 7000
!    WRITE (*,*) IT,RESID,q(2,4,2),q(1,4,2)
    write (filename,'("chan",I5,".plt")')itt
    open (unit=itt,file=filename,status='unknown')
    WRITE (itt,*)'Title="FLATE PLANE BY k-e MODEL SIMULATION"'
    WRITE (itt,*)'Variables = "X","Y","U","V","d","p","M","z" ' 
    WRITE (itt,*)'zone I=',JXm1,',J=',Kxm1,',f=point'

    DO  K=1,Kxm1
      DO  J=1,jxm1
        dom = dsqrt( UF(J,K)**2 + VF(J,K)**2)
        WRITE (itt,*) X(J,K),Y(J,K),UF(J,K),VF(J,K),Q(1,J,K),PF(j,k),dom,Q(5,J,K)/Q(1,J,K)
      enddo
    enddo
    close(itt)

!=======================================================================
!    timeout=tot-(3d-8)
!
!    if( timeout > 0.0  )then
!      write(*,*) 'time out 1'
!    endif
!
!    timeout=tot-(92d-9)
!
!    if( timeout > 0.0  )then
!      write(*,*) 'time out 2'
!    endif
!
!    timeout=tot-(17d-8)
!
!    if( timeout > 0.0  )then
!      write(*,*) 'time out 3'
!    endif
!=======================================================================

  ENDIF

RETURN

END
