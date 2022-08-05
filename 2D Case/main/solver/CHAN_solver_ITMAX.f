!=======================================================================
! 2-D UPWINDING
! VAN LEER FOR THE THIN LAYER GAS FLOW
! WITH AUSMD AND RNG K-E TURBULENCE MODEL
! PREDICITE PARTICULATE TRAJECTORY BY LAGRANGIAN METHOD
!=======================================================================

      subroutine solver(alpha,beta,ISAVE,ITMAX,imethod,ngrd,Q)

      IMPLICIT NONE

      real*8, intent(in)  :: alpha,beta
      real*8, intent(out) :: Q(6,-1:3250,-1:650)
      INTEGER, intent(in) :: ISAVE,ITMAX,imethod,ngrd

      real*8 :: X(0:3250,0:650),Y(0:3250,0:650)
      REAL*8 :: RMACH,CFL,EPS,CK
      REAL*8 :: XIXCF,XIYCF,ETXCF,ETYCF
      REAL*8 :: UF,VF,PF,SONF,UINF,VINF
      REAL*8 :: TEMPERATURE0,TEMPERATURE
      REAL*8 :: HALF,ONE,TWO,ZERO,TEN,HH
      REAL*8 :: DDD,R,RGM,RGM1,RGGM1,RG2M1,GM12,GM1R4,GRGM1
      REAL*8 :: CSP,CPG,C13,C43,PRL,PRT,VEO
      REAL*8 :: Z,CASEE,Q1Z,PZ,PZ2,Z2
      REAL*8 :: PINF
      REAL*8 :: VOL,RJ,W
      REAL*8 :: TOT,DTMIN,DTMAX,DT
      REAL*8 :: XIXC,XIYC,ETXC,ETYC
      REAL*8 :: DMIN1
      REAL*8 :: QBF,F,S,GM
      REAL*8 :: RKX,RKY,RSXY,PB,UB,VB,ZF
      REAL*8 :: DB,BKX,BKY,UR,VR,DR,PR
      REAL*8 :: TERM,ZFF,DFF,UFF,VFF,PFF,P
      REAL*8 :: SOURCE,SOURCE2
      REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS,RESID
      REAL*8 :: HE1DT
      REAL*8 :: GM1,DEL,CT0
      REAL*8 :: VISL,VIST
      REAL*8 :: XIX,XIY,ETX,ETY
      REAL*8 :: Q0,UUR,UUL,UUU1
      REAL*8 :: AKP1,AKY,AJP1,AJX,U,Q1Z2,U2,DOM

      INTEGER :: J,K,N,JK,JP1,KP1
      INTEGER :: JMAX,KMAX,JMP1,JXM1,JXM2,KMP1,KXM1,KXM2
      INTEGER :: IT
      INTEGER :: inlet_frequence
      INTEGER :: splitting
      INTEGER :: ITT

      INTEGER :: split,inlet_closed,outlet_closed,subsonic_outflow
      INTEGER :: muscl,hmsthh
      INTEGER :: ghostpoint_open,INLET_SUBSONIC_INFLOW

      CHARACTER*80 :: filename

c-----------------------------------------------------------------------

      PARAMETER(RMACH=340.1,CFL=0.1,EPS=1.E-6,CK=0.33333)

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

      COMMON/GMD1/DEL,UINF,VINF,PINF

      COMMON/SD/S(6,0:3250,0:650)

      COMMON/XEXY/XIX(0:3250,0:650),XIY(0:3250,0:650),

     & ETX(0:3250,0:650),ETY(0:3250,0:650)

      COMMON/JACO/RJ(-1:3250,-1:650)

      COMMON/Qw/w(6,-1:3250,-1:650)

      COMMON/VISC/CT0,C43,C13,PRL,PRT

      COMMON/VISC1/VISL(-1:3250,-1:650),VIST(-1:3250,-1:650)

      DIMENSION DT(0:3250,0:650),QBF(6,2),F(6,3250)
      DIMENSION p(-1:3250,-1:650)
      DIMENSION temperature(-1:3250,-1:650)
      DIMENSION hh(-1:3250,-1:650)
      DIMENSION source(-1:3250,-1:650,5)
      DIMENSION source2(-1:3250,-1:650,5)
      DIMENSION Qold1ds(-1:3250,-1:650,5)
      DIMENSION Qold1s(-1:3250,-1:650,5)
      DIMENSION Qold1d(-1:3250,-1:650,5)
      DIMENSION Qold1(-1:3250,-1:650,5)
      DIMENSION Qold2d(-1:3250,-1:650,5)
      DIMENSION Qold2ss(-1:3250,-1:650,5)
      DIMENSION Qolddss(-1:3250,-1:650,5)
      DIMENSION uuu1(-1:3250,5),uur(-1:3250,5),uul(-1:3250,5)

c-----------------STATEMENT FUNCTION------------------------------------

      XIXCF(J,K)=.5*(XIX(J,K)+XIX(J+1,K))
      XIYCF(J,K)=.5*(XIY(J,K)+XIY(J+1,K))
      ETXCF(J,K)=.5*(ETX(J,K)+ETX(J,K+1))
      ETYCF(J,K)=.5*(ETY(J,K)+ETY(J,K+1))
      UF(J,K)=Q(2,J,K)/Q(1,J,K)
      VF(J,K)=Q(3,J,K)/Q(1,J,K)
      PF(J,K)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-
     #        q0*Q(5,J,K))


      SONF(J,K)=DSQRT(DABS(GM*PF(J,K)/Q(1,J,K)))


c-----------------PARAMETER SETUP---------------------------------------

      q0=0.5196d10
      temperature0=0.1155d10
      ddd=1.0/0.5825d10
      r=1.4
      half=5d-1
      one=10d-1
      two=one+one
      zero=0.0
      ten=1d1
      GM=1.4
      RGM=1./GM
      GM1=GM-1.
      RGM1=1./GM1
      GM12=2.*GM1
      RGGM1=1./(GM*GM1)
      GM1R4=.25*GM1
      RG2M1=1./(GM*GM-1.)
      GRGM1=GM*RGM1
      CSP=1.294
      CPG=CSP*RGM1
      C13=1./3.
      C43=4./3.
      PRT=.9
      PRL=.72
      VEO=7.0

c-----------------switch------------------------------------------------

      casee=3

      muscl=0
      hmsthh=1

      splitting=1

      inlet_closed=0
      inlet_subsonic_inflow=0   
      inlet_frequence=0

      outlet_closed=0
      subsonic_outflow=0

      ghostpoint_open=1

c-----------------READ GRID DATA----------------------------------------

      CALL GRD(X,Y,JMAX,KMAX,ngrd)

c-----------------INITIALIZE DATA---------------------------------------

c-----------------SET UP INDEX------------------------------------------
      
      JMP1=JMAX+1
      JXM1=JMAX-1
      JXM2=JXM1-1
      KMP1=KMAX+1
      KXM1=KMAX-1
      KXM2=KXM1-1
      JK=JXM1*KXM1

c-----------------Q TERM------------------------------------------------

      akp1=-2.0

      DO  K=-1,kMP1 

        akp1=akp1+1.0
        aky=dabs( akp1*x(2,1) -0.0025)

        if( aky .GE. 0.001 )then
          aky=0.004
        elseif(aky .LT. 0.001)then
          aky=0.005-aky
        endif

        ajp1=-2.0

        DO  J=-1,JMP1 
          ajp1=ajp1+1.0
          ajx=ajp1*x(2,1)

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
          endif ! if(ajx .LE. aky)then

        enddo
      enddo

c-----------------------------------------------------------------------

      do n=1,5
        DO  K=-1,KMP1
          DO  J=-1,JMP1
            Qold1ds(J,K,n)=Q(n,J,K)
            Qold1s(J,K,n)=Q(n,J,K)
            Qold1d(J,K,n)=Q(n,J,K)
            Qold1(J,K,n)=Q(n,J,K)
            Qold2d(J,K,n)=Q(n,J,K)
            Qold2ss(J,K,n)=Q(n,J,K)
            Qolddss(J,K,n)=Q(n,J,K)
          enddo
        enddo
      enddo

c-----------------ETX,ETY-----------------------------------------------

      DO K=1,KMAX
        DO J=1,JXM1
          JP1=J+1
          ETX(J,K)=-Y(JP1,K)+Y(J,K) !9-16
          ETY(J,K)=X(JP1,K)-X(J,K)  !9-17
          ETX(0,K)=ETX(1,K)
          ETY(0,K)=ETY(1,K)
          ETX(JMAX,K)=ETX(JXM1,K)
          ETY(JMAX,K)=ETY(JXM1,K)
        enddo
      enddo

c-----------------XIX,XIY-----------------------------------------------     

      DO K=1,KXM1
        DO J=1,JMAX
          KP1=K+1
          XIX(J,K)=Y(J,KP1)-Y(J,K)  !9-14 
          XIY(J,K)=-X(J,KP1)+X(J,K) !9-15
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

c-----------------RJ=1./VOL---------------------------------------------        

      DO K=1,KXM1
        DO J=1,JXM1
          KP1=K+1
          JP1=J+1
          VOL=((X(J,K)-X(JP1,KP1))*(Y(JP1,K)-Y(J,KP1))-
     &        (X(JP1,K)-X(J,KP1))*(Y(J,K)-Y(JP1,KP1)))/2. ! what 2
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
      
c-----------------TIME MARCHING-----------------------------------------   

      tot =0.
      
      DO IT=1,ITMAX

c-----------------SPACE VARYING DT BASED ON VARIABLE CFL.---------------   

        DTMIN=100.
        DTMAX=0.

        DO K=1,KXM1
          DO J=1,JXM1
            XIXC=XIXCF(J,K)
            XIYC=XIYCF(J,K)
            ETXC=ETXCF(J,K)
            ETYC=ETYCF(J,K)

            DT(J,K)=CFL/((DABS(XIXC*UF(J,K)+XIYC*VF(J,K))+
     &            DABS(ETXC*UF(J,K)+ETYC*VF(J,K))+SONF(J,K)*
     &            DSQRT(DABS(XIXC**2+XIYC**2+ETXC**2+ETYC**2)))*RJ(J,K))
            DTMIN=DMIN1(DTMIN,DT(J,K))
          enddo
        enddo

        DO  K=1,KXM1
          DO  J=1,JXM1
            DT(j,k)=DTMIN
          enddo
        enddo

        tot=tot+dtmin

c-----------------FAR FIELD NUMERICAL B.C.------------------------------

        if(splitting .EQ. 1)then

          do split=1,4

            if(split .EQ. 1)then

c-----------------Ssource ONE-S(k/2)------------------------------------   

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                p(j,k)=GM1*(Q(4,J,K)-0.5*(Q(2,J,K)**2+Q(3,J,K)**2)/
     #                 Q(1,J,K)-q0*Q(5,J,K))
                enddo
              enddo

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  temperature(j,k)=p(j,k)/Q(1,J,K)
                enddo
              enddo

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  if ((temperature(j,k)-temperature0)>zero)then
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
                  source(j,k,5)=-(one/ddd)*Q(5,J,K)*hh(j,k)
                enddo
              enddo

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  source2(j,k,1)=0.0
                  source2(j,k,2)=0.0
                  source2(j,k,3)=0.0
                  source2(j,k,4)=0.0
                  source2(j,k,5)=-one*hh(j,k)
                enddo
              enddo

              do  n=1,5
                DO  K=-1,KMP1
                  DO  J=-1,JMP1
                    Qold1ds(j,k,n)=half*DTMIN*source(j,k,n)/ 
     #                             (one-half*half*DTMIN*source2(j,k,n))
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
              
c-----------------end--S(k/2)-------------------------------------------   

            elseif(split .EQ. 2)then

c-----------------MAIN ONE-SF(k)PACE VARYING DT BASED ON VARIABLE CFL.--

c-----------------EXPLICIT-1--------------------------------------------   

              do k=-1,kmp1
                do j=-1,jmp1
                  w(1,j,k)=Q(1,j,k)
                  w(2,j,k)=Q(2,J,K)/Q(1,J,K)
                  w(3,j,k)=Q(3,J,K)/Q(1,J,K)
                  w(4,j,k)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)
     #                    /Q(1,J,K)-q0*Q(5,J,K))
                  w(5,j,k)=Q(5,J,K)/Q(1,J,K)
                enddo
              enddo

c-----------------E TERM------------------------------------------------ 

              DO K=1,KXM1
                RKX=XIX(1,K)
                RKY=XIY(1,K)
                DO N=1,5
                  QBF(N,1)=w(N,0,K)
                  QBF(N,2)=w(N,1,K)
                enddo
                CALL FVSW(1,RKX,RKY,QBF,F) !!!!!!!e
              
                if(hmsthh .eq. 1)then

                  DO J=1,JXM2+2
                    DO N=1,5
                      uuu1(J,N)=W(N,J,K)
                    ENDDO
                  ENDDO

                  CALL hmsth(uuu1,uur,uul,JXM2+2,alpha,beta,imethod)

                  DO J=1,JXM2
                    JP1=J+1

                    RKX=XIX(JP1,K)
                    RKY=XIY(JP1,K)

                    DO N=1,5
                      QBF(N,1)=uul(J,N)
                      QBF(N,2)=uur(J,N)
                    ENDDO
                    CALL FVSW(JP1,RKX,RKY,QBF,F) !!!!!!!e

                  enddo

                endif !	if(hmsthh .eq. 1)then

                RKX=XIX(JMAX,K)
                RKY=XIY(JMAX,K)
                DO N=1,5
                  QBF(N,1)=W(N,JXM1,K)
                  QBF(N,2)=W(N,JMAX,K)
                enddo
                CALL FVSW(JMAX,RKX,RKY,QBF,F)!!!!!!!e

c-----------------RHS=-HE1*DT*(E(J+.5)-E(J-.5))+HE2*DQ(N-1)-------------

                DO J=1,JXM1
                  HE1DT=DT(J,K)*RJ(J,K)
                  DO N=1,5
                    S(N,J,K)=-HE1DT*(F(N,J+1)-F(N,J)) !!!!!!!e
                  ENDDO
                ENDDO
              ENDDO

c-----------------F TERM------------------------------------------------

      DO J=1,JXM1
        RKX=ETX(J,1)
        RKY=ETY(J,1)
        DO N=1,5
          QBF(N,1)=W(N,J,0)
          QBF(N,2)=W(N,J,1)
        ENDDO

        CALL FVSW(1,RKX,RKY,QBF,F) !!!!!!!f

        if(hmsthh .eq. 1)then

          DO  K=1,KXM2+2
            DO N=1,5
              uuu1(K,N)=W(N,J,K)
            ENDDO
          ENDDO

          CALL hmsth(uuu1,uur,uul,KXM2+2,alpha,beta,imethod)

          DO K=1,KXM2
            KP1=K+1

            RKX=ETX(J,KP1)
            RKY=ETY(J,KP1)
            DO N=1,5
              QBF(N,1)=uul(K,N)
              QBF(N,2)=uur(K,N)
            ENDDO
            CALL FVSW(KP1,RKX,RKY,QBF,F) !!!!!!!e

          ENDDO

        endif !	if(hmsthh .eq. 1)then

        RKX=ETX(J,KMAX)
        RKY=ETY(J,KMAX)
        DO N=1,5
          QBF(N,1)=W(N,J,KXM1)
          QBF(N,2)=W(N,J,KMAX)
        ENDDO

        CALL FVSW(KMAX,RKX,RKY,QBF,F) !!!!!!!f

c-----------------RHS=RHS-HE1*DT*((F(K+.5)-F(K-.5))---------------------

        DO K=1,KXM1
          HE1DT=DT(J,K)*RJ(J,K)
          DO N=1,5
            KP1=K+1
            S(N,J,K)=S(N,J,K)-HE1DT*(F(N,KP1)-F(N,K) ) !!!!!!!f
          ENDDO
        ENDDO
      ENDDO


c-----------------------------------------------------------------------  

      do n=1,5
        DO  K=-1,KMP1
          DO  J=-1,JMP1
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

c-----------------end--S(k/2)-------------------------------------------   

            elseif(split .EQ. 3)then

c-----------------MAIN TWO-SF(k)----------------------------------------   

c-----------------F TERM------------------------------------------------   

c-----------------EXPLICIT--2-------------------------------------------   

              do k=-1,kmp1
                do j=-1,jmp1
                  w(1,j,k)=Q(1,j,k)
                  w(2,j,k)=Q(2,J,K)/Q(1,J,K)
                  w(3,j,k)=Q(3,J,K)/Q(1,J,K)
                  w(4,j,k)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)
     #                    /Q(1,J,K)-q0*Q(5,J,K))
                  w(5,j,k)=Q(5,J,K)/Q(1,J,K)
                enddo
              enddo

c-----------------E TERM------------------------------------------------        

              DO K=1,KXM1
                RKX=XIX(1,K)
                RKY=XIY(1,K)
                DO N=1,5
                  QBF(N,1)=W(N,0,K)
                  QBF(N,2)=W(N,1,K)
                enddo

                CALL FVSW(1,RKX,RKY,QBF,F)

                if(hmsthh .eq. 1)then

                  DO J=1,JXM2+2
                    DO N=1,5
                      uuu1(J,N)=W(N,J,K)
                    ENDDO
                  ENDDO

                  CALL hmsth(uuu1,uur,uul,JXM2+2,alpha,beta,imethod)

                  DO J=1,JXM2
                    JP1=J+1
                    RKX=XIX(JP1,K)
                    RKY=XIY(JP1,K)
                      DO N=1,5
                        QBF(N,1)=uul(J,N)
                        QBF(N,2)=uur(J,N)
                      ENDDO
                    CALL FVSW(JP1,RKX,RKY,QBF,F)
                  ENDDO

                endif ! if(hmsthh .eq. 1)then

                RKX=XIX(JMAX,K)
                RKY=XIY(JMAX,K)
                DO N=1,5
                  QBF(N,1)=W(N,JXM1,K)
                  QBF(N,2)=W(N,JMAX,K)
                ENDDO
                
                CALL FVSW(JMAX,RKX,RKY,QBF,F)!!!!!!!e

c-----------------RHS=-HE1*DT*(E(J+.5)-E(J-.5))+HE2*DQ(N-1)-------------

                DO J=1,JXM1
                  HE1DT=DT(J,K)*RJ(J,K)
                  DO N=1,5
                    S(N,J,K)=-HE1DT*(F(N,J+1)-F(N,J)) !!!!!!!e
                  ENDDO
                ENDDO
              enddo


c-----------------F TERM------------------------------------------------   

              DO J=1,JXM1
                RKX=ETX(J,1)
                RKY=ETY(J,1)
                DO N=1,5
                  QBF(N,1)=W(N,J,0)
                  QBF(N,2)=W(N,J,1)
                ENDDO

                CALL FVSW(1,RKX,RKY,QBF,F) !!!!!!!f

                if(hmsthh .eq. 1)then

                  DO  K=1,KXM2+2
                    DO N=1,5
                      uuu1(K,N)=W(N,J,K)
                    ENDDO
                  ENDDO

                  CALL hmsth(uuu1,uur,uul,KXM2+2,alpha,beta,imethod)

                  DO K=1,KXM2
                    
                    KP1=K+1
                    RKX=ETX(J,KP1)
                    RKY=ETY(J,KP1)
                    DO N=1,5
                      QBF(N,1)=uul(K,N)
                      QBF(N,2)=uur(K,N)
                    ENDDO

                    CALL FVSW(KP1,RKX,RKY,QBF,F) !!!!!!!e

                  ENDDO
        
                endif ! if(hmsthh .eq. 1)then

                RKX=ETX(J,KMAX)
                RKY=ETY(J,KMAX)
                DO N=1,5
                  QBF(N,1)=W(N,J,KXM1)
                  QBF(N,2)=W(N,J,KMAX)
                enddo

                CALL FVSW(KMAX,RKX,RKY,QBF,F) !!!!!!!f

c-----------------RHS=RHS-HE1*DT*((F(K+.5)-F(K-.5))---------------------

                DO K=1,KXM1
                  HE1DT=DT(J,K)*RJ(J,K)
                  DO N=1,5
                    KP1=K+1
                    S(N,J,K)=S(N,J,K)-HE1DT*(F(N,KP1)-F(N,K)) !!!!!!!f
                  enddo
                enddo
              enddo              

c----------------------------------------------------------------------- 

              DO N=1,5
                DO K=-1,KMP1
                  DO J=-1,JMP1
                    Qold2d(J,K,N)=S(N,J,K)
                    Qold2ss(J,K,N)=Qold1s(J,K,N)+half*(Qold1d(J,K,N)+
     #                             Qold2d(J,K,N))!---3
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

c-----------------------------------------------------------------------

            elseif(split .EQ. 4)then

c-----------------source TWO-S(k/2)-------------------------------------


              do  n=1,5
                DO  K=-1,KMP1
                  DO  J=-1,JMP1
                    Qold2ss(j,k,n)=Q(n,J,k)
                  enddo
                enddo
              enddo

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  p(J,K)=GM1*(Qold2ss(J,K,4)- 
     #                   half*(Qold2ss(J,K,2)**2+Qold2ss(J,K,3)**2)/
     #                   Qold2ss(J,K,1)-q0*Qold2ss(J,K,5))
                enddo
              ENDDO

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  temperature(J,K)=p(J,K)/Qold2ss(J,K,1)
                enddo
              ENDDO

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  if ((temperature(J,K)-temperature0)>zero)then
                    hh(J,K)=one
                  else
                    hh(J,K)=zero  
                  endif
                enddo
              ENDDO

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  source(J,K,1)=0.0
                  source(J,K,2)=0.0
                  source(J,K,3)=0.0
                  source(J,K,4)=0.0
                  source(J,K,5)=-(one/ddd)*Qold2ss(J,K,5)*hh(J,K)
                enddo
              ENDDO

              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  source2(J,K,1)=0.0
                  source2(J,K,2)=0.0
                  source2(J,K,3)=0.0
                  source2(J,K,4)=0.0
                  source2(J,K,5)=-one*hh(J,K)
                enddo
              ENDDO

              do N=1,5
                DO  K=-1,KMP1
                  DO  J=-1,JMP1
                    Qolddss(J,K,N)=(half*DTMIN*source(J,K,N))/    
     #                             (one-half*half*DTMIN*source2(J,K,N))
                    Q(N,J,K)=Qold2ss(J,K,N)+Qolddss(J,K,N)    !---4
                  enddo
                enddo
              ENDDO

C-----------------SOLID SURFACE PHYSICAL B.C.-1-------------------------

c----------------LOWER WALL---------------------------------------------

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
                  Q(4,J,1)=PB*RGM1+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
                  Q(5,J,1)=DB*ZF
                enddo

                DO  N=1,5 
                  DO  J=1,jxm1
                    Q(N,J,0)=Q(N,J,1)
                    Q(N,J,-1)=Q(N,J,0)
                  enddo
                enddo

c----------------UP WALL------------------------------------------------

                DO  J=1,jxm1     
                  RKX=ETX(J,KXm2)
                  RKY=ETY(J,KXm2)
                  RSXY=1./DSQRT(RKX*RKX+RKY*RKY)
                  BKX=RKX*RSXY
                  BKY=RKY*RSXY

                  DR=Q(1,J,KXm2)
                  PR=PF(J,KXm2)
                  UR=UF(J,KXm2)
                  VR=VF(J,KXm2)
                  ZF=Q(5,J,KXm2)/DR
                  TERM=BKX*UR+BKY*VR

                  PB=PR
                  DB=DR
                  UB=UR-BKX*TERM
                  VB=VR-BKY*TERM
                  Q(1,J,kxm1)=DB
                  Q(2,J,kxm1)=DB*UB
                  Q(3,J,kxm1)=DB*VB
                  Q(4,J,kxm1)=PB*RGM1+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
                  Q(5,J,kxm1)=DB*ZF
                enddo

                DO  N=1,5 
                  DO  J=1,jxm1
                    Q(N,J,kmax)=Q(N,J,kxm1)
                    Q(N,J,kmp1)=Q(N,J,kmax)
                  enddo
                enddo

c----------------ghost point REFERENCE----------------------------------

                if(ghostpoint_open == 1)then
                  DO  N=1,5    
                    DO  K=2,KXM2
                      Q(N,2,K)=Q(N,3,K)
                      Q(N,1,K)=Q(N,2,K)
                      Q(N,0,K)=Q(N,1,K)
                      Q(N,-1,K)=Q(N,0,K)
                      Q(N,JXM2,K)=Q(N,JXM2-1,K)  
                      Q(N,JXM1,K)=Q(N,JXM2,K)  
                      Q(N,JMAX,K)=Q(N,JXM1,K)  
                      Q(N,Jmp1,K)=Q(N,jmax,K)
                    enddo
                  enddo
                endif

c----------------control z between 0 and 1------------------------------

                DO  K=-1,KMP1  
                  DO  J=-1,JMP1
                    ZFF=Q(5,J,K)/Q(1,J,K)
                    ZFF= min(1.0,ZFF)  
                    ZFF=max(0.0,ZFF)
                    DFF=Q(1,J,K)
                    UFF=Q(2,J,K)/Q(1,J,K)
                    VFF=Q(3,J,K)/Q(1,J,K)
                    PFF=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/
     #                  Q(1,J,K)-q0*Q(5,J,K))
                    Q(4,J,K)=PFF/(r-1)+0.5*DFF*(UFF**2+VFF**2)+
     #                       q0*DFF*ZFF
                    Q(5,J,K)=DFF*ZFF
                  enddo
                enddo

c-----------------------------------------------------------------------

              endif !if(split .EQ. 1)then
                    !if(split .EQ. 2)then
                    !if(split .EQ. 3)then
                    !if(split .EQ. 4)then

            enddo  !do split=1,4

          endif !if(splitting.EQ.1)then

!=======================================================================
! RESIDUAL
!=======================================================================

          RESID=0.

          DO K=1,KXM1
            DO J=1,JXM1
              RESID=RESID+S(1,J,K)*S(1,J,K)
            enddo
          enddo
          
          RESID=DSQRT(DABS(RESID/JK))

c-----------------------------------------------------------------------
c          IF (MOD(IT,ISAVE).EQ.0.or.it.eq.1) THEN
c            WRITE (*,*) IT,RESID,tot,dtmin
c          ENDIF
c-----------------------------------------------------------------------

          IF (MOD(IT,ISAVE).EQ.0.or.it.eq.1) THEN
c          IF (MOD(IT,ISAVE).EQ.0.) THEN  
            itt =it/ISAVE + 7000
c            WRITE (*,*) IT,RESID,q(2,4,2),q(1,4,2)
            write (filename,'("chan",I5,".plt")')itt
            open (unit=itt,file=filename,status='unknown')
            WRITE (itt,*)'Title="FLATE PLANE BY k-e MODEL SIMULATION"'
            WRITE (itt,*)'Variables = "X","Y","U","V","d","p","M","z" ' 
            WRITE (itt,*)'zone I=',JXm1,',J=',Kxm1,',f=point'

          DO  K=1,Kxm1
            DO  J=1,jxm1
              dom = dsqrt( UF(J,K)**2 + VF(J,K)**2)
              WRITE (itt,*) X(J,K),Y(J,K),UF(J,K),VF(J,K),Q(1,J,K),
     #              PF(j,k),dom,Q(5,J,K)/Q(1,J,K)
            enddo
          enddo

          close(itt)

c-----------------------------------------------------------------------
c          timeout=tot-(3d-8)
c          if( timeout > 0.0  )then
c            write(*,*) 'time out 1'
c          endif
c          timeout=tot-(92d-9)
c          if( timeout > 0.0  )then
c            write(*,*) 'time out 2'
c          endif
c          timeout=tot-(17d-8)
c          if( timeout > 0.0  )then
c            write(*,*) 'time out 3'
c          endif
c-----------------------------------------------------------------------

        ENDIF
      enddo

!=======================================================================
! PRINT/PLOT FINAL DATA
!=======================================================================
c-----------------------------------------------------------------------
c      OPEN (6,FILE='flow.dat')
c      DO N=1,5
c        DO J=-1,JMp1
c          DO K=-1,KMp1
c            WRITE (6,*) Q(N,J,K)
c          enddo
c        enddo
c      enddo
c      CLOSE(6)
c-----------------------------------------------------------------------

      return

      END

!=======================================================================
! Subroutine AFVS
! AUSMDV FLUX DIFFERENCE SPLITTING
!=======================================================================

      SUBROUTINE FVSW(I,xiw,eta,QA,F)  

      IMPLICIT NONE

      REAL*8 :: Q0,A1,A2,A12,F
      REAL*8 :: TXY,XIW,ETA
      REAL*8 :: D1,D2,QA,U1,U2,V1,V2,P1,P2,Z1,Z2,C1,C2
      REAL*8 :: PD1,PD2,UU1,UU2,CU1,CU2
      REAL*8 :: ALPHAL,ALPHAR
      REAL*8 :: CML,CMR,FLR
      REAL*8 :: AUSMDU,AUSMDV,AUSMVU,AUSMVV,DP,COEF,ST
      REAL*8 :: GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      INTEGER :: N,I,K

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      DIMENSION QA(6,2),F(6,3250),fLR(2,6)

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

c----------------2------------------------------------------------------

      d2=qa(1,2)
      u2=qa(2,2)
      v2=qa(3,2)
      p2=qa(4,2)
      z2=qa(5,2)

      c2=dsqrt(dabs(p2*1.4/qa(1,2)))
      pd2=p2/qa(1,2)

      uu2=(xiw*u2+eta*v2)

c----------------interface----------------------------------------------

      a12=(c1+c2)/2.
      cu1=uu1/txy/a12
      cu2=uu2/txy/a12
      alphaL=2*pd1/(pd1+pd2)
      alphaR=2*pd2/(pd1+pd2)

c----------------for particle based on density--------------------------

      IF (ABS(cu1).Le.1.0) THEN
        CML=alphaL*((cu1+1.)**2/4-(cu1+abs(cu1))/2.)
     #     +(cu1+abs(cu1))/2.
        a1=a1+CML*d1
        a2=a2+p1*(cu1+1.)**2*(2-cu1)/4.
      ELSE
        CML=(cu1+abs(cu1))/2.
        a1=a1+CML*d1
        a2=a2+p1*(cu1+abs(cu1))/(2*cu1)
      endif

c----------------2------------------------------------------------------

      IF (ABS(cu2).Le.1.0) THEN
        CMR=alphaR*(-(cu2-1.)**2/4-(cu2-abs(cu2))/2.)
     #     +(cu2-abs(cu2))/2.
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
      fLR(1,4)=(p1*3.5+.5*d1*
     $         (u1*u1+v1*v1) +q0*z1*d1)/d1
      fLR(1,5)=z1

      fLR(2,1)=1.
      fLR(2,2)=u2
      fLR(2,3)=v2
      fLR(2,4)=(p2*3.5+.5*d2*
     $         (u2*u2+v2*v2)+q0*z2*d2)/d2
      fLR(2,5)=z2
     
      do  k=1,1
        f(k,I)=.5*txy*a1*a12
     $        *(fLR(1,k)+fLR(2,k))
     $        -abs(.5*txy*a1)*a12
     $        *(fLR(2,k)-fLR(1,k))
      enddo

      do  k=5,5
        f(k,I)=.5*txy*a1*a12
     $        *(fLR(1,k)+fLR(2,k))
     $        -abs(.5*txy*a1)*a12
     $        *(fLR(2,k)-fLR(1,k))
      enddo

      ausmdU=.5*txy*a1*a12
     $      *(fLR(1,2)+fLR(2,2))
     $      -abs(.5*txy*a1)*a12
     $      *(fLR(2,2)-fLR(1,2))

      ausmdV=.5*txy*a1*a12
     $      *(fLR(1,3)+fLR(2,3))
     $      -abs(.5*txy*a1)*a12
     $      *(fLR(2,3)-fLR(1,3))

      ausmvU=(CML*d1*u1+CMR*d2*u2)*a12*txy 
      ausmvV=(CML*d1*v1+CMR*d2*v2)*a12*txy

      dp = p1-p2
      coef = 10.*abs(dp)/(p1+p2)
      st = .5*min(1.d0,coef)

      f(2,i)=(.5+st)*ausmvU+(.5-st)*ausmdU
      f(3,i)=(.5+st)*ausmvv+(.5-st)*ausmdv

c-----------------------------------------------------------------------

      do  k=4,4
        f(k,I)=.5*txy*a1*a12
     $        *(fLR(1,k)+fLR(2,k))
     $        -abs(.5*txy*a1)*a12
     $        *(fLR(2,k)-fLR(1,k))
      enddo

      f(2,I)=f(2,I)+xiw*a2
      f(3,I)=f(3,I)+eta*a2

      RETURN

c-----------------------------------------------------------------------

      END

!=======================================================================
! Subroutine GRD
!=======================================================================

      SUBROUTINE GRD(X,Y,JMAX,KMAX,ngrd)

      IMPLICIT NONE

      real*8 :: xini
      real*8 :: xend
      real*8 :: yini
      real*8 :: yend
      real*8 :: dx
      real*8 :: dy
      real*8 :: X(0:3250,0:650),Y(0:3250,0:650)      

      integer :: mx,i,JMAX
      integer :: my,j,KMAX
      integer :: ngrd

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
      
      RETURN

c-----------------------------------------------------------------------

      END 

!=======================================================================
! Subroutine ORI
!=======================================================================

      SUBROUTINE hmsth(u1,ur,ul,imax,alpha,beta,imethod)
      
      IMPLICIT NONE

      REAL*8 :: CK,HI,EPS,DTP,DTM,U1,SS,VPL,VPR,beta,EPSVOF,CFL,CBVOF
      REAL*8 :: DCOSH,TBVOF,DTANH,DIP0,DDP,UR,UL
      REAL*8 :: AVOFM,AVOF,SVOF,BB,CC,TBVOFP,TVOFP,TVOFM,VL,VR,ALPHA
      REAL*8 :: dvnm,dvnp,l,r,qmm1,qmm1e,qmp1,qmp1e,zeta
      
      INTEGER :: IMM1,IMAX,I,K,J,JP1,JM1,JP,JM,imm2,imethod

      DIMENSION vpr(-1:3250,5) ,vpl(-1:3250,5)
      DIMENSION vr(-1:3250,5) ,vl(-1:3250,5)
      DIMENSION ur(-1:3250,5) ,ul(-1:3250,5)
      DIMENSION u1(-1:3250,5)
      
c-----------------MUSCL-------------------------------------------------

      CK=1./3.
      imm1 = imax-1
      imm2 = imm1-1
      hi=1.
      eps=1e-14

      DO k=1,5
        DO J=1,imm1
          JP1=J+1

          DTP=u1(JP1,k)-u1(J,k)
          DTM=u1(J,k)-u1(J-1,k)

          SS=(2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)

          vpl(j,k)=u1(J,k)+hi*.25*SS*((1.-CK*SS)*DTM+(1.+CK*SS)*DTP)

          DTM=DTP
          DTP=u1(JP1+1,k)-u1(JP1,k)
          SS=(2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)

          vpr(j,k)=u1(JP1,k)-hi*.25*SS*((1.-CK*SS)*DTP+(1.+CK*SS)*DTM)

        enddo
      enddo  

c-----------------thincem-----------------------------------------------    

      epsvof = 1e-14
      CFL=1e-6

      cbvof = dcosh(beta)
      tbvof = dtanh(beta)

      DO k=1,5
        DO j=1,imm1

          jp=j+1
          jm=j-1

          jp = min(jp,imax)
          jm = max(jm,1)

          avofm = min(u1(jp,k),u1(jm,k))
          avof =  max(u1(jp,k),u1(jm,k))-avofm
          svof = -1.0

          if(avofm.eq.u1(jm,k)) svof = 1.0
            BB = exp(svof*beta*(2.0*(u1(j,k)-avofm+epsvof)/
     #         (avof+epsvof)-1.0))

            CC = (BB/cbvof - 1.0)/tbvof

            tbvofp = (CC+tbvof)/(1.0+CC*tbvof)
            tvofp = log( cosh(beta*cfl) - sinh(beta*cfl)*tbvofp)/
     #              (beta*(CFL+1e-12))
            tvofm = log( cosh(beta*cfl) + sinh(beta*cfl)*CC)/
     #            (beta*(CFL+1e-12))

            vl(j,k) = avofm + 0.5*avof*(1. - svof*tvofp)
            vr(j-1,k) = avofm + 0.5*avof*(1. + svof*tvofm)

          if(sign( 1d0,(u1(jp,k)-u1(j,k))*(u1(j,k)-u1(jm,k)) )
     #      .eq. -1d0) then
            vr(j-1,k) = u1(j,k)
            vl(j,k) = u1(j,k)
          endif

        ENDDO
      ENDDO

c-----------------combine-----------------------------------------------

c-----------------ORI---------------------------------------------------

      if (imethod .eq. 1) then

        eps =1e-12

        DO  i=1,imm1
          Jm1=i-1
          Jp1=i+1
          DTP= ( u1(Jm1,4) - 2*u1(i,4) + u1(Jp1,4) )
          DTM= ( u1(Jm1,4) + 2*u1(i,4) + u1(Jp1,4) )
          dip0 =  alpha * abs( ( eps+dtp)/(dtm+eps))
          ddp = max( (1d0 - dip0), 0. )
          DO  k=1,5

             ur(i,k)  = ddp * vpr(i,k)  + ( 1.- ddp ) * vr(i,k)
             ul(i,k) =  ddp * vpL(i,k)  + ( 1.- ddp ) * vl(i,k) 

          enddo
        enddo

c-----------------JCP---------------------------------------------------

      elseif (imethod .eq. 2) then

        epsvof = 1e-10

        do k=1,5
          do i=1,imm1
            Jm1=i-1
            Jp1=i+1

            r = vpl(i,k)-vpr(Jm1,k)+epsvof
            qmp1 = u1(Jp1,k)-u1(i,k)
            qmp1e = u1(Jp1,k)-u1(i,k)+epsvof
            dvnp = r/qmp1e
            
            l = vpl(i,k)-vpr(Jm1,k)+epsvof
            qmm1 = u1(i,k)-u1(Jm1,k)
            qmm1e = u1(i,k)-u1(Jm1,k)+epsvof
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

c-----------------ATM---------------------------------------------------

      elseif (imethod .eq. 3) then

        eps = 1e-12

        DO  k=1,5
          DO  i=1,imm1
            
            Jm1=i-1
            Jp1=i+1
            DTP= ( u1(Jm1,k) - 2*u1(i,k) + u1(Jp1,k) )
            DTM= ( u1(Jm1,k) + 2*u1(i,k) + u1(Jp1,k) )
            dip0 =  alpha * abs( ( eps+dtp)/(dtm+eps))
            ddp = max( (1d0 - dip0), 0. )
            
            ur(i,k) =  ddp * vpr(i,k)  + ( 1.- ddp ) * vr(i,k)
            ul(i,k) =  ddp * vpL(i,k)  + ( 1.- ddp ) * vl(i,k)

          enddo
        enddo

c-----------------MUSCL-------------------------------------------------

      elseif (imethod .eq. 4) then

        DO  k=1,5
          DO  i=1,imm1
            
            ur(i,k) =  vpr(i,k)
            ul(i,k) =  vpL(i,k)

          enddo
        enddo

c-----------------THINCEM-----------------------------------------------

      elseif (imethod .eq. 5) then

        DO  k=1,5
          DO  i=1,imm1
            
            ur(i,k) =  vr(i,k)
            ul(i,k) =  vl(i,k)

          enddo
        enddo

c-----------------------------------------------------------------------

      endif

c-----------------------------------------------------------------------

      RETURN

c-----------------------------------------------------------------------

      END
      
