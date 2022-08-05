!=======================================================================
!  hybrid  ausm-Roe  
!=======================================================================
      
      subroutine solver(alpha,beta,CFL,imethod,ISAVE,ITMAX,Q)

      IMPLICIT NONE

      real*8, intent(in)  :: alpha,beta,CFL

      INTEGER, intent(in) :: imethod,ISAVE,ITMAX

      real*8, intent(out) :: Q(6,-1:260,-1:110)

      real*8 :: RMACH,EPS,CK
      real*8 :: XIXCF,XIYCF,ETXCF,ETYCF
      real*8 :: XIXC,XIYC,ETXC,ETYC
      real*8 :: UF,VF,PF,VOL,RJ
      real*8 :: SONF,PAI,RGM,RGM1,GM12,C13,C43,PRL,PRT,VEO,SLT,FI,HI
      real*8 :: HE1,HE2,DEL,DT
      real*8 :: FMACH,VINF,UINF,PINF,FMACH2,EINF,TP
      real*8 :: Q1I,Q2I,Q3I,Q4I,ALT,PREF,TREF,DREF,CREF,VISREF,REFL
      real*8 :: X,Y,XIX,XIY,ETX,ETY
      real*8 :: DTMIN,DTMAX
      REAL*8 :: RKX,RKY,RSXY,PB,UB,VB
      REAL*8 :: DB,BKX,BKY,UR,VR,DR,PR
      REAL*8 :: C0,TERM
      REAL*8 :: QBF,F,S,GM
      REAL*8 :: DTP,DTM,DDP,SS
      REAL*8 :: Fr,Fa,XMU
      REAL*8 :: HE1DT,HIDT,AVJ,AFB,BFB,A,B,C,D,U,V,BV,RESID
      REAL*8 :: GM1,RGGM1,RG2M1

      INTEGER :: IT,IAUSM,ITT
      INTEGER :: KMAX,KMP1,KXM1,KXM2,KP1,KM1
      INTEGER :: JMAX,JMP1,JXM1,JXM2,JP1,JM1
      INTEGER :: J,K,N,M,JK

      CHARACTER*80 :: filename

c-----------------------------------------------------------------------

      PARAMETER (RMACH=340.1,EPS=1.E-6,CK=.333333)

      COMMON/GMD/ GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1
      COMMON/IJK/ JMAX,JXM1,JMP1,KMAX,KXM1,KMP1
      COMMON/iau/ iausm
      COMMON/GMD1/ DEL,UINF,VINF,PINF
      COMMON/SD/ S(6,0:260,0:110)
      COMMON/XY/ X(0:260,0:110),Y(0:260,0:110)
      COMMON/XEXY/ XIX(0:260,0:110),XIY(0:260,0:110),
     &             ETX(0:260,0:110),ETY(0:260,0:110)
      COMMON/JACO/ RJ(-1:260,-1:110)
     
      DIMENSION DT(0:260,0:110),QBF(6,2),F(6,1000)
      DIMENSION AFB(4,4),BFB(4,4),A(4,4),B(4),BV(4,4)
      DIMENSION Fa(6,1000),Fr(6,1000)

c-----------------STATEMENT FUNCTION------------------------------------

      XIXCF(J,K)=.5*(XIX(J,K)+XIX(J+1,K))
      XIYCF(J,K)=.5*(XIY(J,K)+XIY(J+1,K))
      ETXCF(J,K)=.5*(ETX(J,K)+ETX(J,K+1))
      ETYCF(J,K)=.5*(ETY(J,K)+ETY(J,K+1))
      UF(J,K)=Q(2,J,K)/Q(1,J,K)
      VF(J,K)=Q(3,J,K)/Q(1,J,K)
      PF(J,K)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K))
      SONF(J,K)=DSQRT(DABS(GM*PF(J,K)/Q(1,J,K)))
     
c-----------------PARAMETER SETUP---------------------------------------

      PAI=3.1415926535
      GM=1.4
      RGM=1./GM
      GM1=GM-1.
      RGM1=1./GM1
      GM12=2.*GM1

      C13=1./3.
      C43=4./3.
      PRL=.72
      PRT=.9
      VEO=13.2
      SLT=1.
      FI=0
      HE1=1.
      HE2=0.
      HI=HE1

c-----------------FREE STREAM CONDINTION--------------------------------

      DEL=1.
      FMACH=8.  

c-----------------FLUID-------------------------------------------------

      VINF=0.
      UINF=FMACH
      PINF=RGM
      FMACH2=FMACH*FMACH
      EINF=PINF*RGM1+.5*DEL*(UINF*UINF+VINF*VINF)
      TP=GM*PINF/DEL
      Q1I=DEL
      Q2I=DEL*UINF
      Q3I=DEL*VINF
      Q4I=EINF

c-----------------PARTICLE PARAMETER------------------------------------
      
      ALT=0.
      CALL SA (ALT,PREF,TREF,DREF,CREF,VISREF)
      REFL=0.1
        
c-----------------READ GRID DATA----------------------------------------

      OPEN (11,file='c2.grd',status='unknown')
      READ (11,*) JMAX,KMAX
      
      DO K=1,KMAX
        DO J=1,JMAX
           Read (11,*) X(J,K),Y(J,K)
        enddo
      enddo

      close(11)

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
      
      DO K=-1,KMp1
        DO J=-1,JMp1
          Q(1,J,K)=Q1I
          Q(2,J,K)=Q2I
          Q(3,J,K)=Q3I
          Q(4,J,K)=Q4I
        enddo
      enddo

c-----------------ETX,ETY-----------------------------------------------

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

c-----------------FAR FIELD NUMERICAL B.C.------------------------------

c-----------------OUTLET BOUNDARY---------------------------------------

        DO K=0,kmp1
          DO N=1,4
            Q(N,JMAX,K)=Q(N,JXM1,K)
            Q(N,JMP1,K)=Q(N,JMAX,K)
          ENDDO
        ENDDO

        DO K=0,kmp1
          DO N=1,4
            Q(N,0,K)=Q(N,1,K)
            Q(N,-1,K)=Q(N,0,K)
          ENDDO
        ENDDO

c-----------------SOLID SURFACE-----------------------------------------

c-----------------Lower BOUNDARY----------------------------------------

        DO  J = 0,jmp1
          RKX=ETX(J,1)
          RKY=ETY(J,1)
          RSXY=1./DSQRT(RKX*RKX+RKY*RKY)
          BKX=RKX*RSXY
          BKY=RKY*RSXY
          C0=SONF(J,0)
          DR=Q(1,J,1)
          PR=PF(J,1)
          UR=UF(J,1)
          VR=VF(J,1)

          UB=Ur-RKX*(BKX*Ur+BKY*Vr)
          VB=Vr-RKY*(BKX*Ur+BKY*Vr)
          TERM=BKX*UR+BKY*VR
          PB=PR-Q(1,J,0)*C0*TERM
          DB=DR+(PB-PR)/C0**2

          Q(1,J,0)=Db
          Q(2,J,0)=DB*Ub
          Q(3,J,0)=DB*Vb
          Q(4,J,0)=pb/.4+.5*db*(ur*ur+vr*vr)
        enddo

        DO N=1,4
          DO  J=-1,jmp1
            Q(N,J,-1)=Q(N,J,0)
          ENDDO
        ENDDO

c-----------------EXPLICIT---------------------------------------------- 
     
        DO K=1,KXM1
          RKX=XIX(1,K)
          RKY=XIY(1,K)
          DO N=1,4
            QBF(N,1)=Q(N,0,K)
            QBF(N,2)=Q(N,1,K)
          enddo

          CALL ausm(1,RKX,RKY,QBF,F)

          DO J=1,JXM2
            JP1=J+1
            RKX=XIX(JP1,K)
            RKY=XIY(JP1,K)
            DO N=1,4
              DTP=Q(N,JP1,K)-Q(N,J,K)
              DTM=Q(N,J,K)-Q(N,J-1,K)
              SS=1.+SLT*((2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)-1.)
              QBF(N,1)=Q(N,J,K)+FI*.25*SS*((1.-CK*SS)
     #                *DTM+(1.+CK*SS)*DTP)
              DTM=DTP
              DTP=Q(N,JP1+1,K)-Q(N,JP1,K)
              SS=1.+SLT*((2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)-1.)
              QBF(N,2)=Q(N,JP1,K)-FI*.25*SS*((1.-CK*SS)
     #                *DTP+(1.+CK*SS)*DTM)         
            enddo

c-----------------hybrid------------------------------------------------ 

            CALL fVS(JP1,RKX,RKY,QBF,Fr,beta)
            CALL ausm(JP1,RKX,RKY,QBF,Fa)

            DO N=1,4
              DTP= ( pf(jp1-1,k) - 2*pf(jp1,k) + pf(jp1+1,k) )
              DTM= ( pf(jp1-1,k) + 2*pf(jp1,k) + pf(jp1+1,k) )
              ddp =  alpha * abs( dtp/dtm) 
              xmu = max( (1d0 - ddp), 0. )

              if (imethod .eq. 1) then
                f(N,jp1)  = xmu * fa(n,jp1)  + (1-xmu) * fr(n,jp1)
              elseif (imethod .eq. 2) then
                f(N,jp1)  = fr(n,jp1)
              elseif (imethod .eq. 3) then
                f(N,jp1)  = fa(n,jp1)
              endif

            enddo
          enddo

          RKX=XIX(JMAX,K)
          RKY=XIY(JMAX,K)

          DO N=1,4
            QBF(N,1)=Q(N,JXM1,K)
            QBF(N,2)=Q(N,JMAX,K)
          enddo

          CALL ausm(JMAX,RKX,RKY,QBF,F)

c-----------------RHS=-HE1*DT*(E(J+.5)-E(J-.5))-------------------------

          DO J=1,JXM1
            HE1DT=HE1*DT(J,K)*RJ(J,K)
            DO N=1,4
              S(N,J,K)=-HE1DT*(F(N,J+1)-F(N,J))
            enddo
          enddo
        enddo

c-----------------F TERM------------------------------------------------

        DO J=1,JXM1
          RKX=ETX(J,1)
          RKY=ETY(J,1)
          DO N=1,4
            QBF(N,1)=Q(N,J,0)
            QBF(N,2)=Q(N,J,1)
          enddo

          CALL ausm(1,RKX,RKY,QBF,F)

          DO K=1,KXM2
            KP1=K+1
            RKX=ETX(J,KP1)
            RKY=ETY(J,KP1)
            DO N=1,4
              DTP=Q(N,J,KP1)-Q(N,J,K)
              DTM=Q(N,J,K)-Q(N,J,K-1)
              SS=1.+SLT*((2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)-1.)
              QBF(N,1)=Q(N,J,K)+FI*.25*SS*((1.-CK*SS)
     #                *DTM+(1.+CK*SS)*DTP)
              DTM=DTP
              DTP=Q(N,J,KP1+1)-Q(N,J,KP1)
              SS=1.+SLT*((2.*DTP*DTM+EPS)/(DTP*DTP+DTM*DTM+EPS)-1.)
              QBF(N,2)=Q(N,J,KP1)-FI*.25*SS*((1.-CK*SS)
     #                *DTP+(1.+CK*SS)*DTM)
            enddo
 
c-----------------hybrid------------------------------------------------

            CALL fVS(KP1,RKX,RKY,QBF,Fr,beta)
            CALL ausm(KP1,RKX,RKY,QBF,Fa)

            DO N=1,4
              DTP= ( pf(j,KP1-1) - 2*pf(j,KP1) + pf(j,KP1+1) )
              DTM= ( pf(j,KP1-1) + 2*pf(j,KP1) + pf(j,KP1+1) )
              ddp =  alpha * abs( dtp/dtm) 
              xmu = max( (1d0 - ddp), 0. )

              if (imethod .eq. 1) then
                f(N,Kp1)  = xmu * fa(n,kp1)  + (1-xmu) * fr(n,kp1)
              elseif (imethod .eq. 2) then
                f(N,Kp1)  = fr(n,kp1)
              elseif (imethod .eq. 3) then
                f(N,Kp1)  = fa(n,kp1)
              endif

            enddo
          enddo

          RKX=ETX(J,KMAX)
          RKY=ETY(J,KMAX)

          DO N=1,4
            QBF(N,1)=Q(N,J,KXM1)
            QBF(N,2)=Q(N,J,KMAX)
          enddo

          CALL ausm(KMAX,RKX,RKY,QBF,F)

c-----------------RHS=-HE1*DT*(E(J+.5)-E(J-.5))-------------------------

          DO K=1,KXM1
            HE1DT=HE1*DT(J,K)*RJ(J,K)
            DO N=1,4
              KP1=K+1
              S(N,J,K)=S(N,J,K)-HE1DT*(F(N,KP1)-F(N,K))
            enddo
          enddo
        enddo

c-----------------------------------------------------------------------  

        DO K=1,KXM1
          KM1=K-1
          KP1=K+1
          AVJ=RJ(1,K)
          RKX=AVJ*XIX(1,K)
          RKY=AVJ*XIY(1,K)

          CALL JMAB(RKX,RKY,Q(1,0,K),
     #         UF(0,K),VF(0,K),SONF(0,K),1,AFB)
          
          DO J=1,JXM1
            JM1=J-1
            JP1=J+1
            AVJ=.5*(RJ(J,KM1)+RJ(J,K))
            RKX=AVJ*ETX(J,K)
            RKY=AVJ*ETY(J,K)

            CALL JMAB(RKX,RKY,Q(1,J,KM1),
     #           UF(J,KM1),VF(J,KM1),SONF(J,KM1),1,BFB)
            
            HIDT=HI*DT(J,K)

            DO N=1,4
              B(N)=S(N,J,K)+HIDT*
     &            (AFB(N,1)*S(1,JM1,K)+(BFB(N,1)-.5*BV(N,1))*S(1,J,KM1)+
     &             AFB(N,2)*S(2,JM1,K)+(BFB(N,2)-.5*BV(N,2))*S(2,J,KM1)+
     &             AFB(N,3)*S(3,JM1,K)+(BFB(N,3)-.5*BV(N,3))*S(3,J,KM1)+
     &             AFB(N,4)*S(4,JM1,K)+(BFB(N,4)-.5*BV(N,4))*S(4,J,KM1))
            enddo

            D=Q(1,J,K)
            U=UF(J,K)
            V=VF(J,K)
            C=SONF(J,K)
            AVJ=.5*(RJ(J,K)+RJ(JP1,K))
            RKX=AVJ*XIX(JP1,K)
            RKY=AVJ*XIY(JP1,K)

            CALL JMAB(RKX,RKY,D,U,V,C,1,AFB)

            AVJ=.5*(RJ(J,K)+RJ(J,KP1))
            RKX=AVJ*ETX(J,KP1)
            RKY=AVJ*ETY(J,KP1)

            CALL JMAB(RKX,RKY,D,U,V,C,1,BFB)
           
            DO N=1,4
              DO M=1,4
                A(N,M)=HIDT*(AFB(N,M)+BFB(N,M))
              enddo
              A(N,N)=1.+A(N,N)
            enddo

            CALL SLES(J,K,A,B)

          enddo
        enddo

c----------------------------------------------------------------------- 

        DO K=KXM1,1,-1
          KM1=K-1
          KP1=K+1
          AVJ=RJ(JXM1,K)
          RKX=AVJ*XIX(JMAX,K)
          RKY=AVJ*XIY(JMAX,K)

          CALL JMAB(RKX,RKY,Q(1,JMAX,K),
     #         UF(JMAX,K),VF(JMAX,K),SONF(JMAX,K),2,AFB)

          DO J=JXM1,1,-1
            JM1=J-1
            JP1=J+1
            AVJ=.5*(RJ(J,K)+RJ(J,KP1))
            RKX=AVJ*ETX(J,KP1)
            RKY=AVJ*ETY(J,KP1)

            CALL JMAB(RKX,RKY,Q(1,J,KP1),
     #           UF(J,KP1),VF(J,KP1),SONF(J,KP1),2,BFB)
        
            HIDT=HI*DT(J,K)

            DO N=1,4
              B(N)=S(N,J,K)-HIDT*
     &            (AFB(N,1)*S(1,JP1,K)+(BFB(N,1)-.5*BV(N,1))*S(1,J,KP1)+
     &             AFB(N,2)*S(2,JP1,K)+(BFB(N,2)-.5*BV(N,2))*S(2,J,KP1)+
     &             AFB(N,3)*S(3,JP1,K)+(BFB(N,3)-.5*BV(N,3))*S(3,J,KP1)+
     &             AFB(N,4)*S(4,JP1,K)+(BFB(N,4)-.5*BV(N,4))*S(4,J,KP1))
            enddo

            D=Q(1,J,K)
            U=UF(J,K)
            V=VF(J,K)
            C=SONF(J,K)
            AVJ=.5*(RJ(JM1,K)+RJ(J,K))
            RKX=AVJ*XIX(J,K)
            RKY=AVJ*XIY(J,K)

            CALL JMAB(RKX,RKY,D,U,V,C,2,AFB)

            AVJ=.5*(RJ(J,KM1)+RJ(J,K))
            RKX=AVJ*ETX(J,K)
            RKY=AVJ*ETY(J,K)

            CALL JMAB(RKX,RKY,D,U,V,C,2,BFB)

            DO N=1,4
              DO M=1,4
                A(N,M)=-HIDT*(AFB(N,M)+BFB(N,M))
              enddo
              A(N,N)=1.+A(N,N)
            enddo

            CALL SLES(J,K,A,B)

          enddo
        enddo

c-----------------Q(N+1)=Q(N)+DQ----------------------------------------  

        DO N=1,4
          DO K=1,KXM1
            DO J=1,JXM1
              Q(N,J,K)=Q(N,J,K)+S(N,J,K)
            enddo
          enddo
        enddo

!=======================================================================
!  RESIDUAL
!=======================================================================

        RESID=0.

        DO K=1,KXM1
          DO J=1,JXM1
            RESID=RESID+S(1,J,K)*S(1,J,K)
          enddo
        enddo
        
        JK=jxm1*kxm1
        RESID=DSQRT(DABS(RESID/JK))
        
c        IF (MOD(IT,ISAVE).EQ.0.or.it.eq.1) THEN
c          WRITE (*,*) IT,RESID,q(2,30,2)
c        ENDIF
        
c        IF (RESID.LT.1.E-7) then
c          OPEN (6,FILE='flow-3b.dat',status='unknown')
c          DO N=1,4
c            DO J=-1,JMP1
c              DO K=-1,KMP1
c                WRITE (6,*) Q(N,J,K)
c              enddo
c            enddo
c          enddo
c          close(6)
c        endif


!=======================================================================
!  PRINT/PLOT FINAL DATA
!=======================================================================
      
        IF (MOD(IT,ISAVE).EQ.0.) THEN
          itt = it/ISAVE
          write (filename,'("HAR",I5,".plt")')itt
          open (unit=itt,file=filename,status='unknown')
          WRITE (itt,*)'Title="FLATE PLANE BY k-e MODEL SIMULATION"'
          WRITE (itt,*)'Variables = "X","Y","UF","VF","del","vel","pre"'
          WRITE (itt,*)'zone I=',JMAX,',J=',KMAX,',f=point'

          DO K=1,KMAX
            DO J=1,JMAX
              WRITE (itt,*) X(J,K),Y(J,K),UF(J,K),VF(J,K),Q(1,J,K),
     &                     dsqrt(uf(j,k)**2+vf(j,k)**2),PF(J,K)
            enddo
          enddo

        close(itt)

        ENDIF
        
c-----------------------------------------------------------------------

      enddo

c      write(*,*) 'converging!',it

      RETURN

      END

!=======================================================================
! Subroutine SA
! THE ARDC MODEL ATMOSPHERE,1959(.LE.91.KM)
!
!  INPUT:
!  ALT:ALTITUDE(KM)
!
!  OUTPUT:
!  P:PRESSURE(N/M2)
!  T:TEMPERATURE(K)
!  D:DENSITY(KG/M3)
!  C:SOUND SPEED(M/S)
!  VIS:VISCOSITY(KG/M-S)
!
!=======================================================================

      SUBROUTINE SA(ALT,P,T,D,C,VIS)

      IMPLICIT NONE

      real*8 :: ALT,P,T,D,C,VIS
      real*8 :: B1,B2,B3,CH,CHI
      real*8 :: DI,HI,RI,TI,TM,TMI

      INTEGER :: I

      DIMENSION HI(6),DI(6),CHI(6),TMI(6),TI(6),
     #          B1(6),B2(6),B3(6),RI(6)
      
      DATA HI/2.4,17.8,41.4,48.2,72.,80./,
     #     DI/.9667,.1225,.3293E-2,.1342E-2,.7542E-4,.212E-5/,
     #     CHI/7.984,6.3772,7.8596,8.3994,8.3416,4.972/,
     #     TMI/272.57,216.66,265.06,282.66,200.79,165.7/,
     #     TI/272.57,216.66,265.06,282.66,200.79,165.7/,
     #     B1/-.187708,.002,.090063,.002624,-.1296385,.00154545/,
     #     B2/-6.48034,0.,2.96406,0.,-4.404321,0./,
     #     B3/-6.48034,0.,2.9641,0.,-4.39464,0./,
     #     RI/11.,25.,47.4,53.,80.,91./

      IF(ALT.GT.91.)THEN
        WRITE(*,*)'I/O ERROR:ALT > 91.(KM)'
        STOP
      ENDIF
      
      DO I=1,6
        IF(ALT.LE.RI(I)) THEN
          CH=CHI(I)+B1(I)*(ALT-HI(I))
          TM=TMI(I)+B2(I)*(ALT-HI(I))
          D=DI(I)*(TMI(I)/TM)*(CHI(I)/CH)**(1./B1(I))
          T=TI(I)+B3(I)*(ALT-HI(I))
          P=.2870396E3*D*TM
          C=20.0463*DSQRT(DABS(TM))
          VIS=(1.458E-6*DSQRT(DABS(T**3)))/(T+110.4)
        ENDIF
      ENDDO

      RETURN

      END

!=======================================================================
! Subroutine JMAB
! JACOBIAN MATRICES: AF / AB / BF / BB
!=======================================================================

      SUBROUTINE JMAB(RKX,RKY,D,U,V,C,L,AB)
 
      IMPLICIT NONE

      real*8 :: RKX,RKY,BKX,BKXC,BKY,BKYC,CSXY
      real*8 :: D,U,V,C,AB,BT,C2,CBST
      real*8 :: AEV1,AEV3,AEV4
      real*8 :: AF,ED1,EV,EV1,EV2,EV3,EV4
      real*8 :: FC2RG1,FI2
      real*8 :: GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1
      real*8 :: GM1U,GM1V,RC2,RD,RSXY,S2
      real*8 :: SXY,T1,T2
      
      INTEGER :: M,N,L

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1
      
      DIMENSION EV(3),T1(4,4),T2(4,4),AB(4,4)

c-----------------EV(+/-)-----------------------------------------------

c-----------------EV----------------------------------------------------

      EV1=RKX*U+RKY*V
      SXY=DSQRT(DABS(RKX*RKX+RKY*RKY))
      CSXY=C*SXY
      EV3=EV1+CSXY
      EV4=EV1-CSXY
      AEV1=DABS(EV1)
      AEV3=DABS(EV3)
      AEV4=DABS(EV4)

      IF(L.EQ.1)THEN

c-----------------(EV+)-------------------------------------------------

        EV(1)=.5*(EV1+AEV1)
        EV(2)=.5*(EV3+AEV3)
        EV(3)=.5*(EV4+AEV4)

      ELSE

c-----------------(EV-)-------------------------------------------------

        EV(1)=.5*(EV1-AEV1)
        EV(2)=.5*(EV3-AEV3)
        EV(3)=.5*(EV4-AEV4)

      ENDIF

c-----------------EV*TINV-----------------------------------------------

      FI2=.5*GM1*(U*U+V*V)
      C2=C*C
      RC2=1./C2
      RD=1./D
      RSXY=1./SXY
      BKX=RKX*RSXY
      BKY=RKY*RSXY
      ED1=2.
      S2=DSQRT(ED1)
      BT=1./(S2*D*C)
      CBST=C*(BKX*U+BKY*V)
      BKXC=BKX*C
      BKYC=BKY*C
      GM1U=GM1*U
      GM1V=GM1*V
      EV1=EV(1)
      T2(1,1)=EV1*(1.-FI2*RC2)
      T2(1,2)=EV1*GM1*RC2*U
      T2(1,3)=EV1*GM1*RC2*V
      T2(1,4)=-EV1*GM1*RC2
      EV2=EV1
      T2(2,1)=-EV2*RD*(BKY*U-BKX*V)
      T2(2,2)=EV2*BKY*RD
      T2(2,3)=-EV2*BKX*RD
      T2(2,4)=0.
      EV3=EV(2)
      T2(3,1)=EV3*BT*(FI2-CBST)
      T2(3,2)=EV3*BT*(BKXC-GM1U)
      T2(3,3)=EV3*BT*(BKYC-GM1V)
      T2(3,4)=EV3*BT*GM1
      EV4=EV(3)
      T2(4,1)=EV4*BT*(FI2+CBST)
      T2(4,2)=-EV4*BT*(BKXC+GM1U)
      T2(4,3)=-EV4*BT*(BKYC+GM1V)
      T2(4,4)=EV4*BT*GM1

c-----------------T-----------------------------------------------------

      AF=D/(S2*C)
      FC2RG1=(FI2+C2)*RGM1
      T1(1,1)=1.
      T1(1,2)=0.
      T1(1,3)=AF
      T1(1,4)=AF
      T1(2,1)=U
      T1(2,2)=BKY*D
      T1(2,3)=AF*(U+BKXC)
      T1(2,4)=AF*(U-BKXC)
      T1(3,1)=V
      T1(3,2)=-BKX*D
      T1(3,3)=AF*(V+BKYC)
      T1(3,4)=AF*(V-BKYC)
      T1(4,1)=FI2*RGM1
      T1(4,2)=D*(BKY*U-BKX*V)
      T1(4,3)=AF*(FC2RG1+CBST)
      T1(4,4)=AF*(FC2RG1-CBST)

c-----------------AB=T*EV*TINV------------------------------------------

      DO N=1,4
        DO M=1,4
          AB(N,M)=T1(N,1)*T2(1,M)+T1(N,2)*T2(2,M)
     &           +T1(N,3)*T2(3,M)+T1(N,4)*T2(4,M)
        enddo
      enddo

c-----------------------------------------------------------------------

      RETURN
      END

!=======================================================================
! Subroutine SLES
! SOLVE LINEAR EQS BY LU METHOD
!=======================================================================

      SUBROUTINE SLES(J,K,A,B)

      IMPLICIT NONE

      real*8 :: A,B,S
      real*8 :: D1,D2,D3,D4
      real*8 :: RL11,RL21,RL31,RL41
      real*8 :: RL22,RL32,RL42,RL33
      real*8 :: RL43,RL44
      real*8 :: RU12,RU13,RU14,RU23
      real*8 :: RU24,RU34

      INTEGER :: J,K

      COMMON/SD/S(6,0:260,0:110)

      DIMENSION A(4,4),B(4)

      RL11=1./A(1,1)
      RL21=A(2,1)
      RU12=A(1,2)*RL11
      RL22=1./(A(2,2)-RL21*RU12)
      RU13=A(1,3)*RL11
      RU14=A(1,4)*RL11
      RL31=A(3,1)
      RL32=A(3,2)-RL31*RU12
      RU23=(A(2,3)-RL21*RU13)*RL22
      RL33=1./(A(3,3)-RU13*RL31-RU23*RL32)
      RU24=(A(2,4)-RL21*RU14)*RL22
      RL41=A(4,1)
      RL42=A(4,2)-RL41*RU12
      RL43=A(4,3)-RL41*RU13-RL42*RU23
      RU34=(A(3,4)-RL31*RU14-RL32*RU24)*RL33
      RL44=1./(A(4,4)-RU14*RL41-RU24*RL42-RU34*RL43)
      D1=RL11*B(1)
      D2=RL22*(B(2)-RL21*D1)
      D3=RL33*(B(3)-RL31*D1-RL32*D2)
      D4=RL44*(B(4)-RL41*D1-RL42*D2-RL43*D3)
      S(4,J,K)=D4
      S(3,J,K)=D3-RU34*D4
      S(2,J,K)=D2-RU23*S(3,J,K)-RU24*D4
      S(1,J,K)=D1-RU12*S(2,J,K)-RU13*S(3,J,K)-RU14*D4

c-----------------------------------------------------------------------

      RETURN
      END

!=======================================================================
! Subroutine AUSM
! flux splitting  
!=======================================================================

      SUBROUTINE ausm(I,rkx,rky,QA,F)

      IMPLICIT NONE

      real*8 :: RKX,RKY,QA,F,AHALF
      real*8 :: ALL,ALR,AR,AL,BTL,BTR
      real*8 :: DEL,UINF,VINF,PINF
      real*8 :: DELTA,DENK2,DENKL
      real*8 :: DER,DKL,DKP,DKR,DL,DR
      real*8 :: FML,FMR,HL,HR,omemins,omeplus
      real*8 :: pdlr,PL,PR,pnet,ppl,ppr,Tgp
      real*8 :: UCL,UCR,UL,UR,VL,VR
      real*8 :: xmaddm,xmaddp,xmbsl,xmbsr,xmface,xmlp
      real*8 :: xmmins,xmml,xmmr,xmplus,xmrp

      INTEGER :: N,I
      
      COMMON /GMD1/DEL,UINF,VINF,PINF

      DIMENSION QA(6,2),F(6,1000)

      data delta, dkp/2.,1./

      DO N=1,6
        F(N,I)=0.
      enddo

      tgp=dsqrt(rkx*rkx+rky*rky)

      dL=qa(1,1)
      uL=qa(2,1)/qa(1,1)
      vL=qa(3,1)/qa(1,1)
      dkL=qa(5,1)/qa(1,1)
      deL=qa(6,1)/qa(1,1)
      denkl=.66667* qa(5,1)
      pL=0.4*(qa(4,1)-.5*qa(1,1)*(uL*uL+vL*vL))
      aL=dsqrt(dabs(pL*1.4/qa(1,1)))
      hL=pL*3.5/dL+.5*(uL*uL+vL*vL)
      ucL=(rkx*uL+rky*vL)/tgp

c----------------2------------------------------------------------------

      dr=qa(1,2)
      ur=qa(2,2)/qa(1,2)
      vr=qa(3,2)/qa(1,2)
      dkr=qa(5,2)/qa(1,2)
      der=qa(6,2)/qa(1,2)
      denk2=.66667* qa(5,2)
      pr=0.4*(qa(4,2)-.5*qa(1,2)*(ur*ur+vr*vr))
      Ar=dsqrt(dabs(pr*1.4/qa(1,2)))
      hr=pr*3.5/dr+.5*(ur*ur+vr*vr)
      ucR=(rkx*ur+rky*vr)/tgp

c----------------interface----------------------------------------------

      ahalf = 0.5*(ar + aL)

      xmlp = ucl/ahalf
      xmrp = ucr/ahalf

c----------------ausmd--------------------------------------------------

c----------------switching functions------------------------------------

      all = 0.5*(1. + dsign(1.d0,xmlp))  
      alr = 0.5*(1. - dsign(1.d0,xmrp))  

      btl = -dmax1(0.d0,1.d0-dint(abs(xmlp)))
      btr = -dmax1(0.d0,1.d0-dint(abs(xmrp)))

c----------------split and interface Mach numbers-----------------------

      xmml = 0.25*(xmlp+1.0)**2 
      xmmr = -0.25*(xmrp-1.0)**2 
      xmbsl = 0.5*(xmlp + dabs(xmlp))
      xmbsr = 0.5*(xmrp - dabs(xmrp))
      xmaddp = btl*(xmbsl - xmml)
      xmaddm = btr*(xmbsr - xmmr)
      
      pdlr=.5*(pL/dL+ pR/dR)

      omeplus = pL/dL/pdlr
      omemins = pR/dR/pdlr
     
      xmplus = xmbsl + omeplus*xmaddp 
      xmmins = xmbsr + omemins*xmaddm 
      xmface = ahalf*(dl*xmplus + dr*xmmins)


c----------------interface mass flux------------------------------------

      fml = Tgp*0.5*(xmface + dabs(xmface))
      fmr = Tgp*0.5*(xmface - dabs(xmface))

c----------------pressure splitting-------------------------------------

      ppl = 0.25*(xmlp+1.)**2*(2.-xmlp)
      ppr = 0.25*(xmrp-1.)**2*(2.+xmrp)
      pnet = (all*(1.+btl) - btl*ppl)*PL + (alr*(1.+btr) - btr*ppr)*PR

      f(1,i) =(fml + fmr)
      f(2,i) = (fml*uL + fmr*ur + rkx*pnet)
      f(3,i) = (fml*vL + fmr*vr + rky*pnet)
      f(4,i) = (fml*hL + fmr*hr)

c-----------------------------------------------------------------------

      RETURN
      END

!=======================================================================
! Subroutine FSV
! ROE FLUX DIFFERENCE SPLITTING
!=======================================================================

      SUBROUTINE fvs(I,RKX,RKY,QLR,F,beta)

      IMPLICIT NONE

      real*8 :: rkx,rky,qlr,f,ab,AEV
      real*8 :: beta,BKX,bky,c,d,dl,dr,du,dv,e,epp
      real*8 :: ev,FLPR,h,hlr,p,qrml,RSXY,sdr,sig
      real*8 :: GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1
      real*8 :: sxy,u,uvc,v

      INTEGER :: I,L,N

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      DIMENSION QLR(6,2),FLPR(4),HLR(2),QRML(4),EV(3),AB(4,4),F(6,1000)

      DO N=1,4
        FLPR(N)=0.
        QRML(N)=0.
      enddo

c---------------- FL+FR,HLR,QR-QL---------------------------------------

      SIG=1.
      DO L=1,2
        D=QLR(1,L)
        U=QLR(2,L)/D
        V=QLR(3,L)/D
        P=0.4*(QLR(4,L)-.5*D*(u*u+v*v))
        UVC=RKX*U+RKY*V
        DU=QLR(2,L)
        DV=QLR(3,L)
        FLPR(1)=FLPR(1)+D*UVC
        FLPR(2)=FLPR(2)+DU*UVC+RKX*P
        FLPR(3)=FLPR(3)+DV*UVC+RKY*P
        E=P*RGM1+.5*D*(U*U+V*V)
        EPP=E+P
        FLPR(4)=FLPR(4)+EPP*UVC

        HLR(L)=EPP/D

        SIG=-1.*SIG
        QRML(1)=QRML(1)+SIG*D
        QRML(2)=QRML(2)+SIG*DU
        QRML(3)=QRML(3)+SIG*DV
        QRML(4)=QRML(4)+SIG*E
      enddo

c----------------ROE AVERAGE--------------------------------------------

      DL=QLR(1,1)
      DR=QLR(1,2)
      SDR=DSQRT(DR/DL)
      D=(QLR(1,2)+SDR*QLR(1,1))/(1.+SDR)
      U=(QLR(2,1)/DL+SDR*QLR(2,2)/DR)/(1.+SDR)
      V=(QLR(3,1)/DL+SDR*QLR(3,2)/DR)/(1.+SDR)
      H=(HLR(1)+SDR*HLR(2))/(1.+SDR)
      C=DSQRT(GM1*(H-.5*(U*U+V*V)))

c----------------EV-----------------------------------------------------

      SXY=DSQRT(RKX*RKX+RKY*RKY)
      RSXY=1./SXY
      BKX=RKX*RSXY
      BKY=RKY*RSXY
      UVC=BKX*U+BKY*V
      EV(1)=UVC
      EV(2)=UVC+C
      EV(3)=UVC-C


!=======================================================================
! ENTROPY CONDITION:
! IF AEV.LT.E EV=(AEV**2+E**2)/(2E) WHERE E=.125
!=======================================================================

c----------------EV/RJ--------------------------------------------------

      DO N=1,3

        AEV=DABS(EV(N))

c----------------beta---------------------------------------------------

        IF(AEV.LT.beta)THEN
          EV(N)=SXY*(4.*EV(N)*EV(N)+.0625)
        ELSE
          EV(N)=SXY*AEV
        ENDIF

      enddo

c----------------JACOBIAN MATRICES A / B--------------------------------

      CALL ROEM(BKX,BKY,D,U,V,C,EV,AB)

c----------------FLUX---------------------------------------------------

      DO N=1,4
        F(N,I)=.5*(FLPR(N)-(AB(N,1)*QRML(1)+AB(N,2)*QRML(2)
     #        +AB(N,3)*QRML(3)+AB(N,4)*QRML(4)))
      enddo

c-----------------------------------------------------------------------

      RETURN
      END

!=======================================================================
! Subroutine ROEM
! JACOBIAN MATRICES A / B
!=======================================================================

      SUBROUTINE ROEM(BKX,BKY,D,U,V,C,EV,AB)

      IMPLICIT NONE

      real*8 :: ab
      real*8 :: BKX,bky,c,d
      real*8 :: ev
      real*8 :: GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1
      real*8 :: u,v

      real*8 :: af,BKXC,BKYC,bt,c2,cbst,ed1,ev1,ev2,ev3,ev4
      real*8 :: FC2RG1,FI2,GM1U,GM1V,RC2,rd,s2,t1,t2

      INTEGER :: M,N

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      DIMENSION EV(3),T1(4,4),T2(4,4),AB(4,4)
      
c-----------------EV*TINV-----------------------------------------------

      FI2=.5*GM1*(U*U+V*V)
      C2=C*C
      RC2=1./C2
      RD=1./D
      ED1=2.
      S2=DSQRT(ED1)
      BT=1./(S2*D*C)
      CBST=C*(BKX*U+BKY*V)
      BKXC=BKX*C
      BKYC=BKY*C
      GM1U=GM1*U
      GM1V=GM1*V
      EV1=EV(1)
      T2(1,1)=EV1*(1.-FI2*RC2)
      T2(1,2)=EV1*GM1*RC2*U
      T2(1,3)=EV1*GM1*RC2*V
      T2(1,4)=-EV1*GM1*RC2
      EV2=EV1
      T2(2,1)=-EV2*RD*(BKY*U-BKX*V)
      T2(2,2)=EV2*BKY*RD
      T2(2,3)=-EV2*BKX*RD
      T2(2,4)=0.
      EV3=EV(2)
      T2(3,1)=EV3*BT*(FI2-CBST)
      T2(3,2)=EV3*BT*(BKXC-GM1U)
      T2(3,3)=EV3*BT*(BKYC-GM1V)
      T2(3,4)=EV3*BT*GM1
      EV4=EV(3)
      T2(4,1)=EV4*BT*(FI2+CBST)
      T2(4,2)=-EV4*BT*(BKXC+GM1U)
      T2(4,3)=-EV4*BT*(BKYC+GM1V)
      T2(4,4)=EV4*BT*GM1
      
c-----------------T-----------------------------------------------------

      AF=D/(S2*C)
      FC2RG1=(FI2+C2)*RGM1
      T1(1,1)=1.
      T1(1,2)=0.
      T1(1,3)=AF
      T1(1,4)=AF
      T1(2,1)=U
      T1(2,2)=BKY*D
      T1(2,3)=AF*(U+BKXC)
      T1(2,4)=AF*(U-BKXC)
      T1(3,1)=V
      T1(3,2)=-BKX*D
      T1(3,3)=AF*(V+BKYC)
      T1(3,4)=AF*(V-BKYC)
      T1(4,1)=FI2*RGM1
      T1(4,2)=D*(BKY*U-BKX*V)
      T1(4,3)=AF*(FC2RG1+CBST)
      T1(4,4)=AF*(FC2RG1-CBST)
      
c-----------------AB=T*EV*TINV------------------------------------------

      DO N=1,4
        DO M=1,4
          AB(N,M)=T1(N,1)*T2(1,M)+T1(N,2)*T2(2,M)
     &           +T1(N,3)*T2(3,M)+T1(N,4)*T2(4,M)
        enddo
      enddo

c-----------------------------------------------------------------------

      RETURN
      END