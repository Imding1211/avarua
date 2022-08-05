      subroutine solver(ITMAX,isave,ini_t,ini_time,max_thrust)
!-----------------------------------------------------------------------
!     a  2-d multi-block for compressible flow code  6/20,2019 
!
!     6 blocks
!
!     with  6 type BCs
!
!     1 --x inlet (subsonic)
!     2 --y outlet
!     3  lower wall
!     4  upper wall
!     5  left wall
!     6  right wall
!     7  lower outlet (subsonic)
!     8  upper outlet (subsonic)
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
      !INCLUDE 'omp_lib.h'
      INTEGER :: n_threads = 6

      INTEGER, intent(in) :: ITMAX
      INTEGER, intent(in) :: isave

      REAL, intent(in) :: ini_t
!     REAL, intent(in) :: ini_p
      REAL, intent(in) :: ini_time
      REAL, intent(out) :: max_thrust

!-----------------------------------------------------------------------

      REAL*8 :: RMACH,CFL,EPS,CK,HI,W
      REAL*8 :: XIXCF,XIYCF,ETXCF,ETYCF
      REAL*8 :: UF,VF,PF,SONF,UINF,VINF
      REAL*8 :: TEMPERATURE0,TEMPERATURE,TEMPER
      REAL*8 :: HALF,ONE,TWO,ZERO,TEN,HH
      REAL*8 :: DDD,R,PAI,RGM,RGM1,RGGM1,RG2M1,GM12,GM1R4,GRGM1
      REAL*8 :: CSP,CPG,C13,C43,PRL,PRT,FI,SLT,VEO
      REAL*8 :: XT,YT,X,Y,Z,CASEE,Q1Z,PZ,PZ2,Z2
      REAL*8 :: PINF,QT,Q1I,Q2I,Q3I,Q4I,Q5I
      REAL*8 :: VOL,RJ,RJJ,XXIX,XXIY,XETX,XETY
      REAL*8 :: TOT,DTMIN,DTMAX,DT
      REAL*8 :: XIXC,XIYC,ETXC,ETYC
      REAL*8 :: DTI,DMIN1
      REAL*8 :: Q1A,Q2A,Q3A,Q4A,Q5A
      REAL*8 :: Q1B,QBF,F,S,GM,QD,Q,UP,VP,DM,DS
      REAL*8 :: Q1E,Q2E,Q3E,Q4E,Q5E,RL,RM,DP
      REAL*8 :: UA,VA,UE,VE,QA,PA,CA,QE,PE,CE,UN,UM
      REAL*8 :: RKX,RKY,RSXY,ABK,PB,AKK,UB,VB,ZF
      REAL*8 :: DA,DB,ZZ,BKX,BKY,UR,VR,C0,DR,PR
      REAL*8 :: TERM,PSEC,ZFF,DFF,UFF,VFF,PFF,P
      REAL*8 :: SOURCE,SOURCE2
      REAL*8 :: QOLD1DS,QOLD1S,QOLD1,QOLD1D,QOLD2D,QOLD2SS,QOLDDSS,RESID
      REAL*8 :: VPR,VPL,HE1DT
      REAL*8 :: DDGMAX,DDX,DDY
      REAL*8 :: DDX1,DDX2,DG,DD
      REAL*8 :: UVEL,VVEL,PRE,AA
      REAL*8 :: DPHI,GM1,DEL,FMACH,CT0
      REAL*8 :: VISL,VIST,UB1,UB2,UB3,UB4
      REAL*8 :: XIX,XIY,ETX,ETY,AFB,BFB,A,B,BV
      REAL*8 :: UVR,UVL,DU,UV,DMD,DAST,DASS,VPR1,VPL1,VNR,VNL,VNR1,VNL1
      REAL*8 :: QM,TP,EINF,Q0,UUR,UUL,UUU1
      REAL*8 :: alltime,thrust,temper_i

      INTEGER :: I,J,K,N,M,IST,IAUSM,IFD,IO,NB,MB,JK,JP1,KP1
      INTEGER :: NJMM,NKMM,JMAX,KMAX,JMP1,JXM1,JXM2,KMP1,KXM1,KXM2
      INTEGER :: IIC,IT,IFREE,ISUBIN,IWALLL,ITEACHER
      INTEGER :: NKM3,NJM1,NJM5,NKM2,NJM6

      INTEGER :: LCH,LCHT
      INTEGER :: split, inlet_closed, outlet_closed,subsonic_outflow
      INTEGER :: inlet_frequence
      INTEGER :: rkkkk,splitting
      INTEGER :: jj,M1,M2,N1,N2
      INTEGER :: ITTT,ITTTTT,ITTTTTTTT

      character*80 :: filename

!-----------------------------------------------------------------------

      PARAMETER (RMACH=340.1,CFL=0.01,EPS=1.E-7,CK=0.33333,hi=1.)

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      COMMON/IJK/JMAX,JXM1,JMP1,KMAX,KXM1,KMP1

      COMMON/GMD1/DEL,UINF(9),VINF(9),PINF(9)

      COMMON/iau/iausm

      COMMON/GMD2/FMACH

      COMMON/SD/S(9,0:420,0:105)

      COMMON/XY/X(-1:420,-1:105),Y(-1:420,-1:105)

      COMMON/JACO/RJ(-1:420,-1:105)

      COMMON/QD/Q(9,-1:420,-1:105)

      COMMON/Qw/w(9,-1:420,-1:105)

      COMMON/PART1/UP(0:420,0:105),VP(0:420,0:105)

      COMMON/PART2/DM,DS,PAI,RL,RM,VEO,DP

      COMMON/VISC/CT0(9),C43,C13,PRL,PRT

      COMMON/VISC1/VISL(-1:420,-1:105),VIST(-1:420,-1:105)

      COMMON/CCC/UB1(-1:2401,-1:401,5),UB2(-1:2401,-1:401,5),UB3(-1:2401,-1:401,5),UB4(-1:2401,-1:401,5)

      COMMON/XEXY/XIX(0:420,0:105),XIY(0:420,0:105),ETX(0:420,0:105),ETY(0:420,0:105)

      DIMENSION DT(0:420,0:105),QBF(6,2),F(6,-1:1000),AFB(4,4),BFB(4,4),A(4,4),B(4),BV(4,4)
      DIMENSION uvr(10,-1:2401),uvl(10,-1:2401),du(10,-1:2401),uv(10,-1:2401),dmd(10,-1:2401)

      DIMENSION DAST(9,0:420,0:105),DASS(9,0:420,0:105)     
      DIMENSION vpr1(-1:2401,9),vpl1(-1:2401,9)
      DIMENSION vpr(-1:2401,9) ,vpl(-1:2401,9)
      DIMENSION vnr(-1:2401,9) ,vnl(-1:2401,9)
      DIMENSION vnr1(-1:2401,9) ,vnl1(-1:2401,9)
      DIMENSION un(-1:2401,9)

!------------------MULTIBLOCK VARIABLES---------------------------------

      DIMENSION QT(9,6,-1:420,-1:105),RJJ(9,-1:420,-1:105),QM(9,6,-1:420,-1:105),TP(9), DG(9,-1:420,-1:105)
      DIMENSION XXIX(9,0:420,0:105),XXIY(9,0:420,0:105),XETX(9,0:420,0:105),XETY(9,0:420,0:105)
      DIMENSION XT(9,-1:420,-1:105),YT(9,-1:420,-1:105),NJMM(9),NKMM(9)     
      DIMENSION QD(9,6,0:420,0:105)
      DIMENSION IFD(9,4),IO(9,4),LCH(9)
      DIMENSION EINF(9),Q1I(9),Q2I(9),Q3I(9),Q4I(9),Q5I(9)

!------------------CHEM STRANG------------------------------------------

      DIMENSION p(-1:501,-1:201)
      DIMENSION temperature(-1:501,-1:201)
      DIMENSION hh(-1:501,-1:201)
      DIMENSION source(-1:501,-1:201,5)
      DIMENSION source2(-1:501,-1:201,5)
      DIMENSION Qold1ds(-1:501,-1:201,5)
      DIMENSION Qold1s(-1:501,-1:201,5)
      DIMENSION Qold1d(-1:501,-1:201,5)   
      DIMENSION Qold1(-1:501,-1:201,5)
      DIMENSION Qold2d(-1:501,-1:201,5)
      DIMENSION Qold2ss(-1:501,-1:201,5)
      DIMENSION Qolddss(-1:501,-1:201,5)
      DIMENSION uuu1(-1:501,5),uur(-1:501,5),uul(-1:501,5)

!------------------STATEMENT FUNCTION-----------------------------------

      XIXCF(J,K)=.5*(XIX(J,K)+XIX(J+1,K))
      XIYCF(J,K)=.5*(XIY(J,K)+XIY(J+1,K))
      ETXCF(J,K)=.5*(ETX(J,K)+ETX(J,K+1))
      ETYCF(J,K)=.5*(ETY(J,K)+ETY(J,K+1))
      UF(J,K)=Q(2,J,K)/Q(1,J,K)
      VF(J,K)=Q(3,J,K)/Q(1,J,K)
      PF(J,K)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
      SONF(J,K)=DSQRT(DABS(GM*PF(J,K)/Q(1,J,K)))

!------------------PARAMETER SETUP--------------------------------------

      q0=0.5196d10
      temperature0=0.1155d10
      ddd=1.0/0.5825d10
      r=1.4

      half=5d-1
      one=10d-1
      two=one+one
      zero=0.0
      ten=1d1
      ist=1
    
      PAI=3.1415926535
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
      PRL=.72
      PRT=.9
      FI=1.
      SLT=1.
      VEO=0.
      max_thrust = 0.

      NB=6 

!------------------READ GRID DATA---------------------------------------

      OPEN(997,FILE='pdejet6.grd')
      OPEN(993,FILE='pdemb6.dat')

!------------------INTERFACE DATA---------------------------------------

      DO MB=1,nb

        READ(997,*)NJMM(MB),NKMM(MB)

        DO K=1,NKMM(MB)
          DO J=1,NJMM(MB)
            READ(997,*)XT(MB,J,K),YT(MB,J,K)
          ENDDO
        ENDDO

        READ(993,*)IO(MB,1),IO(MB,2),IO(MB,3),IO(MB,4)

      ENDDO

      close(997)
      close(993)

!-----------------------------------------------------------------------   

      DO I=1,nb
        LCH(I)=I
      ENDDO

      DO MB=1,nb
      
        JMAX=NJMM(MB)
        KMAX=NKMM(MB)
        JMP1=JMAX+1
        JXM1=JMAX-1
        JXM2=JXM1-1
        KMP1=KMAX+1
        KXM1=KMAX-1
        KXM2=KXM1-1
        JK=JXM1*KXM1

!------------------Q TERM-----------------------------------------------

!$OMP parallel do
        DO K=1,KMAX
          DO J=1,JMAX
            X(J,K)=XT(MB,J,K)
            Y(J,K)=YT(MB,J,K)
          ENDDO
        ENDDO
!$OMP end parallel do

!-----------------------------------------------------------------------

        splitting=1

!-----------------------------------------------------------------------

        q1z = 0.29 * 1.29d-3   ! dnesity
        pz  = 0.29 * 1.01325d6 ! 0.29atm
        z   = 0.0               ! z

        temper_i = ini_t

        UINF(MB)=0.0
        VINF(MB)=0.0
        PINF(MB)=pz

!$OMP parallel do
        DO  K=-1,KMP1
          DO  J=-1,JMP1
            QT(MB,1,J,K)=q1z
            QT(MB,2,J,K)=q1z*UINF(MB)
            QT(MB,3,J,K)=q1z*VINF(MB)
            QT(MB,5,J,K)=q1z*z
            QT(MB,4,J,K)=(temper_i*(QT(MB,1,J,K)*287d4))/GM1+q0*Qt(mb,5,J,K)+.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)
          enddo
        enddo
!$OMP end parallel do

!$OMP parallel do
        DO  K=-1,KMP1
          DO  J=-1,JMP1
            Q1I(MB)=QT(MB,1,J,K)
            Q2I(MB)=QT(MB,2,J,K)
            Q3I(MB)=QT(MB,3,J,K)
            Q4I(MB)=QT(MB,4,J,K)
            Q5I(MB)=QT(MB,5,J,K)
          enddo
        enddo
!$OMP end parallel do

!------------------ETX ETY----------------------------------------------

!$OMP parallel do
        DO K=1,KMAX
          DO J=1,JXM1
             JP1=J+1
             ETX(J,K)=-Y(JP1,K)+Y(J,K)
             ETY(J,K)=X(JP1,K)-X(J,K)
             ETX(0,K)=ETX(1,K)
             ETY(0,K)=ETY(1,K)
             ETX(JMAX,K)=ETX(JXM1,K)
             ETY(JMAX,K)=ETY(JXM1,K)
          ENDDO
        ENDDO
!$OMP end parallel do

!$OMP parallel do
        DO J=1,JMAX
           ETX(J,0)=ETX(J,1)
           ETY(J,0)=ETY(J,1)
           ETX(J,KMAX)=ETX(J,KXM1)
           ETY(J,KMAX)=ETY(J,KXM1)
        ENDDO
!$OMP end parallel do

!------------------XIX XIY----------------------------------------------

!$OMP parallel do
        DO K=1,KXM1
          DO J=1,JMAX
           KP1=K+1
           XIX(J,K)=Y(J,KP1)-Y(J,K)
           XIY(J,K)=-X(J,KP1)+X(J,K)
          ENDDO
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        DO J=1,JMAX
           XIX(J,0)=XIX(J,1)
           XIY(J,0)=XIY(J,1)
           XIX(J,KMAX)=XIX(J,KXM1)
           XIY(J,KMAX)=XIY(J,KXM1)
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        DO K=1,KMAX
           XIX(0,K)=XIX(1,K)
           XIY(0,K)=XIY(1,K)
           XIX(JMAX,K)=XIX(JXM1,K)
           XIY(JMAX,K)=XIY(JXM1,K)
        ENDDO
!$OMP end parallel do

!------------------RJ=1./VOL--------------------------------------------

!$OMP parallel do
        DO K=1,KXM1
          DO J=1,JXM1
            KP1=K+1
            JP1=J+1
            VOL=((X(J,K)-X(JP1,KP1))*(Y(JP1,K)-Y(J,KP1))-(X(JP1,K)-X(J,KP1))*(Y(J,K)-Y(JP1,KP1)))/2.
            RJ(J,K)=1./VOL
          ENDDO
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        DO J=1,JXM1
           RJ(J,0)=RJ(J,1)
           RJ(J,-1)=RJ(J,0)
           RJ(J,KMAX)=RJ(J,KXM1)
           RJ(J,KMP1)=RJ(J,KMAX)
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        DO K=1,KXM1
           RJ(0,K)=RJ(1,K)
           RJ(-1,K)=RJ(0,K)
           RJ(JMAX,K)=RJ(JXM1,K)
           RJ(JMP1,K)=RJ(JMAX,K)
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        DO K=-1,KMP1
          DO J=-1,JMP1
            RJJ(MB,J,K)=RJ(J,K)
          ENDDO
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        DO K=0,KMAX
          DO J=0,JMAX
            XXIX(MB,J,K)=XIX(J,K)
            XXIY(MB,J,K)=XIY(J,K)
            XETX(MB,J,K)=ETX(J,K)
            XETY(MB,J,K)=ETY(J,K)
          ENDDO
        ENDDO
!$OMP end parallel do

      ENDDO

!-----------------------------------------------------------------------
!     TIME MARCHING
!-----------------------------------------------------------------------

      tot=0.
      dtmin =0. 
      iic = 0 

      DO IT=1,ITMAX

        tot=tot+dtmin
      
!------------------VARIABLES CHANGE TO VARIABLES------------------------

        DO MB=1,nb
          JMAX=NJMM(MB)
          KMAX=NKMM(MB)
          JMP1=JMAX+1
          JXM1=JMAX-1
          JXM2=JXM1-1
          KMP1=KMAX+1
          KXM1=KMAX-1
          KXM2=KXM1-1
          JK=JXM1*KXM1

!$OMP parallel do
          DO K=0,KMAX
            DO J=0,JMAX
              XIX(J,K)=XXIX(MB,J,K)
              XIY(J,K)=XXIY(MB,J,K)
              ETX(J,K)=XETX(MB,J,K)
              ETY(J,K)=XETY(MB,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do
!$OMP parallel do
          DO K=-1,KMP1
            DO J=-1,JMP1
              RJ(J,K)=RJJ(MB,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do
!$OMP parallel do
          DO K=-1,KMP1
            DO J=-1,JMP1
              Q(1,J,K)=QT(MB,1,J,K)
              Q(2,J,K)=QT(MB,2,J,K)
              Q(3,J,K)=QT(MB,3,J,K)
              Q(4,J,K)=QT(MB,4,J,K)
              Q(5,J,K)=QT(MB,5,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do
!$OMP parallel do     
          DO K=1,KMAX
            DO J=1,JMAX
              X(J,K)=XT(MB,J,K)
              Y(J,K)=YT(MB,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do

!------------------SPACE VARYING DT BASED ON VARIABLE CFL.--------------

          DTMIN=100.
          DTMAX=0.

!------------------SPACE VARYING DT BASED ON VARIABLE CFL.--------------

          if(ist.eq.0) then

!$OMP parallel do
            DO  K=1,KXM1
              DO  J=1,JXM1
                XIXC=XIXCF(J,K)
                XIYC=XIYCF(J,K)
                ETXC=ETXCF(J,K)
                ETYC=ETYCF(J,K)

                DT(j,k)=CFL/((DABS(XIXC*UF(J,K)+XIYC*VF(J,K))&
                       +DABS(ETXC*UF(J,K)+ETYC*VF(J,K))+SONF(J,K)&
                       *DSQRT(DABS(XIXC**2+XIYC**2+ETXC**2+ETYC**2)))&
                       *RJ(J,K))
               enddo
            enddo
!$OMP end parallel do

!-----------------------------------------------------------------------
        
          else

!$OMP parallel do            
            DO K=1,KXM1
              DO J=1,JXM1
                XIXC=XIXCF(J,K)
                XIYC=XIYCF(J,K)
                ETXC=ETXCF(J,K)
                ETYC=ETYCF(J,K)

                DTI=CFL/((DABS(XIXC*UF(J,K)+XIYC*VF(J,K))&
                   +DABS(ETXC*UF(J,K)+ETYC*VF(J,K))+SONF(J,K)&
                   *DSQRT(DABS(XIXC**2+XIYC**2+ETXC**2+ETYC**2)))&
                   *RJ(J,K))

                DTMIN=DMIN1(DTMIN,DTI)
              enddo
            enddo
!$OMP end parallel do

          endif
        
        ENDDO

!-----------------------------------------------------------------------

        DO MB=1,nb

!$OMP parallel do
          DO  K=1,KXM1
            DO  J=1,JXM1
              DT(J,K)=DTmin
            ENDDO
          ENDDO
!$OMP end parallel do

!------------------VARIABLES CHANGE TO VARIABLES------------------------

          JMAX=NJMM(MB)
          KMAX=NKMM(MB)
          JMP1=JMAX+1
          JXM1=JMAX-1
          JXM2=JXM1-1
          KMP1=KMAX+1
          KXM1=KMAX-1
          KXM2=KXM1-1
          JK=JXM1*KXM1

!$OMP parallel do
          DO K=0,KMAX
            DO J=0,JMAX
              XIX(J,K)=XXIX(MB,J,K)
              XIY(J,K)=XXIY(MB,J,K)
              ETX(J,K)=XETX(MB,J,K)
              ETY(J,K)=XETY(MB,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do
!$OMP parallel do
          DO K=-1,KMP1
            DO J=-1,JMP1
              RJ(J,K)=RJJ(MB,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do
!$OMP parallel do
          DO K=-1,KMP1
            DO J=-1,JMP1
              Q(1,J,K)=QT(MB,1,J,K)
              Q(2,J,K)=QT(MB,2,J,K)
              Q(3,J,K)=QT(MB,3,J,K)
              Q(4,J,K)=QT(MB,4,J,K)
              Q(5,J,K)=QT(MB,5,J,K)
            ENDDO
          ENDDO
!$OMP end parallel do

!-----------------------------------------------------------------------
!     FAR FIELD NUMERICAL B.!.
!-----------------------------------------------------------------------

!------------------CHECK THE OUTLET-------------------------------------

          DO I=1,4
            IF(IO(MB,I).EQ.2)THEN
!$OMP parallel do
              DO  K=1,KMAX
                M1=JXM1
                M2=JMAX
                N1=K
                N2=K
                Q1A=Q(1,M1,N1)
                Q2A=Q(2,M1,N1)
                Q3A=Q(3,M1,N1)
                Q4A=Q(4,M1,N1)
                Q5A=Q(5,M1,N1)

                Q1E=Q(1,M2,N2)
                Q2E=Q(2,M2,N2)
                Q3E=Q(3,M2,N2)
                Q4E=Q(4,M2,N2)
                Q5E=Q(5,M2,N2)

                UA=Q2A/Q1A
                VA=Q3A/Q1A
                UE=Q2E/Q1E
                VE=Q3E/Q1E

                QA=0.5*Q1A*(UA*UA+VA*VA)-q0*Q5A
                PA=.4*(Q4A-QA)
                CA=DSQRT(1.4*DABS(PA/Q1A))

                QE=0.5*Q1E*(UE*UE+VE*VE)-q0*Q5E
                PE=.4*(Q4E-QE)
                CE=DSQRT(1.4*DABS(PE/Q1E))

                RKX=XIX(M2,N2)
                RKY=XIY(M2,N2)
                ABK=DSQRT(RKX*RKX+RKY*RKY)

                PB=PINF(MB)

                Q1B=Q1A+(PB-PA)/(CE*CE)
                AKK=ABK*Q1E*CE

                UB=UA+RKX*(PA-PB)/AKK
                VB=VA+RKY*(PA-PB)/AKK

                ZF=Q5A/Q1A

                Q(1,JMAX,K)=Q1A
                Q(2,JMAX,K)=UA*Q1A
                Q(3,JMAX,K)=VA*Q1A
                Q(4,JMAX,K)=Pa*2.5+0.5*q1a*(UA*UA+VA*VA)+q0*q1a*ZF !test
                Q(5,JMAX,K)=Q1A*ZF

                Q(1,jmp1,K)=Q(1,JMAX,K)
                Q(2,jmp1,K)=Q(2,JMAX,K)
                Q(3,jmp1,K)=Q(3,JMAX,K)
                Q(4,jmp1,K)=Q(4,JMAX,K)
                Q(5,jmp1,K)=Q(5,JMAX,K)
              enddo
!$OMP end parallel do
            ENDIF
          ENDDO

!-----------------------------------------------------------------------
!     INLET BOUNDARY
!-----------------------------------------------------------------------

!------------------lower nozzle tube------------------------------------

          alltime = 0.00003

          DO I=1,4
            IF(IO(MB,I).EQ.1)THEN
              psec = mod(tot,alltime)
              if(Psec .le. (alltime*ini_time))then
                iic = iic +1

!------------------refill----------------------------------------------- 

                if( iic .gt. 0  .and.  iic .le. 1) then   !refill and fire
!$OMP parallel do
                  DO K=-1,KMP1
                    DO J=-1,JMP1-5
                      z2=1.0       
                      q1z=Q(1,J,K)
                      UA=Q(2,J,K)/Q(1,J,K)
                      VA=Q(3,J,K)/Q(1,J,K)
                      pz=PF(j,k)

                      Q(1,J,K)=q1z
                      Q(2,J,K)=q1z*UA
                      Q(3,J,K)=q1z*VA
                      Q(4,J,K)=pz/(r-1)+0.5*q1z*(UA**2+VA**2)+q0*q1z*z2
                      Q(5,J,K)=q1z*z2
                    enddo
                  enddo
!$OMP end parallel do

!------------------fire-------------------------------------------------

                  PB= 30.*1.01325d6
                  DB=1.29d-3        ! density

!$OMP parallel do
                  DO  J=5,10
                    DO  K=1,Kxm1
                      M1=j
                      N1=K

                      Q1A=Q(1,M1,N1)
                      Q2A=Q(2,M1,N1)
                      Q3A=Q(3,M1,N1)
                      ZF=Q(5,M1,N1)/Q1A

                      UA=Q2A/Q1A
                      VA=Q3A/Q1A

                      Q(4,J,K)=PB*2.5+0.5*DB*(UA*UA+VA*VA)+q0*DB*ZF
                    enddo
                  enddo
!$OMP end parallel do

                endif

!------------------wall-------------------------------------------------

!$OMP parallel do              
                DO  K=1,KXM1    
                  RKX=XIX(1,K)
                  RKY=XIY(1,K)
                  RSXY=1./DSQRT(RKX*RKX+RKY*RKY)
                  BKX=RKX*RSXY
                  BKY=RKY*RSXY
                  C0=SONF(0,K)
                  DR=Q(1,1,k)
                  PR=PF(1,k)
                  UR=UF(1,k)
                  VR=VF(1,k)
                  ZF=Q(5,1,k)/DR

                  TERM=BKX*UR+BKY*VR
                  UB=UR-BKX*TERM
                  VB=VR-BKY*TERM
                  PB=PR
                  DB=DR

                  Q(1,0,k)=DB
                  Q(2,0,k)=db*ub
                  Q(3,0,k)=db*vb
                  Q(4,0,k)=PB*RGM1+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
                  Q(5,0,K)=DB*ZF
                enddo
!$OMP end parallel do


!-----------------------------------------------------------------------

              elseif (Psec .gt. (alltime*ini_time))then
                iic = 0
                pB=2.2*1.01325d6  ! wind
                DB=1.29d-3        ! density

!$OMP parallel do                
                DO  K=1,KMAX
                  M1=1
                  M2=0
                  N1=K
                  N2=K
                  Q1A=Q(1,M1,N1)
                  Q2A=Q(2,M1,N1)
                  Q3A=Q(3,M1,N1)
                  Q4A=Q(4,M1,N1)
                  Q5A=Q(5,M1,N1)

                  da =q1a
                  UA=Q2A/Q1A
                  VA=Q3A/Q1A
                  QA=0.5*Q1A*(UA*UA+VA*VA)
                  PA=0.4*(Q4A-QA)
                  CA=DSQRT(1.4*DABS(PA/Q1A))
                  RKX=XIX(1,K)
                  RKY=XIY(1,K)
                  RSXY=1./DSQRT(RKX*RKX+RKY*RKY)
                  BKX=RKX*RSXY
                  BKY=RKY*RSXY
                  UB=UR-BKX*TERM
                  VB=VR-BKY*TERM
                  ZF=Q5A/DB

                  Q(1,0,K)= db
                  Q(2,0,K)=UA*db
                  Q(3,0,K)=VA*db
                  Q(4,0,K)=PB*2.5+0.5*da*(UA*UA+VA*VA)+q0*da*ZF
                  Q(5,0,K)=DB*ZF
                enddo
!$OMP end parallel do

              endif
            ENDIF
          ENDDO

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SOLID SURFACE
!-----------------------------------------------------------------------

!------------------BOTTOM BOUNDARY--------------------------------------

        DO I=1,4
          IF(IO(MB,I).EQ.3)THEN

!$OMP parallel do
            DO J=1,JXM1
              RKX=ETX(J,1)
              RKY=ETY(J,1)
              RSXY=1./DSQRT(DABS(RKX*RKX+RKY*RKY))
              BKX=RKX*RSXY
              BKY=RKY*RSXY
              C0=SONF(J,0)
              DR=Q(1,J,1)
              PR=PF(J,1)
              UR=UF(J,1)
              VR=VF(J,1)
              TERM=BKX*UR+BKY*VR
              PB=PR-Q(1,J,0)*C0*TERM
              ZF=Q(5,J,1)/DR

              DB=DR+(PB-PR)/C0**2
              UB=UR-BKX*TERM
              VB=VR-BKY*TERM
              Q(1,J,0)=DB
              Q(2,J,0)=DB*UB
              Q(3,J,0)=DB*VB
              Q(4,J,0)=PB*RGM1+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
              Q(5,J,0)=DB*ZF
            enddo
!$OMP end parallel do

          ENDIF
        ENDDO

!------------------333-open--------------------------------------------

        DO I=1,4
          IF(IO(MB,I).EQ.333)THEN

!$OMP parallel do
            do j=1,jmax
              M1=J
              M2=J
              N1=1
              N2=0

              Q1A=Q(1,M1,N1)
              Q2A=Q(2,M1,N1)
              Q3A=Q(3,M1,N1)
              Q4A=Q(4,M1,N1)
              Q5A=Q(5,M1,N1)

              Q1E=Q(1,M2,N2)
              Q2E=Q(2,M2,N2)
              Q3E=Q(3,M2,N2)
              Q4E=Q(4,M2,N2)
              Q5E=Q(5,M2,N2)

              UA=Q2A/Q1A
              VA=Q3A/Q1A
              UE=Q2E/Q1E
              VE=Q3E/Q1E

              QA=0.5*Q1A*(UA*UA+VA*VA)-q0*Q5A
              PA=.4*(Q4A-QA)
              CA=DSQRT(1.4*DABS(PA/Q1A))

              QE=0.5*Q1E*(UE*UE+VE*VE)-q0*Q5E
              PE=.4*(Q4E-QE)
              CE=DSQRT(1.4*DABS(PE/Q1E))

              RKX=ETX(M2,N2)
              RKY=ETY(M2,N2)
              ABK=DSQRT(RKX*RKX+RKY*RKY)

              PB=PINF(MB)

              Q1B=Q1A+(PB-PA)/(CE*CE)
              AKK=ABK*Q1E*CE

              UB=UA+RKX*(PA-PB)/AKK
              VB=VA+RKY*(PA-PB)/AKK

              ZF=Q5A/Q1A

              Q(1,J,0)=Q1A
              Q(2,J,0)=UA*Q1A
              Q(3,J,0)=VA*Q1A
              Q(4,J,0)=Pa*2.5+0.5*q1a*(UA*UA+VA*VA)+q0*q1a*ZF
              Q(5,J,0)=Q1A*ZF

              Q(1,J,-1)=Q(1,J,0)
              Q(2,J,-1)=Q(2,J,0)
              Q(3,J,-1)=Q(3,J,0)
              Q(4,J,-1)=Q(4,J,0)
              Q(5,J,-1)=Q(5,J,0)
            enddo
!$OMP end parallel do

          ENDIF
        ENDDO

!------------------TOP BOUNDARY-----------------------------------------

        DO I=1,4
          IF(IO(MB,I).EQ.4)THEN

!$OMP parallel do
            DO J=1,JXM1
              RKX=ETX(J,KXM2)
              RKY=ETY(J,KXM2)
              RSXY=1./DSQRT(RKX*RKX+RKY*RKY)
              BKX=RKX*RSXY
              BKY=RKY*RSXY
              C0=SONF(J,KXm1)
              DR=Q(1,J,KXM2)
              PR=PF(J,KXM2)
              UR=UF(J,KXM2)
              VR=VF(J,KXM2)
              TERM=BKX*UR+BKY*VR
              PB=PR+Q(1,J,KXm1)*C0*TERM
              ZF=Q(5,J,KXM2)/DR
              
              DB=DR+(PB-PR)/C0**2
              UB=UR-BKX*TERM
              VB=VR-BKY*TERM
              Q(1,J,KXm1)=DB
              Q(2,J,kxm1)=DB*UB
              Q(3,J,kxm1)=DB*VB
              Q(4,J,kxm1)=PB*RGM1+.5*DB*(UB*UB+VB*VB)+q0*DB*ZF
              Q(5,J,kxm1)=DB*ZF
            enddo
!$OMP end parallel do

          ENDIF
        ENDDO

!------------------444 open---------------------------------------------

        DO I=1,4
          IF(IO(MB,I).EQ.444)THEN

!$OMP parallel do
            DO J=1,Jmax
              M1=J
              M2=J
              N1=KXM2
              N2=kxm1
              Q1A=Q(1,M1,N1)
              Q2A=Q(2,M1,N1)
              Q3A=Q(3,M1,N1)
              Q4A=Q(4,M1,N1)
              Q5A=Q(5,M1,N1)

              Q1E=Q(1,M2,N2)
              Q2E=Q(2,M2,N2)
              Q3E=Q(3,M2,N2)
              Q4E=Q(4,M2,N2)
              Q5E=Q(5,M2,N2)

              UA=Q2A/Q1A
              VA=Q3A/Q1A
              UE=Q2E/Q1E
              VE=Q3E/Q1E

              QA=0.5*Q1A*(UA*UA+VA*VA)-q0*Q5A
              PA=.4*(Q4A-QA)
              CA=DSQRT(1.4*DABS(PA/Q1A))

              QE=0.5*Q1E*(UE*UE+VE*VE)-q0*Q5E
              PE=.4*(Q4E-QE)
              CE=DSQRT(1.4*DABS(PE/Q1E))

              RKX=ETX(M2,N2)
              RKY=ETY(M2,N2)
              ABK=DSQRT(RKX*RKX+RKY*RKY)

              PB=PINF(MB)

              Q1B=Q1A+(PB-PA)/(CE*CE)
              AKK=ABK*Q1E*CE

              UB=UA+RKX*(PA-PB)/AKK
              VB=VA+RKY*(PA-PB)/AKK

              ZF=Q5A/Q1A

              Q(1,J,kxm1)=Q1A
              Q(2,J,kxm1)=UA*Q1A
              Q(3,J,kxm1)=VA*Q1A
              Q(4,J,kxm1)=Pa*2.5+0.5*q1a*(UA*UA+VA*VA)+q0*q1a*ZF
              Q(5,J,kxm1)=Q1A*ZF

              Q(1,J,kmax)=Q(1,J,kxm1)
              Q(2,J,kmax)=Q(2,J,kxm1)
              Q(3,J,kmax)=Q(3,J,kxm1)
              Q(4,J,kmax)=Q(4,J,kxm1)
              Q(5,J,kmax)=Q(5,J,kxm1)
            enddo
!$OMP end parallel do

          ENDIF
        ENDDO

!-----------------------------------------------------------------------

        DO I=1,4
          IF(IO(MB,I).EQ.555)THEN

!$OMP parallel do
            do k=1,kmax
              M1=1
              M2=0
              N1=K
              N2=K
              Q1A=Q(1,M1,N1)
              Q2A=Q(2,M1,N1)
              Q3A=Q(3,M1,N1)
              Q4A=Q(4,M1,N1)
              Q5A=Q(5,M1,N1)

              Q1E=Q(1,M2,N2)
              Q2E=Q(2,M2,N2)
              Q3E=Q(3,M2,N2)
              Q4E=Q(4,M2,N2)
              Q5E=Q(5,M2,N2)

              UA=Q2A/Q1A
              VA=Q3A/Q1A
              UE=Q2E/Q1E
              VE=Q3E/Q1E

              QA=0.5*Q1A*(UA*UA+VA*VA)-q0*Q5A
              PA=.4*(Q4A-QA)
              CA=DSQRT(1.4*DABS(PA/Q1A))

              QE=0.5*Q1E*(UE*UE+VE*VE)-q0*Q5E
              PE=.4*(Q4E-QE)
              CE=DSQRT(1.4*DABS(PE/Q1E))

              RKX=XIX(M2,N2)
              RKY=XIY(M2,N2)
              ABK=DSQRT(RKX*RKX+RKY*RKY)

              PB=PINF(MB)

              Q1B=Q1A+(PB-PA)/(CE*CE)
              AKK=ABK*Q1E*CE

              UB=UA+RKX*(PA-PB)/AKK
              VB=VA+RKY*(PA-PB)/AKK

              ZF=Q5A/Q1A

              Q(1,JMAX,K)=Q1A
              Q(2,JMAX,K)=UA*Q1A
              Q(3,JMAX,K)=VA*Q1A
              Q(4,JMAX,K)=Pa*2.5+0.5*q1a*(UA*UA+VA*VA)+q0*q1a*ZF !test
              Q(5,JMAX,K)=Q1A*ZF

              Q(1,jmp1,K)=Q(1,JMAX,K)
              Q(2,jmp1,K)=Q(2,JMAX,K)
              Q(3,jmp1,K)=Q(3,JMAX,K)
              Q(4,jmp1,K)=Q(4,JMAX,K)
              Q(5,jmp1,K)=Q(5,JMAX,K)
            enddo
!$OMP end parallel do

          ENDIF 
        ENDDO

!-----------------------------------------------------------------------

!$OMP parallel do
        DO N=1,5
          DO M=1,JXM1
            Q(N,M,-1)=Q(N,M,0)
            Q(N,M,KMax)=Q(N,M,Kxm1)
            Q(N,M,KMP1)=Q(N,M,KMAX)
          ENDDO
        ENDDO
!$OMP end parallel do
!$OMP parallel do
        do n=1,5
          do m=1,kxm1
            Q(N,-1,M)=Q(N,0,M)
            Q(N,JMP1,M)=Q(N,JMAX,M)
          ENDDO
        ENDDO
!$OMP end parallel do

!-----------------------------------------------------------------------
!     MODIFIED THE INLET INTERFACE
!-----------------------------------------------------------------------

!------------------1----------------------------------------------------

        IF(MB.EQ.1)THEN

          ! 1 to 2
!$OMP parallel do
          DO N=1,5
            DO K=1,kxm1
              Q(N,JMAX,K)=QT(2,N,1,K)
              Q(N,JMP1,K)=QT(2,N,2,K)
            enddo
          enddo 
!$OMP end parallel do

!------------------2----------------------------------------------------

        ELSEIF(MB.EQ.2)THEN

          ! 2 to 4
!$OMP parallel do
          DO N=1,5
            DO J=1,JMAX
              Q(N,J,KMAX)=QT(4,N,J,1)
              Q(N,J,KMP1)=QT(4,N,J,2)
            enddo
          enddo 
!$OMP end parallel do
          
          ! 2 to 3  
!$OMP parallel do
          DO N=1,5
            DO J=1,JMAX                    
              Q(N,J,0)=QT(3,N,J,NKMM(3)-1)
              Q(N,J,-1)=QT(3,N,J,NKMM(3)-1-1)
            ENDDO
          ENDDO
!$OMP end parallel do

          ! 2 get 1
!$OMP parallel do          
          DO N=1,5
            DO K=1,KMAX
              Q(N,0,K)=QT(1,N,NJMM(1)-1,K) !1
              Q(N,-1,K)=QT(1,N,NJMM(1)-1-1,K)
            ENDDO
          ENDDO
!$OMP end parallel do

!------------------3----------------------------------------------------

        ELSEIF(MB.EQ.3)THEN

          ! 3 get 2
!$OMP parallel do
          DO N=1,5
            DO J=1,JMAX
              Q(N,J,KMAX)=QT(2,N,J,1)
              Q(N,J,KMP1)=QT(2,N,J,2)
            ENDDO
          ENDDO
!$OMP end parallel do

          ! 3 to 5 
!$OMP parallel do
          DO N=1,5
            DO K=1,KMAX
              Q(N,0,K)=QT(5,N,NJMM(5)-1,K)
              Q(N,-1,K)=QT(5,N,NJMM(5)-1-1,K)
            ENDDO
          ENDDO
!$OMP end parallel do

!------------------4----------------------------------------------------

        ELSEIF(MB.EQ.4)THEN

          ! 4 get 2
!$OMP parallel do
          DO N=1,5
            DO J=1,JMAX
              Q(N,J,0)=QT(2,N,J,NKMM(2)-1)
              Q(N,J,-1)=QT(2,N,J,NKMM(2)-1-1)
            ENDDO
          ENDDO
!$OMP end parallel do

          ! 4 to 6
!$OMP parallel do          
          DO N=1,5      
            DO K=1,KMAX
              Q(N,0,K)=QT(6,N,NJMM(6)-1,K)
              Q(N,-1,K)=QT(6,N,NJMM(6)-1-1,K)
            END DO
          ENDDO
!$OMP end parallel do

!------------------5----------------------------------------------------

        ELSEIF(MB.EQ.5)THEN

          ! 5 get 3
!$OMP parallel do
          DO N=1,5
            DO K=1,Kxm1
              Q(N,JMAX,K)=QT(3,N,1,K)
              Q(N,JMP1,K)=QT(3,N,2,K)
            enddo
          enddo
!$OMP end parallel do

!------------------6----------------------------------------------------

        ELSEIF(MB.EQ.6)THEN

          ! 6 get 4
!$OMP parallel do
          DO N=1,5
            DO K=1,KMAX
              Q(N,JMAX,k)=QT(4,N,1,K)
              Q(N,JMP1,k)=QT(4,N,2,K)
            enddo
          enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

        ENDIF

!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!     control z between 0 and 1
!-----------------------------------------------------------------------

!$OMP parallel do
      DO K=-1,KMP1  
        DO J=-1,JMP1
          ZFF=Q(5,J,K)/Q(1,J,K)
          ZFF= min(1.0,ZFF)  
          ZFF=max(0.0,ZFF)
          DFF=Q(1,J,K)
          UFF=Q(2,J,K)/Q(1,J,K)
          VFF=Q(3,J,K)/Q(1,J,K)
          PFF=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
          Q(4,J,K)=PFF/(r-1)+0.5*DFF*(UFF**2+VFF**2)+q0*DFF*ZFF
          Q(5,J,K)=DFF*ZFF
        enddo
      enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

      if(splitting .EQ. 1)then

        do split=1,4

          if(split .EQ. 1)then

!------------------source ONE S(k/2)------------------------------------

!$OMP parallel do
            DO K=-1,KMP1
              DO J=-1,JMP1
                p(j,k)=GM1*(Q(4,J,K)-0.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            DO K=-1,KMP1
              DO J=-1,JMP1
                temperature(j,k)=p(j,k)/Q(1,J,K)
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            DO K=-1,KMP1
              DO J=-1,JMP1
                if ((temperature(j,k)-temperature0)>zero)then
                  hh(j,k)=1.0
                else
                  hh(j,k)=0.0  
                endif
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            DO K=-1,KMP1
              DO J=-1,JMP1
                source(j,k,1)=0.0
                source(j,k,2)=0.0
                source(j,k,3)=0.0
                source(j,k,4)=0.0
                source(j,k,5)=-(one/ddd)*Q(5,J,K)*hh(j,k)
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            DO K=-1,KMP1
              DO J=-1,JMP1
                source2(j,k,1)=0.0
                source2(j,k,2)=0.0
                source2(j,k,3)=0.0
                source2(j,k,4)=0.0
                source2(j,k,5)=-one*hh(j,k)
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            do n=1,5
              DO K=-1,KMP1
                DO J=-1,JMP1
                  Qold1ds(j,k,n)=half*DTMIN*source(j,k,n)/(one-half*half*DTMIN*source2(j,k,n) )
                  Qold1s(j,k,n)=Q(n,J,K)+Qold1ds(j,k,n)
                enddo
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            do n=1,5
              DO K=-1,KMP1
                DO J=-1,JMP1
                  Q(n,J,k)=Qold1s(j,k,n)
                enddo
              enddo
            enddo
!$OMP end parallel do

!------------------end S(k/2)-------------------------------------------

          elseif(split .EQ. 2)then

!------------------MAIN ONE SF(k)---------------------------------------

!-----------------------------------------------------------------------
!     EXPLICIT--1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  EXPLICIT  up
!-----------------------------------------------------------------------

!$OMP parallel do
            do k=-1,kmp1
              do j=-1,jmp1
                w(1,j,k)=Q(1,j,k)
                w(2,j,k)=Q(2,J,K)/Q(1,J,K)
                w(3,j,k)=Q(3,J,K)/Q(1,J,K)
                w(4,j,k)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
                w(5,j,k)=Q(5,J,K)/Q(1,J,K)
              enddo
            enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

            DO K=1,KXM1  
              RKX=XIX(1,K)
              RKY=XIY(1,K)
              DO N=1,5
                QBF(N,1)=w(N,0,K)
                QBF(N,2)=w(N,1,K)
              enddo

              CALL fvsw(1,RKX,RKY,QBF,F)

              DO  N=1,5
                DO J = -1,Jmp1
                  un(j,n) = w(N,j,K)
                enddo
              enddo

              call hmsth(un,vpr,vpl,jmax)

              DO J=1,JXM2
                JP1=J+1
                RKX=XIX(JP1,K)
                RKY=XIY(JP1,K)
                DO  N=1,5
                  QBF(N,1)= vpl(j,n)
                  QBF(N,2)= vpr(j,n)
                enddo
                CALL fVSw(JP1,RKX,RKY,QBF,F)
              enddo

              RKX=XIX(JMAX,K)
              RKY=XIY(JMAX,K)
              DO N=1,5
                QBF(N,1)=w(N,JXM1,K)
                QBF(N,2)=w(N,JMAX,K)
              enddo

              CALL FVSw(JMAX,RKX,RKY,QBF,F)

              DO J=1,JXM1
                HE1DT=DT(J,K)*RJ(J,K)
                DO N=1,5
                  S(N,J,K)=-HE1DT*(F(N,J+1)-F(N,J))
                enddo
              enddo
            enddo

!-----------------------------------------------------------------------

            DO J=1,JXM1
              RKX=ETX(J,1)
              RKY=ETY(J,1)
              DO N=1,5
                QBF(N,1)=w(N,J,0)
                QBF(N,2)=w(N,J,1)
              enddo

              CALL fVSw(1,RKX,RKY,QBF,F)

              DO N=1,5
                DO k = -1,kmp1
                  un(k,n) = w(N,j,K)
                enddo
              enddo

              call hmsth(un,vpr,vpl,kmax)

              DO K=1,KXM2
                KP1=K+1
                RKX=ETX(J,KP1)
                RKY=ETY(J,KP1)
                DO N=1,5
                  QBF(N,1)= vpl(k,n)
                  QBF(N,2)= vpr(k,n)
                enddo
                CALL fVSw(KP1,RKX,RKY,QBF,F)
              enddo

              RKX=ETX(J,KMAX)
              RKY=ETY(J,KMAX)

              DO N=1,5
                QBF(N,1)=w(N,J,KXM1)
                QBF(N,2)=w(N,J,KMAX)
              enddo
            
              CALL fVSw(KMAX,RKX,RKY,QBF,F)

              DO K=1,KXM1
                HE1DT=DT(J,K)*RJ(J,K)
                DO N=1,5
                  KP1=K+1
                  S(N,J,K)=S(N,J,K)-HE1DT*(F(N,KP1)-F(N,K))
                enddo
              enddo
            enddo

!-----------------------------------------------------------------------
!     EXPLICIT  dowm
!-----------------------------------------------------------------------

!$OMP parallel do
              do n=1,5
                DO K=-1,KMP1
                  DO J=-1,JMP1 
                    Qold1d(j,k,n)=S(N,J,K)
                    Qold1(j,k,n)=Qold1s(j,k,n)+Qold1d(j,k,n)
                  enddo
                enddo
              enddo
!$OMP end parallel do
!$OMP parallel do
              do n=1,5
                DO K=-1,KMP1
                  DO J=-1,JMP1
                    Q(n,J,k)=Qold1(j,k,n)
                  enddo
                enddo
              enddo
!$OMP end parallel do

!------------------end S(k/2)------------------------------------------- 

          elseif(split .EQ. 3)then

!------------------MAIN TWO SF(k)---------------------------------------

!-----------------------------------------------------------------------
!     EXPLICIT--2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  EXPLICIT  up
!-----------------------------------------------------------------------

!$OMP parallel do
            do k=-1,kmp1
              do j=-1,jmp1
                w(1,j,k)=Q(1,j,k)
                w(2,j,k)=Q(2,J,K)/Q(1,J,K)
                w(3,j,k)=Q(3,J,K)/Q(1,J,K)
                w(4,j,k)=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
                w(5,j,k)=Q(5,J,K)/Q(1,J,K)
              enddo
            enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

            DO K=1,KXM1  
              RKX=XIX(1,K)
              RKY=XIY(1,K)
              DO N=1,5
                QBF(N,1)=w(N,0,K)
                QBF(N,2)=w(N,1,K)
              enddo

              CALL fvsw(1,RKX,RKY,QBF,F)

              DO  N=1,5
                DO J = -1,Jmp1
                  un(j,n) = w(N,j,K)
                enddo
              enddo

              call hmsth(un,vpr,vpl,jmax)

              DO J=1,JXM2
                JP1=J+1
                RKX=XIX(JP1,K)
                RKY=XIY(JP1,K)
                DO  N=1,5
                  QBF(N,1)= vpl(j,n)
                  QBF(N,2)= vpr(j,n)
                enddo
                CALL fVSw(JP1,RKX,RKY,QBF,F)
              enddo

              RKX=XIX(JMAX,K)
              RKY=XIY(JMAX,K)
              DO N=1,5
                QBF(N,1)=w(N,JXM1,K)
                QBF(N,2)=w(N,JMAX,K)
              enddo

              CALL FVSw(JMAX,RKX,RKY,QBF,F)

              DO J=1,JXM1
                HE1DT=DT(J,K)*RJ(J,K)
                DO N=1,5
                  S(N,J,K)=-HE1DT*(F(N,J+1)-F(N,J))
                enddo
              enddo
            enddo

!-----------------------------------------------------------------------

            DO J=1,JXM1
              RKX=ETX(J,1)
              RKY=ETY(J,1)
              DO N=1,5
                QBF(N,1)=w(N,J,0)
                QBF(N,2)=w(N,J,1)
              enddo

              CALL fVSw(1,RKX,RKY,QBF,F)

              DO  N=1,5
                DO k = -1,kmp1
                  un(k,n) = w(N,j,K)
                enddo
              enddo

              call hmsth(un,vpr,vpl,kmax)

              DO K=1,KXM2
                KP1=K+1
                RKX=ETX(J,KP1)
                RKY=ETY(J,KP1)
                DO  N=1,5
                  QBF(N,1)= vpl(k,n)
                  QBF(N,2)= vpr(k,n)
                enddo
                CALL fVSw(KP1,RKX,RKY,QBF,F)
              enddo

              RKX=ETX(J,KMAX)
              RKY=ETY(J,KMAX)
              DO N=1,5
                QBF(N,1)=w(N,J,KXM1)
                QBF(N,2)=w(N,J,KMAX)
              enddo

              CALL fVSw(KMAX,RKX,RKY,QBF,F)

              DO K=1,KXM1
              HE1DT=DT(J,K)*RJ(J,K)
                DO N=1,5
                  KP1=K+1
                  S(N,J,K)=S(N,J,K)-HE1DT*(F(N,KP1)-F(N,K))
                enddo
              enddo
            enddo

!-----------------------------------------------------------------------
!     EXPLICIT  dowm
!-----------------------------------------------------------------------

!$OMP parallel do
            DO N=1,5
              DO K=-1,KMP1
                DO J=-1,JMP1
                  Qold2d(J,K,N)=S(N,J,K)
                  Qold2ss(J,K,N)=Qold1s(J,K,N)+half*(Qold1d(J,K,N)+Qold2d(J,K,N))
                enddo
              enddo
            ENDDO
!$OMP end parallel do
!$OMP parallel do
            do  n=1,5
              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  Q(n,J,k)=Qold2ss(j,k,n)
                enddo
              enddo
            enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

          elseif(split .EQ. 4)then

!------------------source TWO S(k/2)------------------------------------

!$OMP parallel do
            do  n=1,5
              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  Qold2ss(j,k,n)=Q(n,J,k)
                enddo
              enddo
            enddo
!$OMP end parallel do
!$OMP parallel do
            DO  K=-1,KMP1
              DO  J=-1,JMP1
                p(J,K)=GM1*(Qold2ss(J,K,4)-half*(Qold2ss(J,K,2)**2+Qold2ss(J,K,3)**2)/Qold2ss(J,K,1)-q0*Qold2ss(J,K,5))
              enddo
            ENDDO
!$OMP end parallel do
!$OMP parallel do
            DO  K=-1,KMP1
              DO  J=-1,JMP1
                temperature(J,K)=p(J,K)/Qold2ss(J,K,1)
              enddo
            ENDDO
!$OMP end parallel do
!$OMP parallel do
            DO  K=-1,KMP1
              DO  J=-1,JMP1
                if ((temperature(J,K)-temperature0)>zero)then
                  hh(J,K)=one
                else
                  hh(J,K)=zero  
                endif
              enddo
            ENDDO
!$OMP end parallel do
!$OMP parallel do
            DO  K=-1,KMP1
              DO  J=-1,JMP1
                source(J,K,1)=0.0
                source(J,K,2)=0.0
                source(J,K,3)=0.0
                source(J,K,4)=0.0
                source(J,K,5)=-(one/ddd)*Qold2ss(J,K,5)*hh(J,K)
              enddo
            ENDDO
!$OMP end parallel do
!$OMP parallel do
            DO  K=-1,KMP1
              DO  J=-1,JMP1
                source2(J,K,1)=0.0
                source2(J,K,2)=0.0
                source2(J,K,3)=0.0
                source2(J,K,4)=0.0
                source2(J,K,5)=-one*hh(J,K)
              enddo
            ENDDO
!$OMP end parallel do
!$OMP parallel do
            do N=1,5
              DO  K=-1,KMP1
                DO  J=-1,JMP1
                  Qolddss(J,K,N)=(half*DTMIN*source(J,K,N))/(one-half*half*DTMIN*source2(J,K,N))
                  Q(N,J,K)=Qold2ss(J,K,N)+Qolddss(J,K,N)
                enddo
              enddo
            ENDDO
!$OMP end parallel do

!-----------------------------------------------------------------------
!     control z between 0 and 1
!-----------------------------------------------------------------------

!$OMP parallel do
            DO  K=-1,KMP1  
              DO  J=-1,JMP1
                ZFF=Q(5,J,K)/Q(1,J,K)
                ZFF=min(1.0,ZFF)  
                ZFF=max(0.0,ZFF)
                DFF=Q(1,J,K)
                UFF=Q(2,J,K)/Q(1,J,K)
                VFF=Q(3,J,K)/Q(1,J,K)
                PFF=GM1*(Q(4,J,K)-.5*(Q(2,J,K)**2+Q(3,J,K)**2)/Q(1,J,K)-q0*Q(5,J,K))
                Q(4,J,K)=PFF/(r-1)+0.5*DFF*(UFF**2+VFF**2)+q0*DFF*ZFF
                Q(5,J,K)=DFF*ZFF
              enddo
            enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

          endif   
        enddo  
      endif

!-----------------------------------------------------------------------
!     RESIDUAL
!-----------------------------------------------------------------------

        RESID=0.

!$OMP parallel do
        DO K=1,KXM1
          DO J=1,JXM1
            RESID=RESID+S(1,J,K)*S(1,J,K)
          enddo
        enddo
!$OMP end parallel do

        RESID=DSQRT(DABS(RESID/JK))

        IF (MOD(IT,isave).EQ.0.or.it.eq.1) THEN
          !WRITE (*,*) MB,IT,RESID,tot,dtmin, iic
        ENDIF

!-----------------------------------------------------------------------
!     STORAGE GEOMETRY DATA
!-----------------------------------------------------------------------

!$OMP parallel do
        DO K=-1,KMP1
          DO J=-1,JMP1
            QT(MB,1,J,K)=Q(1,J,K)
            QT(MB,2,J,K)=Q(2,J,K)
            QT(MB,3,J,K)=Q(3,J,K)
            QT(MB,4,J,K)=Q(4,J,K)
            QT(MB,5,J,K)=Q(5,J,K)
          enddo
        enddo
!$OMP end parallel do

      enddo
        
      
      IF (MOD(IT,isave).EQ.0.or.IT.EQ.1.) THEN
        ittt=it/isave+1000
        ittttt=it/isave+3000
        itttttttt=it/isave+5000

        ddgmax= 0.

!$OMP parallel do
        DO MB=1,nb
          DO K=1,NKMM(MB)
            DO J=1,NJMM(MB)

              ddx = ( Qt(mb,1,j+1,k)-Qt(mb,1,j-1,k) )/2.
              ddy = ( Qt(mb,1,j,k+1)-Qt(mb,1,j,k-1) )/2.

              ddx1 = ( ddx * xxix(mb,j,k) + ddy * xetx(mb,j,k) )**2
              ddx2 = ( ddx * xxiy(mb,j,k) + ddy * xety(mb,j,k) )**2
              dg(MB,j,k) = dsqrt( ddx1 + ddx2) 
              ddgmax = Dmax1( ddgmax,dg(MB,j,k) )

            enddo
          enddo
        enddo
!$OMP end parallel do

!-----------------------------------------------------------------------

        DO MB=1,nb
          write (filename,'("flow-fluid",I5,".dat")')ittt
          open  (unit=ittt,file=filename,status='unknown')
          WRITE (ittt,*) 'Tittle="MULTIBLOCK"'
          WRITE (ittt,*) 'Variables = "X","Y","u","v","M","d","p","s","si","Z"'
          WRITE (ittt,*) 'zone I=',NJMM(MB),',J=',NKMM(MB),',f=point'
          WRITE (ittt,*) 'SOLUTIONTIME=',tot

!!!!$OMP parallel do
          DO K=1,NKMM(MB)
            DO J=1,NJMM(MB)

              dd=qt(mb,1,j,k)
              uvel=qt(mb,2,j,k)/qt(mb,1,j,k)
              vvel=qt(mb,3,j,k)/qt(mb,1,j,k)
              Pre=GM1*(Qt(mb,4,J,K)-.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)-q0*Qt(mb,5,J,K))

              aa=dsqrt( dabs(pre*gm/dd))
              uM=dsqrt(uvel**2+vvel**2)/aa

              dphi = dexp( -15d0 * dg(MB,j,k)/ddgmax)

              WRITE (ittt,*) XT(MB,J,K),YT(MB,J,K),uvel,vvel,um,dd,pre,dphi,dg(mb,j,k),qt(mb,5,j,k)/qt(mb,1,j,k)
            enddo
          enddo
!!!!$OMP end parallel do
        enddo
        close(ittt)

!------------------centrol line 3000------------------------------------

!!!!$OMP parallel do
        write (filename,'("centrol-line-3000",I5,".dat")')ittttt
        open  (unit=ittttt,file=filename,status='unknown')
        do mb =1,2
          DO  J=1,NJMM(MB)
            k=50
            dd=qt(mb,1,j,k)
            uvel=qt(mb,2,j,k)/qt(mb,1,j,k)
            vvel=qt(mb,3,j,k)/qt(mb,1,j,k)
            Pre=GM1*(Qt(mb,4,J,K)-.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)-q0*Qt(mb,5,J,K))
            aa=dsqrt( dabs(pre*gm/dd))
            uM=dsqrt(uvel**2+vvel**2)/aa
            WRITE (ittttt,*) XT(MB,J,K), uvel,uM
          enddo
        enddo
        close(ittttt)
!!!!$OMP end parallel do

!------------------centrol line 5000------------------------------------

!!!!$OMP parallel do
        write (filename,'("centrol-line-5000",I5,".dat")')itttttttt
        open  (unit=itttttttt,file=filename,status='unknown')
        do mb =1,2
          DO  J=1,NJMM(MB)
            k=50
            dd=qt(mb,1,j,k)
            uvel=qt(mb,2,j,k)/qt(mb,1,j,k)
            vvel=qt(mb,3,j,k)/qt(mb,1,j,k)
            Pre=GM1*(Qt(mb,4,J,K)-.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)-q0*Qt(mb,5,J,K))
            aa=dsqrt( dabs(pre*gm/dd))
            uM=dsqrt(uvel**2+vvel**2)/aa
            temper=Pre/(dd*287d4)
            WRITE (itttttttt,*) XT(MB,J,K),pre,temper
          enddo
        enddo
        close(itttttttt)
!!!!$OMP end parallel do

      endif

!------------------thrust line 7000-------------------------------------

      thrust=0.0

!$OMP parallel do
      do mb =1,1
        DO  K=1,NKMM(MB)
          J=NJMM(MB)
          dd=qt(mb,1,j,k)
          uvel=qt(mb,2,j,k)/qt(mb,1,j,k)
          vvel=qt(mb,3,j,k)/qt(mb,1,j,k)
          Pre=GM1*(Qt(mb,4,J,K)-.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)-q0*Qt(mb,5,J,K))
          thrust=thrust+dd*uvel*0.2*uvel+(Pre-0.29*1.01325d6)*0.2
        enddo
      enddo
!$OMP end parallel do

      if (thrust.GT.max_thrust) then
        max_thrust = thrust
        open  (unit=7000,file='Thrust.dat',status='unknown')
        WRITE (7000,*) tot,thrust
      endif

      

!$OMP parallel do
      do mb=1,1
        k=NkMM(MB)/2
        j=NjmM(MB)/2
        dd=qt(mb,1,j,k)
        uvel=qt(mb,2,j,k)/qt(mb,1,j,k)
        vvel=qt(mb,3,j,k)/qt(mb,1,j,k)
        Pre=GM1*(Qt(mb,4,J,K)-.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)-q0*Qt(mb,5,J,K))
        aa=dsqrt(dabs(pre*gm/dd))
        uM=dsqrt(uvel**2+vvel**2)/aa
      enddo
!$OMP end parallel do

      !write(7001,*)tot,uM

!$OMP parallel do
      do mb=1,1
        k=NkMM(MB)/2
        j=NjmM(MB)
        dd=qt(mb,1,j,k)
        uvel=qt(mb,2,j,k)/qt(mb,1,j,k)
        vvel=qt(mb,3,j,k)/qt(mb,1,j,k)
        Pre=GM1*(Qt(mb,4,J,K)-.5*(Qt(mb,2,J,K)**2+Qt(mb,3,J,K)**2)/Qt(mb,1,J,K)-q0*Qt(mb,5,J,K))
        aa=dsqrt(dabs(pre*gm/dd))
        uM=dsqrt(uvel**2+vvel**2)/aa
      enddo
!$OMP end parallel do

      !write(7002,*)tot,uM
      !write(7003,*)tot,pre
      
    enddo
    close (7000)

!-----------------------------------------------------------------------
!     PRINT/PLOT FINAL DATA
!-----------------------------------------------------------------------
!
!------------------MODIFIED NUMERICAL GRID TO REAL GRID-----------------

      !WRITE(*,*)'ALL BLOCKS CONVERGED!!'
      
      RETURN
      END 

!=======================================================================
!     subroutine fvsw
!=======================================================================
      
      SUBROUTINE fvsw(I,xiw,eta,QA,F)

      IMPLICIT NONE

      COMMON/GMD/GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1
      DIMENSION  QA(6,2),F(6,-1:1000),fLR(2,6)

      REAL*8 :: Q0,A1,A2,A12,D1,D2,QA,U1,U2,V1,V2,P1,P2,Z1,Z2,C1,C2,F
      REAL*8 :: DABS,PD1,PD2,UU1,UU2,TXY,XIW,ETA,CU1,CU2
      REAL*8 :: ALPHAL,ALPHAR,CML,CMR,FLR
      REAL*8 :: AUSMDU,AUSMDV,AUSMVU,AUSMVV,DP,COEF,ST
      REAL*8 :: GM,RGM,GM1,RGM1,GM12,RG2M1,RGGM1

      INTEGER :: N,I,K

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

!------------------2----------------------------------------------------

      d2=qa(1,2)
      u2=qa(2,2)
      v2=qa(3,2)
      p2=qa(4,2)
      z2=qa(5,2)

      c2=dsqrt(dabs(p2*1.4/qa(1,2)))
      pd2=p2/qa(1,2)
      uu2=(xiw*u2+eta*v2)

!------------------interface--------------------------------------------

      a12=(c1+c2)/2.
      cu1=uu1/txy/a12
      cu2=uu2/txy/a12
      alphaL=2*pd1/(pd1+pd2)
      alphaR=2*pd2/(pd1+pd2)

!------------------for particle based on density------------------------

      IF (ABS(cu1).Le.1.0) THEN
        CML=alphaL*((cu1+1.)**2/4-(cu1+abs(cu1))/2.)+(cu1+abs(cu1))/2.
        a1=a1+CML*d1
        a2=a2+p1*(cu1+1.)**2*(2-cu1)/4.
      ELSE
        CML=(cu1+abs(cu1))/2.
        a1=a1+CML*d1
        a2=a2+p1*(cu1+abs(cu1))/(2*cu1)
      endif

!------------------2----------------------------------------------------

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
!     subroutine HMSTH
!=======================================================================

      SUBROUTINE hmsth(u1,ur,ul,imax)

      IMPLICIT NONE

      DIMENSION  vpr(-1:2401,6) ,vpl(-1:2401,6)
      DIMENSION  vr(-1:2401,6) ,vl(-1:2401,6)
      DIMENSION  ur(-1:2401,6) ,ul(-1:2401,6)
      DIMENSION  u1(-1:2401,6)

      REAL*8 :: CK,HI,EPS,DTP,DTM,U1,SS,VPL,VPR,BVOF,EPSVOF,CFL,CBVOF
      REAL*8 :: DCOSH,TBVOF,DTANH,DIP0,DDP,UR,UL
      REAL*8 :: AVOFM,AVOF,SVOF,BB,CC,TBVOFP,TVOFP,TVOFM,VL,VR,ALPHA

      INTEGER :: IMM1,IMAX,I,K,J,JP1,JM1,JP,JM

      CK=1./3.
      imm1 = imax-1

      hi=1.

      eps=1e-14

!$OMP parallel do
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
!$OMP end parallel do

      bvof = 1.3
      epsvof = 1e-12
      CFL=1e-6

      cbvof = dcosh(bvof)
      tbvof = dtanh(bvof)

      imm1 =imax-1

!-----------------------------------------------------------------------

!$OMP parallel do
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

          BB = exp(svof*bvof*(2.0*(u1(j,k)-avofm+epsvof)/(avof+epsvof)-1.0))

          CC = (BB/cbvof - 1.0)/tbvof

          tbvofp = (CC+tbvof)/(1.0+CC*tbvof)
          tvofp = dlog( cosh(bvof*cfl) - dsinh(bvof*cfl)*tbvofp)/(bvof*(CFL+1e-12))
          tvofm = dlog( cosh(bvof*cfl) + dsinh(bvof*cfl)*CC)/(bvof*(CFL+1e-12))

          vl(j,k) = avofm + 0.5*avof*(1. - svof*tvofp)
          vr(j-1,k) = avofm + 0.5*avof*(1. + svof*tvofm)

          if(dsign(1.d0,(u1(jp,k)-u1(j,k))*(u1(j,k)-u1(jm,k))).eq.-1.d0) then
            vr(j-1,k)  = u1(j,k)
            vl(j,k) =  u1(j,k)
          endif
        ENDDO
      ENDDO
!$OMP end parallel do

!-----------------------------------------------------------------------

      eps =1e-12
      alpha =250.

!$OMP parallel do
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
!$OMP end parallel do

      RETURN
      END
