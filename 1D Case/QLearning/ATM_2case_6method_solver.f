!=======================================================================
! Initialization (Conservative variables)
!=======================================================================

      SUBROUTINE ini(icase,mesh_independent,u1,nt,l_x)
  
      integer,intent(in) :: icase
      integer,intent(in) :: mesh_independent

      integer,intent(out) :: nt
      real,intent(out) :: u1(0:13000,4)
      real,intent(out) :: l_x

      real :: pz(0:13000),uz(0:13000),q1z(0:13000),z(0:13000)
      
      integer :: imax,imp1
      real :: q0,r,dx,del,dt
      
c-----------------------------------------------------------------------

      if (icase .EQ. 1) then
        ! TIBW
        nt = 760 * mesh_independent
        q0 = 0.0
        l_x = 1.0
        dx = l_x/(100*mesh_independent)
        del = 0.005
        dt = del*dx

      elseif (icase .EQ. 2) then
        ! DWC
        nt = 200 * mesh_independent
        q0 = 0.5196d10
        l_x = 0.05
        dx = l_x/(100*mesh_independent)
        dt = 5d-10
        del = dt / dx

      endif

      imax = 100*mesh_independent+1
      imp1 = imax + 1
      r = 1.4

c-----------------------------------------------------------------------

      if (icase.eq.1) then
        do i = 0,(10*mesh_independent)
          uz(i) = 0.0
          q1z(i) = 1.0
          pz(i) = 1000.0
          z(i) = 0.0
        enddo
            
        do i = (10*mesh_independent+1),(imax-(10*mesh_independent+1))
          uz(i) = 0.0
          q1z(i) = 1.0
          pz(i) = 0.01
          z(i) = 0.0
        enddo
            
        do i = (imax-(10*mesh_independent)), imp1
          uz(i) = 0.0
          q1z(i) = 1.0 
          pz(i) = 100.0
          z(i) = 0.0
        enddo
      endif

c-----------------------------------------------------------------------  

      if (icase.eq.2) then
        do i = 0, (10*mesh_independent+1)
          uz(i) = 4.162d4 
          q1z(i) = 1.945d-3
          pz(i) = 8.27d6
          z(i) = 0.0
        enddo
       
        do i = (10*mesh_independent+2), imp1
          uz(i) = 0.0
          q1z(i) = 1.201d-3 
          pz(i) = 8.321d5
          z(i) = 1.0
        enddo
      endif
 
c-----------------Primitive -> Conservative-----------------------------  

      do i = 0, imp1
        u1(i, 1) = q1z(i)
        u1(i, 2) = q1z(i)*uz(i)
        u1(i, 3) = pz(i)/(r-1)+0.5*q1z(i)*uz(i)**2+q0*q1z(i)*z(i)
        u1(i, 4) = q1z(i)*z(i)
      enddo
      
      return
      end

!=======================================================================
! solver
!=======================================================================
      
      subroutine solver(icase,
     &               imethod,
     &               mesh_independent,
     &               hmsth_alpha,
     &               hmsth_bvof,
     &               un,
     &               u3)

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0, ddd
     
      real :: fp(0:13000,4),fn(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)
      real :: u0(0:13000,4)
      real :: pz(0:13000),uz(0:13000),q1z(0:13000)
      real :: z(0:13000)

      real :: p(0:13000)
      real :: temperature(0:13000)
      
      real :: source(0:13000,4)
      real :: source2(0:13000,4)
      real :: Qold1ds(0:13000,4)
      real :: Qold1s(0:13000,4)
      real :: Qold1d(0:13000,4)
      real :: Qold1(0:13000,4)
      real :: Qold2d(0:13000,4)
      real :: Qold2ss(0:13000,4)
      real :: Qolddss(0:13000,4)
      
      integer :: nt

      integer, intent(in) :: icase, imethod, mesh_independent
      real, intent(in) :: hmsth_alpha,hmsth_bvof
      real, intent(in) :: un(0:13000,4)
      real, intent(out) :: u3(0:13000,4)
      
      !icase = 1 TIBW
      !icase = 2 DWC
 
      !imethod =  1, MUSCL_PRIMITIVE
      !imethod =  2, MUSCL_CONSERVATIVE
      !imethod =  3, THINCEM_PRIMITIVE
      !imethod =  4, THINCEM_CONSERVATIVE
      !imethod =  5, THINC_PRIMITIVE
      !imethod =  6, THINC_CONSERVATIVE
      !imethod =  7, HMT_PRIMITIVE
      !imethod =  8, HMT_CONSERVATIVE
      !imethod =  9, HMT2_PRIMITIVE
      !imethod = 10, HMT2_CONSERVATIVE
      !imethod = 11, HMT(2021)_PRIMITIVE
      !imethod = 12, HMT(2021)_CONSERVATIVE

c---------------------STEP----------------------------------------------
      
      if (icase .EQ. 1) then
        ! TIBW
        nt = 760 * mesh_independent
        q0 = 0.0
        dx = 1.0/(100*mesh_independent)
        del = 0.005
        dt = del*dx

      elseif (icase .EQ. 2) then
        ! DWC
        nt = 200 * mesh_independent
        q0 = 0.5196d10
        dx = 0.05/(100*mesh_independent)
        dt = 5d-10
        dt = dt / mesh_independent
        del = dt / dx

      endif

c-----------------------------------------------------------------------

      interp =1
      ifrog =10

c---------------------constant------------------------------------------

      half = 5d-1
      one = 10d-1
      two = one + one
      zero = 0.0
      ten = 1d1
      kflux = 2
      pi=4*atan(1.)

c---------------------constant------------------------------------------

      temperature0 = 0.1155d10
      ddd = 1.0/0.5825d10
      r = 1.4
      r1 = 0.4

!=======================================================================
! Flow field setup (Initialize Primitive Variables)
!=======================================================================

      imax = 100*mesh_independent+1
      imp1 = imax + 1
      imm1 = imax - 1
      imm2 = imax - 2
      imm3 = imax - 3
      imm4 = imax - 4
      imm5 = imax - 5    

c---------------------source ONE---------------------S(k/2)-------------

      call SOURCE_TERM(un, Qold1ds, Qold1s, 
     &                   source, source2,
     &                   p, temperature)

      if (icase .EQ. 1) then
        call BC_TIBW(Qold1s, icase)
      elseif (icase .EQ. 2) then
        call BC_DWC(Qold1s, icase)
      endif

c---------------------MAIN ONE-----------------------SF(k)--------------

      if (imethod .EQ. 1) then
        call MUSCL_P(Qold1s,upr,upl,unr,unl)

      elseif (imethod .EQ. 2) then
        call MUSCL_C(Qold1s,upr,upl,unr,unl)

      elseif (imethod .EQ. 3) then
        call THINCEM_P(Qold1s,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 4) then
        call THINCEM_C(Qold1s,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 5) then
        call THINC_P(Qold1s,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 6) then
        call THINC_C(Qold1s,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 7) then
        call HMSTH_P(Qold1s,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 8) then
        call HMSTH_C(Qold1s,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 9) then
        call HMSTH2_P(Qold1s,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 10) then
        call HMSTH2_C(Qold1s,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 11) then
        call HMT_P(Qold1s,upr,upl,unr,unl,hmsth_bvof,imax)

      elseif (imethod .EQ. 12) then
        call HMT_C(Qold1s,upr,upl,unr,unl,hmsth_bvof,imax)

      endif

      if (icase .EQ. 1) then
        call ausm(fp,upr,upl)
        call ausm(fn,unr,unl)
      elseif (icase .EQ. 2) then
        call ausm2(fp,upr,upl)
        call ausm2(fn,unr,unl)
      endif

      do k = 1, 4
        do i = 3, imm2

          Qold1d(i, k) = -(dt/dx)*(fp(i,k)-fn(i,k)) 
          Qold1(i, k) = Qold1s(i,k)+Qold1d(i,k)

        enddo
      enddo

      if (icase .EQ. 1) then
        call BC_TIBW(Qold1, icase)
      elseif (icase .EQ. 2) then
        call BC_DWC(Qold1, icase)
      endif

c---------------------MAIN TWO-----------------------SF(k)--------------

      if (imethod .EQ. 1) then
        call MUSCL_P(Qold1,upr,upl,unr,unl)

      elseif (imethod .EQ. 2) then
        call MUSCL_C(Qold1,upr,upl,unr,unl)

      elseif (imethod .EQ. 3) then
        call THINCEM_P(Qold1,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 4) then
        call THINCEM_C(Qold1,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 5) then
        call THINC_P(Qold1,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 6) then
        call THINC_C(Qold1,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 7) then
        call HMSTH_P(Qold1,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 8) then
        call HMSTH_C(Qold1,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 9) then
        call HMSTH2_P(Qold1,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 10) then
        call HMSTH2_C(Qold1,upr,upl,unr,unl,hmsth_alpha,hmsth_bvof)

      elseif (imethod .EQ. 11) then
        call HMT_P(Qold1,upr,upl,unr,unl,hmsth_bvof)

      elseif (imethod .EQ. 12) then
        call HMT_C(Qold1,upr,upl,unr,unl,hmsth_bvof)

      endif
        
      if (icase .EQ. 1) then
        call ausm(fp,upr,upl)
        call ausm(fn,unr,unl)
      elseif (icase .EQ. 2) then
        call ausm2(fp,upr,upl)
        call ausm2(fn,unr,unl)
      endif

      do k=1,4
        do i=3,imm2

          Qold2d(i,k)=-(dt/dx)*(fp(i,k)-fn(i,k))
          Qold2ss(i,k)=Qold1s(i,k)+half*(Qold1d(i,k)+Qold2d(i,k))

        enddo
      enddo

      if (icase .EQ. 1) then
        call BC_TIBW(Qold2ss, icase)
      elseif (icase .EQ. 2) then
        call BC_DWC(Qold2ss, icase)
      endif

c---------------------source TWO---------------------SF(k/2)------------

      call SOURCE_TERM(Qold2ss, Qolddss, u3, 
     &                   source, source2,
     &                   p, temperature)

      if (icase .EQ. 1) then
        call BC_TIBW(u3, icase)
      elseif (icase .EQ. 2) then
        call BC_DWC(u3, icase)
      endif

c-----------------------------------------------------------------------
 
c      call output(icase,
c     &            imethod,
c     &            mesh_independent,
c     &            hmsth_alpha,
c     &            hmsth_bvof,
c     &            file_name,
c     &            u3,
c     &            nt)

c-----------------------------------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine output
!=======================================================================

      subroutine output(icase,
     &                  imethod,
     &                  mesh_independent,
     &                  hmsth_alpha,
     &                  hmsth_bvof,
     &                  file_name,
     &                  u3,
     &                  nt)

c---------------------output--------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0, ddd
     
      integer :: nt
      integer :: icase, imethod, mesh_independent
      real :: hmsth_alpha,hmsth_bvof
      real :: u3(0:13000,4)

c-----------------------------------------------------------------------

      character(len = 150) file_name
      character(len = 20) file_case
      character(len = 20) file_method
      character(len = 30) file_mesh
      character(len = 40) file_hmsth     
      character(len = 15) mesh_independent_str
      character(len = 15) hmsth_alpha_str
      character(len = 15) hmsth_bvof_str

c---------------------Cases---------------------------------------------

      if (icase.EQ.1) then
        file_case = 'CASE1'
      elseif (icase.EQ.2) then
        file_case = 'CASE2'
      elseif (icase.EQ.3) then
        file_case = 'CASE3'
      elseif (icase.EQ.4) then
        file_case = 'CASE4'
      else
        file_case = 'CASEErr'
      endif
      
      file_case = adjustl(file_case)

c---------------------Reconstruction method-----------------------------

      if (imethod.EQ.1) then
        file_method = "MUSCL"
      elseif (imethod.EQ.2) then
        file_method = "THINCEM"
      elseif (imethod.EQ.3) then
        file_method = "HYBRID"
      else
        file_method = "methodErr"
      endif
      file_method = adjustl(file_method)

c---------------------Build---------------------------------------------

      write(mesh_independent_str, "(I10)") mesh_independent
      mesh_independent_str = adjustl(mesh_independent_str)
      file_mesh = "MESH"//mesh_independent_str

      write(hmsth_alpha_str, "(F12.4)") hmsth_alpha ! float2str
      write(hmsth_bvof_str,  "(F12.4)") hmsth_bvof  ! float2str
      hmsth_alpha_str = adjustl(hmsth_alpha_str) ! move2left
      hmsth_bvof_str = adjustl(hmsth_bvof_str)   ! move2left
       
      if (imethod .EQ. 3) then
        file_hmsth = "a"//trim(hmsth_alpha_str)//"_"//
     &               "b"//trim(hmsth_bvof_str)//"-"
      else
        file_hmsth = ""
      endif

      file_hmsth = adjustl(file_hmsth)

      file_name = trim(file_method)//"-"//
     &            trim(file_hmsth)//
     &            trim(file_case)//"-"//
     &            trim(file_mesh)//".csv"
      
      file_name = trim(adjustl(file_name))
      
      open(21,file=file_name)
      
      do i=1,imax
        xx=(i-1)*dx
        y1=u3(i,1) !rho
        y2=u3(i,2) !rho*u
        y3=u3(i,3)
        dd=u3(i,1) !rho
        ul=u3(i,2)/dd !u   
        pp=0.4*(u3(i,3)-0.5*u3(i,2)*ul-q0*u3(i,4))
        T=pp/dd
        zz=u3(i,4)/u3(i,1)
        write(21,*) xx,y1,pp,T,ul,zz !xx,rho,p,T,u,z
      enddo

      close(21)

      time=dt*nt

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine SOURCE_TERM
!=======================================================================

      SUBROUTINE SOURCE_TERM(Q, QDS, QS, 
     &                       source, source2,
     &                       p, temperature)

c---------------------SOURCE TERM---------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0, ddd
      
      real :: Q(0:13000, 4)
      real :: QDS(0:13000, 4)
      real :: QS(0:13000, 4)
      real :: source(0:13000, 4)
      real :: source2(0:13000, 4)
      
      real :: p(0:13000)
      real :: temperature(0:13000)
      real :: hh(0:13000)

c-----------------------------------------------------------------------

      half = 5d-1
      one = 10d-1
      two = one + one
      zero = 0.0
      ten = 1d1

      do i = 1, imax !change to all Qold WITHOUT P(i)
        p(i) = (r-one)*
     &  (     
     &    Q(i,3)  !p/(r-1) * 0.5*rho*u^2 + q0 rho z 
     &    - half*(Q(i,2)*Q(i,2)) / Q(i,1) ! 0.5(rho^2 u^2)/rho 
     &    - q0*Q(i,4) ! q0*rho*z
     &  )
      enddo

      do i = 1, imax
        temperature(i) = p(i)/Q(i,1)
      enddo
        
      do i = 1, imax
        if ((temperature(i)-temperature0)>zero) then
          hh(i) = one
        else
          hh(i) = zero  
        endif
      enddo

      do i = 1, imax
        source(i, 1) = 0.0
        source(i, 2) = 0.0
        source(i, 3) = 0.0
        source(i, 4) = -(one/ddd)*Q(i,4)*hh(i)
      enddo

      do i = 1, imax
        source2(i, 1) = 0.0
        source2(i, 2) = 0.0
        source2(i, 3) = 0.0
        source2(i, 4) = -one*hh(i)
      enddo

      do k = 1, 4
        do i = 1, imax
          QDS(i,k) = half*dt*source(i,k)/
     #    (one-half*half*dt*source2(i,k) )
          QS(i,k) = Q(i,k)+QDS(i,k)
        enddo
      enddo

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine BC TIBW
!=======================================================================

      SUBROUTINE BC_TIBW(Q, icase)

c---------------------BC_TIBW-------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0, temperature0, ddd
      
      real :: Q(0:13000,4)
      integer :: icase

c-----------------------------------------------------------------------

      if (icase /= 5) then
        Q(2, 2) = -Q(4, 2)   !not tube
      elseif (icase == 5) then 
        Q(2, 2) = -Q(3, 2)  !tube
      endif
        
      Q(1, 2) = -Q(5, 2)
      Q(0, 2) = -Q(6, 2)
      Q(imm1, 2) = -Q(imm3, 2)
      Q(imax, 2) = -Q(imm4, 2)
      Q(imp1, 2) = -Q(imm5, 2)

      Q(2, 1) = Q(4, 1)
      Q(1, 1) = Q(5, 1)
      Q(0, 1) = Q(6, 1)
      Q(imm1, 1) = Q(imm3, 1)
      Q(imax, 1) = Q(imm4, 1)
      Q(imp1, 1) = Q(imm5, 1)

      Q(2, 3) = Q(4, 3)
      Q(1, 3) = Q(5, 3)
      Q(0, 3) = Q(6, 3)
      Q(imm1, 3) = Q(imm3, 3)
      Q(imax, 3) = Q(imm4, 3)
      Q(imp1, 3) = Q(imm5, 3)
          
      Q(2, 4) = Q(3, 4)
      Q(1, 4) = Q(2, 4)
      Q(0, 4) = Q(1, 4)
      Q(imm1, 4) = Q(imm3, 4)
      Q(imax, 4) = Q(imm4, 4)
      Q(imp1, 4) = Q(ima5, 4)
      
c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine BC DWC
!=======================================================================

      SUBROUTINE BC_DWC(Q, icase)

c---------------------BC_DWC--------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0, temperature0, ddd
      
      dimension Q(0:13000,4)
      integer::icase

c-----------------------------------------------------------------------

      if (icase /= 5) then
        Q(2, 2) = Q(3, 2)   !not tube
      elseif (icase == 5) then 
        Q(2, 2) = -Q(3, 2)  !tube
      endif

      Q(1, 2) = Q(2, 2)
      Q(0, 2) = Q(1, 2)
      Q(imm1, 2) = Q(imm2, 2)
      Q(imax, 2) = Q(imm1, 2)
      Q(imp1, 2) = Q(imax, 2)

      Q(2, 1) = Q(3, 1)
      Q(1, 1) = Q(2, 1)
      Q(0, 1) = Q(1, 1)
      Q(imm1, 1) = Q(imm2, 1)
      Q(imax, 1) = Q(imm1, 1)
      Q(imp1, 1) = Q(imax, 1)

      Q(2, 3) = Q(3, 3)
      Q(1, 3) = Q(2, 3)
      Q(0, 3) = Q(1, 3)
      Q(imm1, 3) = Q(imm2, 3)
      Q(imax, 3) = Q(imm1, 3)
      Q(imp1, 3) = Q(imax, 3)
          
      Q(2, 4) = Q(3, 4)
      Q(1, 4) = Q(2, 4)
      Q(0, 4) = Q(1, 4)
      Q(imm1, 4) = Q(imm2, 4)
      Q(imax, 4) = Q(imm1, 4)
      Q(imp1, 4) = Q(imax, 4)
      
c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine AUSM+UP
!=======================================================================

      SUBROUTINE ausm(fp,qr,ql)

c---------------------AUSM+UP-------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0, temperature0, ddd

      real qr(0:13000,4),ql(0:13000,4)
      real qa(2,0:13000,4)
      real fp(0:13000,4)

c-----------------------------------------------------------------------

      do k=1,4
        do  I=3,imm2
          qa(2,i,k)=qr(i,k)
          qa(1,i,k)=ql(i,k)
        enddo
      enddo

      DO  I=3,imm2
        d1=qa(1,i,1)
        u1=qa(1,i,2)/qa(1,i,1)
        pL=0.4*(qa(1,i,3)-.5*qa(1,i,1)*u1*u1)
        c1=sqrt(abs(pL*1.4/qa(1,i,1)))
        hL=pL*3.5/d1+.5*(u1*u1)

        d2= qa(2,i,1)
        u2=qa(2,i,2)/qa(2,i,1)
        pR=0.4*(qa(2,i,3)-.5*qa(2,i,1)*u2*u2)
        c2=sqrt(abs(pR*1.4/qa(2,i,1)))
        hr=pr*3.5/d2+.5*(u2*u2)

c---------------------2-------------------------------------------------

        a12=(c1+c2)/2.
        xmlp = u1/a12
        xmrp = u2/a12

c---------------------Switching Functions-------------------------------

        all = 0.5*(1. + sign(1.0,xmlp))  
        alr = 0.5*(1. - sign(1.0,xmrp))  

        btl = -amax1(0.0,1.0-int(abs(xmlp)))
        btr = -amax1(0.0,1.0-int(abs(xmrp)))

c---------------------Split and interface Mach numbers------------------

        xmml = 0.25*(xmlp+1.0)**2 
        xmmr = -0.25*(xmrp-1.0)**2 
        xmbsl = 0.5*(xmlp + abs(xmlp))
        xmbsr = 0.5*(xmrp - abs(xmrp))
        xmaddp = btl*(xmbsl - xmml)
        xmaddm = btr*(xmbsr - xmmr)
        xmhalf = xmbsl + xmaddp + xmbsr + xmaddm

        xmplus = all*(1.+btl)*xmlp - btl*xmml
        xmmins = alr*(1.+btr)*xmrp - btr*xmmr 

        xm12=.5*(xmlp+xmrp)
        
        dxm=1-xm12**2
        dkp=0.1
        fprec = dkp*amax1(dxm,0.)*(xmaddp - xmaddm)*(PL - PR)/a12
        xmface = a12*(d1*xmplus + d2*xmmins + fprec)

c---------------------Interface mass flux-------------------------------

        fml =  0.5*(xmface + abs(xmface))
        fmr =  0.5*(xmface - abs(xmface))

c---------------------Pressure splitting--------------------------------

        ppl = 0.25*(xmlp+1.)**2*(2.-xmlp)
        ppr = 0.25*(xmrp-1.)**2*(2.+xmrp)

        pnet = (all*(1.+btl)-btl*ppl)*pL+(alr*(1.+btr)-btr*ppr)*pr

c-----------------------------------------------------------------------
   
        fp(i,1) = (fml + fmr)
        fp(i,2) = (fml*u1 + fmr*u2 + pnet)
        fp(i,3) = (fml*hL + fmr*hr)
        fp(i,4) = 0
      
      enddo

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine ausm2
!=======================================================================

      SUBROUTINE ausm2(fp,qr,ql)

c---------------------AUSM+---------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0, temperature0, ddd

      real qr(0:13000,4),ql(0:13000,4)
      real a1(0:13000),a2(0:13000),a3(0:13000)
      real ac(0:13000)
      real qa(2,0:13000,4)
      real fp(0:13000,4)
      real fLR(2,0:13000,4)

      alpha=3/16.0
      beta=1.0/8.0

      do k=1,4
        do  I=3,imm2
          qa(2,i,k)=qr(i,k)
          qa(1,i,k)=ql(i,k)
        enddo
      enddo

      do i=3,imm2
        a1(i)=0.0
        a2(i)=0.0
        a3(i)=1.0
      enddo

      do i=3,imm2
        u1=qa(1,i,2)/qa(1,i,1)
        p1=0.4*(qa(1,i,3)-0.5*qa(1,i,1)*u1*u1-q0*qa(1,i,4))
        c1=sqrt(abs(p1*1.4/qa(1,i,1)))

        u2=qa(2,i,2)/qa(2,i,1)
        p2=0.4*(qa(2,i,3)-0.5*qa(2,i,1)*u2*u2-q0*qa(2,i,4))
     
        c2=sqrt(abs(p2*1.4/qa(2,i,1)))

c---------------------2-------------------------------------------------

        a12=(c1+c2)/2.0
        ac(i)=a12
        cu1=u1/a12
        cu2=u2/a12

        if (ABS(cu1).LT.1.0) then
          a1(i)=a1(i)+a3(i)*0.25*(cu1+a3(i))**2
     #    +beta*a3(i)*(cu1**2-1)**2
          a2(i)=a2(i)+0.25*p1*((cu1+a3(i))**2
     $    *(2.0-a3(i)*cu1)
     $    +alpha*a3(i)*cu1*(cu1*cu1-1)**2)
        else
          a1(i)=a1(i)+0.5*(cu1+a3(i)*abs(cu1))
          a2(i)=a2(i)+0.5*p1*(cu1+a3(i)
     $    *abs(cu1))/cu1
        endif

c---------------------2-------------------------------------------------

        a3(i)=-1.0*a3(i)

        if (ABS(cu2).LT.1.0) then
          a1(i)=a1(i)+a3(i)*0.25*(cu2+a3(i))**2
     #    +beta*a3(i)*(cu2**2-1.0)**2
          a2(i)=a2(i)+0.25*p2*((cu2+a3(i))**2
     $    *(2.0-a3(i)*cu2)
     $    +alpha*a3(i)*cu2*(cu2*cu2-1.0)**2)
        else
          a1(i)=a1(i)+0.5*(cu2+a3(i)*abs(cu2))
          a2(i)=a2(i)+0.5*p2*(cu2+a3(i)
     $    *abs(cu2))/cu2
        endif

c-----------------------------------------------------------------------

       fLR(1,i,1)=qa(1,i,1)
       fLR(1,i,2)=qa(1,i,2)
       fLR(1,i,3)=p1*2.5+.5*qa(1,i,1)*u1*u2+p1+q0*qa(1,i,4)
       fLR(1,i,4)=qa(1,i,4) !--

       fLR(2,i,1)=qa(2,i,1)
       fLR(2,i,2)=qa(2,i,2)
       fLR(2,i,3)=p2*2.5+.5*qa(2,i,1)*u2*u2+p2+q0*qa(2,i,4)
       fLR(2,i,4)=qa(2,i,4) !--
      enddo

c-----------------------------------------------------------------------

      do  k=1,4
        do  I=3,imm2
          fp(i,k)=.5*a1(i)*ac(i)
     $    *(fLR(1,i,k)+fLR(2,i,k))
     $    -abs(.5*a1(i))*ac(i)
     $    *(fLR(2,i,k)-fLR(1,i,k))
        enddo
      enddo

      do i = 3, imm2
        fp(i, 2) = fp(i, 2) + a2(i)
      enddo

c---------------------end-----------------------------------------------

      return
      end

!=======================================================================
! subroutine MUSCL_PRIMITIVE
!=======================================================================

      SUBROUTINE MUSCL_P(u0,upr,upl,unr,unl)

c---------------------MUSCL---------------------------------------------

      real :: u1(0:13000,4)
      real :: u0(0:13000,4)
      real :: vpr(0:13000,4),vpl(0:13000,4)
      real :: vnr(0:13000,4),vnl(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)   
      real :: niu

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      call C2P(u0,u1)

c-----------------------------------------------------------------------
      
      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      call P2C(vpr,vpl,vnr,vnl,upr,upl,unr,unl)

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine MUSCL_CONSERVATIVE
!=======================================================================

      SUBROUTINE MUSCL_C(u1,upr,upl,unr,unl)

c---------------------MUSCL---------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0, ddd

      real :: u1(0:13000,4)
      real :: vpr(0:13000,4),vpl(0:13000,4)
      real :: vnr(0:13000,4),vnl(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)   
      real :: niu

c-----------------------------------------------------------------------

      niu=1./3.

      do k=1,4
        do i=3,imm2
          r1=u1(i+2,k)-u1(i+1,k)
          r2=u1(i+1,k)-u1(i,k)
          r3=u1(i,k)-u1(i-1,k)
          r4=u1(i-1,k)-u1(i-2,k)

          upr(i,k)=u1(i+1,k)
     #    -1./4.*((1-niu)*smd(r1,r2)+(1+niu)*smd(r2,r1))
          upl(i,k)=u1(i,k)
     #    +1./4.*((1-niu)*smd(r3,r2)+(1+niu)*smd(r2,r3))
          unr(i,k)=u1(i,k)
     #    -1./4.*((1-niu)*smd(r2,r3)+(1+niu)*smd(r3,r2))
          unl(i,k)=u1(i-1,k)
     #    +1./4.*((1-niu)*smd(r4,r3)+(1+niu)*smd(r3,r4))
        enddo
      enddo

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! Limited Function
!=======================================================================

      function smd(d1,d2)
      ddd=1.e-6
      smd=(d2*(d1**2+ddd)+d1*(d2**2+ddd))/(d1**2+d2**2+2*ddd)
      return
      end

!=======================================================================
! subroutine THINCEM_PRIMITIVE
!=======================================================================

      SUBROUTINE THINCEM_P(u0,upr,upl,unr,unl,bvof)

c---------------------THINCEM-------------------------------------------
   
      real :: u1(0:13000,4)
      real :: u0(0:13000,4)
      real :: vpr(0:13000,4),vpl(0:13000,4)
      real :: vnr(0:13000,4),vnl(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      call C2P(u0,u1)

c-----------------------------------------------------------------------

      call THINCEM_C(u1,vpr,vpl,vnr,vnl,bvof)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      call P2C(vpr,vpl,vnr,vnl,upr,upl,unr,unl)

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine THINCEM_CONSERVATIVE
!=======================================================================

      SUBROUTINE THINCEM_C(u1,upr,upl,unr,unl,bvof)

c---------------------THINCEM-------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0,ddd
     
      real :: u1(0:13000,4)
      real :: vpr(0:13000,4),vpl(0:13000,4)
      real :: vnr(0:13000,4),vnl(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)
 
c-----------------------------------------------------------------------
      
      epsvof = 10e-20
      CFL = 1e-6
      cbvof = cosh(bvof)
      tbvof = tanh(bvof)

c-----------------------------------------------------------------------

      do k=1,4
        do j=2,imm1

          jp=j+1
          jm=j-1

          jp = min(jp,imax)
          jm = max(jm,1)

          avofm = min(u1(jp,k),u1(jm,k))
          avof =  max(u1(jp,k),u1(jm,k))-avofm
          svof = -1.0

          if(avofm.eq.u1(jm,k)) svof = 1.0
            BB = exp(svof*bvof*(2.0*(u1(j,k)-avofm+epsvof)
     #         /(avof+epsvof)-1.0))

            CC = (BB/cbvof - 1.0)/tbvof

            tbvofp = (CC+tbvof)/(1.0+CC*tbvof)
            tvofp = log( cosh(bvof*cfl) - sinh(bvof*cfl)*tbvofp)
     #            /(bvof*(CFL+1e-12))
            tvofm = log( cosh(bvof*cfl) + sinh(bvof*cfl)*CC)
     #            /(bvof*(CFL+1e-12))

            upl(j,k) = avofm + 0.5*avof*(1. - svof*tvofp)
            upr(j-1,k) = avofm + 0.5*avof*(1. + svof*tvofm)

            if(sign(1.0,(u1(jp,k)-u1(j,k))
     #      *(u1(j,k)-u1(jm,k))).eq.-1.0) then

            upr(j-1,k) = u1(j,k)
            upl(j,k) = u1(j,k)

          endif
        enddo
      enddo

c-----------------------------------------------------------------------

      do k=1,4
        do j=2,imm1

          unr(j,k) = upr(j-1,k)
          unl(j,k) = upl(j-1,k)

        enddo
      enddo

c---------------------end-----------------------------------------------
      
      RETURN
      END


!=======================================================================
! subroutine THINC_PRIMITIVE
!=======================================================================

      SUBROUTINE THINC_P(u0,upr,upl,unr,unl,bvof)

c---------------------THINCEM-------------------------------------------
   
      real :: u1(0:13000,4)
      real :: u0(0:13000,4)
      real :: vpr(0:13000,4),vpl(0:13000,4)
      real :: vnr(0:13000,4),vnl(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      call C2P(u0,u1)

c-----------------------------------------------------------------------

      call THINC_C(u1,vpr,vpl,vnr,vnl,bvof)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      call P2C(vpr,vpl,vnr,vnl,upr,upl,unr,unl)

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine THINC_CONSERVATIVE
!=======================================================================

      SUBROUTINE THINC_C(u1,upr,upl,unr,unl,bvof)

c---------------------THINCEM-------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0,ddd
     
      real :: u1(0:13000,4)
      real :: upr(0:13000,4),upl(0:13000,4)
      real :: unr(0:13000,4),unl(0:13000,4)

      real :: thinclim
      real :: niu
      real :: zbar,zmax,zmin,zeps,zeps1
      real :: dchil,tsgn,beta_test
      real :: dup,ddw
      real :: BB,CC
      real :: qnl,qnr,qpl,qpr

c-----------------------------------------------------------------------
      
      zeps = 1.0e-20
      zeps1 = 1.0e-3
      beta_test = bvof
      thinclim = 1.0
      niu=.333333333

c-----------------------------------------------------------------------

      DO k=1,3
        DO i=3,imm2

c---------------------upr-----------------------------------------------

          zp = u1(i+2 ,k)
          z0 = u1(i+1 ,k)
          zm = u1(i   ,k)

          zmin = min(zm,zp)
          zmax = max(zm,zp)
            
          IF((zp-zm)>0.0)THEN
            tsgn = 1.0    
          ELSE
            tsgn = -1.0
          END IF
         
          zbar  = (z0-zmin+zeps)/(zmax-zmin+zeps)
          BB    = exp(tsgn*beta_test*(2.0*zbar-1.0))
          CC    = (BB/cosh(beta_test)-1.0)/tanh(beta_test)
          dchil = 0.5*log((1.0-CC)/(1.0+CC))/beta_test

          thinclim = dchil
 
          IF(((z0-zm)*(zp-z0) > 0.0))THEN
            qpr = zmin+0.5*(zmax-zmin)
     #          *(1.0+tsgn*tanh(-1.0*beta_test*dchil))
            thinclim = dchil
          ELSE
            qpr = z0
          END IF

c---------------------unr,upl-------------------------------------------

          zp = u1(i+1 ,k)
          z0 = u1(i   ,k)
          zm = u1(i-1 ,k)

          zmin = min(zm,zp)
          zmax = max(zm,zp)
            
          IF((zp-zm)>0.0)THEN
            tsgn = 1.0    
          ELSE
            tsgn = -1.0
          END IF
         
          zbar  = (z0-zmin+zeps)/(zmax-zmin+zeps)
          BB    = exp(tsgn*beta_test*(2.0*zbar-1.0))
          CC    = (BB/cosh(beta_test)-1.0)/tanh(beta_test)
          dchil = 0.5*log((1.0-CC)/(1.0+CC))/beta_test

          thinclim = dchil
 
          IF(((z0-zm)*(zp-z0) > 0.0))THEN
            qnr = zmin+0.5*(zmax-zmin)
     #          *(1.0+tsgn*tanh(-1.0*beta_test*dchil))
            qpl = zmin+0.5*(zmax-zmin)
     #          *(1.0+tsgn*tanh(beta_test*(1.0-dchil)))
            thinclim = 1
          ELSE
            qnr = z0
            qpl = z0
          END IF

c---------------------unl-----------------------------------------------

          zp = u1(i   ,k)
          z0 = u1(i-1 ,k)
          zm = u1(i-2 ,k)

          zmin = min(zm,zp)
          zmax = max(zm,zp)
            
          IF((zp-zm)>0.0)THEN
            tsgn = 1.0    
          ELSE
            tsgn = -1.0
          END IF
         
          zbar  = (z0-zmin+zeps)/(zmax-zmin+zeps)
          BB    = exp(tsgn*beta_test*(2.0*zbar-1.0))
          CC    = (BB/cosh(beta_test)-1.0)/tanh(beta_test)
          dchil = 0.5*log((1.0-CC)/(1.0+CC))/beta_test

          thinclim = dchil
 
          IF(((z0-zm)*(zp-z0) > 0.0))THEN
            qnl = zmin+0.5*(zmax-zmin)
     #          *(1.0+tsgn*tanh(beta_test*(1.0-dchil)))
            thinclim = 1
          ELSE
            qnl = z0
          END IF

c-----------------------------------------------------------------------

            upr(i,k) = qpr
            upl(i,k) = qpl
            unr(i,k) = qnr
            unl(i,k) = qnl

c-----------------------------------------------------------------------

        ENDDO
      ENDDO

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine HMSTH_PRIMITIVE
!=======================================================================

      SUBROUTINE HMSTH_P(u0,ur,ul,unr,unl,alpha,bvof)

c---------------------HMSTH---------------------------------------------

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  urh(0:13000,4) ,ulh(0:13000,4)
      real  unrh(0:13000,4) ,unlh(0:13000,4)

      real  u0(0:13000,4) ,u1(0:13000,4)

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      call C2P(u0,u1)

c---------------------MUSCL---------------------------------------------

      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------THINCEM-------------------------------------------

      call THINCEM_C(u1,vr,vl,vrn,vln,bvof)

c---------------------combine together---------------------------------- 

      call HMSTH_combine(u1,
     &                   vpr,vpl,vnr,vnl,
     &                   vr,vl,vrn,vln,
     &                   urh,ulh,unrh,unlh,
     &                   alpha)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      call P2C(urh,ulh,unrh,unlh,ur,ul,unr,unl)

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine HMSTH_CONSERVATIVE
!=======================================================================

      SUBROUTINE HMSTH_C(u1,ur,ul,unr,unl,alpha,bvof)

c---------------------HMSTH---------------------------------------------

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  u1(0:13000,4)

c---------------------MUSCL---------------------------------------------

      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------THINCEM-------------------------------------------

      call THINCEM_C(u1,vr,vl,vrn,vln,bvof)

c---------------------combine together---------------------------------- 

      call HMSTH_combine(u1,
     &                   vpr,vpl,vnr,vnl,
     &                   vr,vl,vrn,vln,
     &                   ur,ul,unr,unl,
     &                   alpha)

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine HMSTH2_PRIMITIVE
!=======================================================================

      SUBROUTINE HMSTH2_P(u0,ur,ul,unr,unl,alpha,bvof)

c---------------------HMSTH---------------------------------------------

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  urh(0:13000,4) ,ulh(0:13000,4)
      real  unrh(0:13000,4) ,unlh(0:13000,4)

      real  u0(0:13000,4) ,u1(0:13000,4)

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      call C2P(u0,u1)

c---------------------MUSCL---------------------------------------------

      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------THINCEM-------------------------------------------

      call THINCEM_C(u1,vr,vl,vrn,vln,bvof)

c---------------------combine together---------------------------------- 

      call HMSTH2_combine(u1,
     &                    vpr,vpl,vnr,vnl,
     &                    vr,vl,vrn,vln,
     &                    urh,ulh,unrh,unlh,
     &                    alpha)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      call P2C(urh,ulh,unrh,unlh,ur,ul,unr,unl)

c---------------------end-----------------------------------------------
      
      RETURN
      END

!=======================================================================
! subroutine HMSTH2_CONSERVATIVE
!=======================================================================

      SUBROUTINE HMSTH2_C(u1,ur,ul,unr,unl,alpha,bvof)

c---------------------HMSTH---------------------------------------------

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  u1(0:13000,4)

c---------------------MUSCL---------------------------------------------

      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------THINCEM-------------------------------------------

      call THINCEM_C(u1,vr,vl,vrn,vln,bvof)

c---------------------combine together---------------------------------- 

      call HMSTH2_combine(u1,
     &                    vpr,vpl,vnr,vnl,
     &                    vr,vl,vrn,vln,
     &                    ur,ul,unr,unl,
     &                    alpha)

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine HMT_PRIMITIVE
!=======================================================================

      SUBROUTINE HMT_P(u0,ur,ul,unr,unl,bvof)

c---------------------HMT-----------------------------------------------

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  urh(0:13000,4) ,ulh(0:13000,4)
      real  unrh(0:13000,4) ,unlh(0:13000,4)

      real  u0(0:13000,4) ,u1(0:13000,4)

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      call C2P(u0,u1)

c---------------------MUSCL---------------------------------------------

      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------THINCEM-------------------------------------------

      call THINCEM_C(u1,vr,vl,vrn,vln,bvof)

c---------------------combine together----------------------------------

      call HMT_combine(u1,
     &                 vpr,vpl,vnr,vnl,
     &                 vr,vl,vrn,vln,
     &                 urh,ulh,unrh,unlh)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      call P2C(urh,ulh,unrh,unlh,ur,ul,unr,unl)

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine HMT_CONSERVATIVE
!=======================================================================

      SUBROUTINE HMT_C(u1,ur,ul,unr,unl,bvof)

c---------------------HMT-----------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imp1,q0,temperature0, ddd

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  u1(0:13000,4)
      
c---------------------MUSCL---------------------------------------------

      call MUSCL_C(u1,vpr,vpl,vnr,vnl)

c---------------------THINCEM-------------------------------------------

      call THINCEM_C(u1,vr,vl,vrn,vln,bvof)
      
c---------------------combine together----------------------------------

      call HMT_combine(u1,
     &                 vpr,vpl,vnr,vnl,
     &                 vr,vl,vrn,vln,
     &                 ur,ul,unr,unl)       

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine CONSERVATIVE -> PRIMITIVE
!=======================================================================

      SUBROUTINE C2P(u0,u1)

c-----------------------------------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imax,imp1,q0,temperature0, ddd

      real :: u1(0:13000,4)
      real :: u0(0:13000,4)

c---------------------CONSERVATIVE -> PRIMITIVE-------------------------

      do i=1,imm2+2

        u1(i,1)=u0(i,1)
        u1(i,2)=u0(i,2)/u0(i,1)
        u1(i,4)=u0(i,4)/u0(i,1)        
        u1(i,3)=.4*(u0(i,3)-0.5*u0(i,2)*u1(i,2)-q0*u0(i,4))

      enddo

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine PRIMITIVE -> CONSERVATIVE
!=======================================================================

      SUBROUTINE P2C(vpr,vpl,vnr,vnl,upr,upl,unr,unl)

c-----------------------------------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imp1,q0,temperature0, ddd

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  upr(0:13000,4) ,upl(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

c---------------------PRIMITIVE -> CONSERVATIVE-------------------------

      do i=3,imm2
        upr(i,1)=vpr(i,1)
        upr(i,2)=vpr(i,1)*vpr(i,2)
        upr(i,4)=vpr(i,1)*vpr(i,4)
        upr(i,3)=vpr(i,3)*2.5+.5*vpr(i,2)*upr(i,2)+q0*upr(i,4)

        upl(i,1)=vpl(i,1)
        upl(i,2)=vpl(i,1)*vpl(i,2)
        upl(i,4)=vpl(i,1)*vpl(i,4)
        upl(i,3)=vpl(i,3)*2.5+.5*vpl(i,2)*upl(i,2)+q0*upl(i,4)

        unr(i,1)=vnr(i,1)
        unr(i,2)=vnr(i,1)*vnr(i,2)     
        unr(i,4)=vnr(i,1)*vnr(i,4)
        unr(i,3)=vnr(i,3)*2.5+.5*vnr(i,2)*unr(i,2)+q0*unr(i,4)

        unl(i,1)=vnl(i,1)
        unl(i,2)=vnl(i,1)*vnl(i,2)
        unl(i,4)=vnl(i,1)*vnl(i,4)
        unl(i,3)=vnl(i,3)*2.5+.5*vnl(i,2)*unl(i,2)+q0*unl(i,4)
      enddo

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine HMSTH_combine
!=======================================================================

      SUBROUTINE HMSTH_combine(u1,
     &                          vpr,vpl,vnr,vnl,
     &                          vr,vl,vrn,vln,
     &                          ur,ul,unr,unl,
     &                          alpha)

c-----------------------------------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imp1,q0,temperature0, ddd

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  u1(0:13000,4)

c---------------------combine together---------------------------------- 

      eps =1e-12

c---------------------right---------------------------------------------

      do k=1,4
        do i=2,imm1

          Jm1=i-1
          Jp1=i+1
          DTP= ( u1(Jm1,3) - 2.*u1(i,3) + u1(Jp1,3) )
          DTM= ( u1(Jm1,3) + 2.*u1(i,3) + u1(Jp1,3) )
          dip0 =  alpha * abs( ( eps+dtp)/(dtm+eps))
          ddp = max( (1d0 - dip0), 0. )

          ur(i,k) =  ddp * vpr(i,k)  + ( 1.- ddp ) * vr(i,k)
          ul(i,k) =  ddp * vpl(i,k)  + ( 1.- ddp ) * vl(i,k) 

        enddo
      enddo

c---------------------left----------------------------------------------

      do k=1,4
        do i=2,imm1
          
          Jm1=i-1
          Jp1=i+1
          DTP= ( u1(Jm1,3) - 2.*u1(i,3) + u1(Jp1,3) )
          DTM= ( u1(Jm1,3) + 2.*u1(i,3) + u1(Jp1,3) )
          dip0 =  alpha * abs( ( eps+dtp)/(dtm+eps))
          ddp = max( (1d0 - dip0), 0. )

          unr(i,k) = ddp * vnr(i,k)  + ( 1.- ddp ) * vrn(i,k)
          unl(i,k) = ddp * vnl(i,k)  + ( 1.- ddp ) * vln(i,k) 

        enddo
      enddo

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine HMSTH2_combine
!=======================================================================

      SUBROUTINE HMSTH2_combine(u1,
     &                          vpr,vpl,vnr,vnl,
     &                          vr,vl,vrn,vln,
     &                          ur,ul,unr,unl,
     &                          alpha)

c-----------------------------------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imp1,q0,temperature0, ddd

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  u1(0:13000,4)

c---------------------combine together---------------------------------- 

      eps =1e-12

c---------------------right---------------------------------------------

      do k=1,4
        do i=2,imm1

          Jm1=i-1
          Jp1=i+1
          DTP= ( u1(Jm1,k) - 2.*u1(i,k) + u1(Jp1,k) )
          DTM= ( u1(Jm1,k) + 2.*u1(i,k) + u1(Jp1,k) )
          dip0 =  alpha * abs( ( eps+dtp)/(dtm+eps))
          ddp = max( (1d0 - dip0), 0. )

          ur(i,k) =  ddp * vpr(i,k)  + ( 1.- ddp ) * vr(i,k)
          ul(i,k) =  ddp * vpl(i,k)  + ( 1.- ddp ) * vl(i,k) 

        enddo
      enddo

c---------------------left----------------------------------------------

      do k=1,4
        do i=2,imm1
          
          Jm1=i-1
          Jp1=i+1
          DTP= ( u1(Jm1,k) - 2.*u1(i,k) + u1(Jp1,k) )
          DTM= ( u1(Jm1,k) + 2.*u1(i,k) + u1(Jp1,k) )
          dip0 =  alpha * abs( ( eps+dtp)/(dtm+eps))
          ddp = max( (1d0 - dip0), 0. )

          unr(i,k) = ddp * vnr(i,k)  + ( 1.- ddp ) * vrn(i,k)
          unl(i,k) = ddp * vnl(i,k)  + ( 1.- ddp ) * vln(i,k) 

        enddo
      enddo

c---------------------end-----------------------------------------------

      RETURN
      END

!=======================================================================
! subroutine HMT_combine
!=======================================================================

      SUBROUTINE HMT_combine(u1,
     &                       vpr,vpl,vnr,vnl,
     &                       vr,vl,vrn,vln,
     &                       ur,ul,unr,unl)

c-----------------------------------------------------------------------

      common /cdt/dx,dt,del,r,r1,pi,imm1,imm2,imm3,imm4,imm5
     #,imp1,q0,temperature0, ddd

      real  vpr(0:13000,4) ,vpl(0:13000,4)
      real  vnr(0:13000,4) ,vnl(0:13000,4)

      real  vr(0:13000,4) ,vl(0:13000,4)
      real  vrn(0:13000,4) ,vln(0:13000,4)

      real  ur(0:13000,4) ,ul(0:13000,4)
      real  unr(0:13000,4) ,unl(0:13000,4)

      real  u1(0:13000,4)

c---------------------combine together----------------------------------
      
      epsvof = 10e-20
      
c---------------------right---------------------------------------------

      do k=1,4
        do i=2,imm1

          dvnp = ((vpl(j,k)-vpr(j,k)+epsvof)/(u1(J+1,k)-u1(J,k)+epsvof))
          dvnm = ((vpl(j,k)-vpr(j,k)+epsvof)/(u1(J,k)-u1(J-1,k)+epsvof))

          zeta = 1.0-min(dvnp,dvnm)

          if ((u1(i+1,k)-u1(i,k))*(u1(i,k)-u1(i-1,k)).GT.0) then
            ur(i,k) = ( 1.- zeta ) * vpr(i,k) + zeta * vr(i,k)
            ul(i,k) = ( 1.- zeta ) * vpl(i,k) + zeta * vl(i,k)
          else
            ur(i,k) = vpr(i,k)
            ul(i,k) = vpl(i,k)
          endif
        enddo
      enddo
      
c---------------------left----------------------------------------------

      do k=1,4
        do i=2,imm1

          dvnp = ((vnl(j,k)-vnr(j,k)+epsvof)/(u1(J+1,k)-u1(J,k)+epsvof))
          dvnm = ((vnl(j,k)-vnr(j,k)+epsvof)/(u1(J,k)-u1(J-1,k)+epsvof))
          
          zeta = 1.0-min(dvnp,dvnm)

          if ((u1(i+1,k)-u1(i,k))*(u1(i,k)-u1(i-1,k)).GT.0) then
            unr(i,k) = ( 1.- zeta ) * vnr(i,k) + zeta * vrn(i,k)
            unl(i,k) = ( 1.- zeta ) * vnl(i,k) + zeta * vln(i,k)
          else
            unr(i,k) = vnr(i,k)
            unl(i,k) = vnl(i,k)
          endif
        enddo
      enddo         

c---------------------end-----------------------------------------------

      RETURN
      END
