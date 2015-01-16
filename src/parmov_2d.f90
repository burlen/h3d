      subroutine parmov_2d
 
 
      use parameter_mod
      use MESH2D
      integer*8 count_kbq,time_begin(8),time_end(8)
      integer*8 nptotp_kbq,npart_kbq(2),np_ijk,l
      data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
      data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
      data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
      integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                 ,ixep1,iyep1,izep1,ixp1,iyp1,izp1,Storage_Error_p,Storage_Error
      double precision:: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,epsilon,xpart,ypart,zpart,r_particle
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi,twopi,myranf,d_ranf,fluxran,vy_tmp,vmag
      double precision, dimension(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: &
      bx_av,by_av,bz_av,tmp
      double precision:: v_limit,eps2,rx0,ry0,rrat,sqrr,outer_radius,q_p
      double precision:: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min,x_disp,y_disp,z_disp          &
                        ,y_disp_max_p,x_disp_max_p,y_disp_max,x_disp_max,coef_oncenter, coef_offcenter   &
                        ,norm_fac,v_pp,v_pc,v_pm,v_cp,v_cc,v_cm,v_mp,v_mc,v_mm
      integer*8:: Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p,n_spatial
!
! Yuri's smoothing - hxv 4/11/2013
!
      n_spatial      = 2
      coef_oncenter  = 1.
      coef_offcenter = smooth_coef/((3**n_spatial-1)*(1.-smooth_coef))
!
!=======================================================================
!
      call date_and_time(values=time_begin_array(:,19))
 

      Courant_Violation_p = 0
      Field_Diverge_p     = 0
      Storage_Error_p     = 0
      x_disp_max_p        = 0
      y_disp_max_p        = 0


      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

 
      twopi = 2.*acos(-1.d+00)
      d_ranf=1./1001.
      epsilon=buffer_zone
      outer_radius=dipole_sphere_radius+min(hx,hy)
      eps2=1.d-25


!     Uniform mesh - Same as is in version 5.0
!      x_dipole=hx*(i_dipole+0.5)
!      y_dipole=hy*(j_dipole+0.5)
!      z_dipole=hz*(k_dipole+0.5)

        
!      Nonuniform mesh 
       x_dipole=meshX%xc(i_dipole)
       y_dipole=meshY%xc(j_dipole)
       z_dipole=meshZ%xc(k_dipole)
 

!     Uniform mesh - Same as is in version 5.0
!      if (dt /= 0.) then
!        v_limit =min(hx/dt,hy/dt)
!      else
!        v_limit=1.d+10
!      endif
 
!     Nonuniform mesh
      hxmin=meshX%dxc(1)
      hxmax=hxmin
      do i=1,size(meshX%dxc)
        if(meshX%dxc(i) < hxmin) hxmin=meshX%dxc(i)
        if(meshX%dxc(i) > hxmax) hxmax=meshX%dxc(i)
      enddo
      hymin=meshY%dxc(1)
      hymax=hymin
      do i=1,size(meshY%dxc)
        if(meshY%dxc(i) < hymin) hymin=meshY%dxc(i)
        if(meshY%dxc(i) > hymax) hymax=meshY%dxc(i)
      enddo
      hzmin=meshZ%dxc(1)
      hzmax=hzmin
      do i=1,size(meshZ%dxc)
        if(meshZ%dxc(i) < hzmin) hzmin=meshZ%dxc(i)
        if(meshZ%dxc(i) > hzmax) hzmax=meshZ%dxc(i)
      enddo
      cell_size_min = min(hxmin,hymin,hzmin)
      v_limit=(cell_size_min/dtwci)/wpiwci
!
!=======================================================================
!
      bx_av=0.;by_av=0.;bz_av=0.
      DO K = KB-1,KE
        DO J = JB-1,JE
          DO I = 1, NX1
            bx_av(i,j,k)=0.250*( bx(i  ,j  ,k  )+bdipole_x(i  ,j  ,k  )             &
                                +bx(i+1,j  ,k  )+bdipole_x(i+1,j  ,k  )             &
                                +bx(i  ,j+1,k  )+bdipole_x(i  ,j+1,k  )             &
                                +bx(i+1,j+1,k  )+bdipole_x(i+1,j+1,k  )             &
                               )
            by_av(i,j,k)=0.250*( by(i  ,j  ,k  )+bdipole_y(i  ,j  ,k  )             &
                                +by(i+1,j  ,k  )+bdipole_y(i+1,j  ,k  )             &
                                +by(i  ,j+1,k  )+bdipole_y(i  ,j+1,k  )             &
                                +by(i+1,j+1,k  )+bdipole_y(i+1,j+1,k  )             &
                               )
            bz_av(i,j,k)=0.250*( bz(i  ,j  ,k  )+bdipole_z(i  ,j  ,k  )             &
                                +bz(i+1,j  ,k  )+bdipole_z(i+1,j  ,k  )             &
                                +bz(i  ,j+1,k  )+bdipole_z(i  ,j+1,k  )             &
                                +bz(i+1,j+1,k  )+bdipole_z(i+1,j+1,k  )             &
                               )
          enddo
        enddo
      enddo
      CALL XREALBCC_PACK_B_2D(BX_AV,BY_AV,BZ_AV,1,NX,NY,NZ)
      CALL XREALBCC_PACK_B_2D(BX   ,BY   ,BZ   ,1,NX,NY,NZ)
 
 
      nparbuf=nxmax*(nylmax+2)*(nzlmax+2)
!
!=======================================================================
!
!  initalize diagnostic variables that keep track of
!  particle number, injection, and escape
!
      deltime1 = 0.0
      deltime2 = 0.0
      nptotp=0
      nptotp_kbq=0
      npleavingp=0

      DO IS=1,NSPEC
 
  
        call date_and_time(values=time_begin_array(:,13))
        NPTOTP=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                NP=IPHEAD(ixe,iye,ize,is)
                DO WHILE (NP.NE.0)
                  NPTOTP=NPTOTP+1
                  NP=LINK(NP)
                ENDDO
            enddo
          enddo
        enddo

                                 
        call MPI_ALLREDUCE(nptotp,nptot_max,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,IERR)
        if (nptotp  == nptot_max) then 
          write(6,*) " PROCESSOR # ",MYID,", # OF PARTICLES = ",NPTOTP
          write(6,*) " MAXIMUM # OF PARTICLES ALLOWED      = ",NPLMAX
        endif
       
        call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
        IF (MYID.EQ.0) THEN
          WRITE(6,*) " IS = ",IS
          WRITE(6,*) " # OF PARTICLES BEFORE PARMOV = ",NPTOT
        ENDIF
 
 
        wmult=wspec(is)
        h=dt*qspec(is)/wmult
        hh=.5*h
 
 
      call date_and_time(values=time_end)
      clock_time1=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        IF(DT.NE.0) THEN
          npart(is) = 0
          npart_kbq(is) = 0
          DO IIZE = KB-1,KE
           DO IIYE = JB-1,JE
            DO IIXE = 1, NX1
              NP=IPHEAD(IIXE,IIYE,IIZE,IS)
!
!  begin advance of particle position and velocity
!  If dt=0, skip
!
              DO WHILE (NP.NE.0)
                L=NP
                npart_kbq(is) = npart_kbq(is)+1
                nptotp_kbq = nptotp_kbq + 1
 
 
!               Uniform mesh - Same as is in version 5.0
!                rxb=hxi*x(l)+1.9999999999999999d+00
!                ryb=hyi*y(l)+0.9999999999999999d+00
!                rzb=hzi*z(l)+0.9999999999999999d+00
!                ixb=rxb
!                iyb=ryb
!                izb=rzb
!                ixbp1 = ixb+1
!                iybp1 = iyb+1
!                izbp1 = izb+1
!                fxb=rxb-ixb
!                fyb=ryb-iyb
!                fzb=rzb-izb
 
!               Nonuniform mesh - using MESH_UNMAP
                rxb=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
                ryb=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
                rzb=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
                ixb=rxb
                iyb=ryb
                izb=rzb

                IZB=1
                fxb=rxb-ixb
                fyb=ryb-iyb
                fzb=rzb-izb
                iyb=iyb-1             ! integer index in y direction starts at 0
                izb=izb-1             ! integer index in z direction starts at 0

!                if (fxb > 0.5) then
!                  fxb = fxb -0.5
!                  ixb = ixb + 1
!                else
!                  fxb=fxb + 0.5
!                endif
!                if (fyb > 0.5) then
!                  fyb = fyb -0.5
!                  iyb = iyb + 1
!                else
!                  fyb=fyb + 0.5
!                endif
!                iyb=min(iyb,je)


                ixbp1 = ixb+1
                iybp1 = iyb+1
                izbp1 = izb+1

                w1b=(1.-fxb)*(1.-fyb)
                w2b=    fxb *(1.-fyb)
                w3b=(1.-fxb)*    fyb
                w4b=    fxb *    fyb

                bx1=bx(ixb  ,iyb  ,izb  )+bdipole_x(ixb  ,iyb  ,izb  )
                bx2=bx(ixbp1,iyb  ,izb  )+bdipole_x(ixbp1,iyb  ,izb  )
                bx3=bx(ixb  ,iybp1,izb  )+bdipole_x(ixb  ,iybp1,izb  )
                bx4=bx(ixbp1,iybp1,izb  )+bdipole_x(ixbp1,iybp1,izb  )
                by1=by(ixb  ,iyb  ,izb  )+bdipole_y(ixb  ,iyb  ,izb  )
                by2=by(ixbp1,iyb  ,izb  )+bdipole_y(ixbp1,iyb  ,izb  )
                by3=by(ixb  ,iybp1,izb  )+bdipole_y(ixb  ,iybp1,izb  )
                by4=by(ixbp1,iybp1,izb  )+bdipole_y(ixbp1,iybp1,izb  )
                bz1=bz(ixb  ,iyb  ,izb  )+bdipole_z(ixb  ,iyb  ,izb  )
                bz2=bz(ixbp1,iyb  ,izb  )+bdipole_z(ixbp1,iyb  ,izb  )
                bz3=bz(ixb  ,iybp1,izb  )+bdipole_z(ixb  ,iybp1,izb  )
                bz4=bz(ixbp1,iybp1,izb  )+bdipole_z(ixbp1,iybp1,izb  )
 
 
                bxa=w1b*bx1+w2b*bx2+w3b*bx3+w4b*bx4
                bya=w1b*by1+w2b*by2+w3b*by3+w4b*by4
                bza=w1b*bz1+w2b*bz2+w3b*bz3+w4b*bz4
 
 
!               Uniform mesh - Same as is in version 5.0
!                rxe=hxi*x(l)+1.5000000000000001d+00
!                rye=hyi*y(l)+0.5000000000000001d+00
!                rze=hzi*z(l)+0.5000000000000001d+00
!                ixe=rxe
!                iye=rye
!                ize=rze
!                IZE=1
!                ixep1 = ixe+1
!                iyep1 = iye+1
!                izep1 = ize+1
!                fxe=rxe-ixe
!                fye=rye-iye
!                fze=rze-ize

 
!               Nonuniform mesh - using MESH_UNMAP
                rxe=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
                rye=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
                rze=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
                ixe=rxe
                iye=rye
                ize=rze
                IZE=1
                fxe=rxe-ixe
                fye=rye-iye
                fze=rze-ize
                iye=iye-1             ! integer index in y direction starts at 0
                ize=ize-1             ! integer index in z direction starts at 0
 
                ixep1 = ixe+1
                iyep1 = iye+1
                izep1 = ize+1
 

                w1e=(1.-fxe)*(1.-fye)
                w2e=    fxe *(1.-fye)
                w3e=(1.-fxe)*    fye
                w4e=    fxe *    fye
 
 
                ex1=ex(ixe  ,iye  ,ize  )
                ex2=ex(ixep1,iye  ,ize  )
                ex3=ex(ixe  ,iyep1,ize  )
                ex4=ex(ixep1,iyep1,ize  )
 
 
                ey1=ey(ixe  ,iye  ,ize  )
                ey2=ey(ixep1,iye  ,ize  )
                ey3=ey(ixe  ,iyep1,ize  )
                ey4=ey(ixep1,iyep1,ize  )
 
 
                ez1=ez(ixe  ,iye  ,ize  )
                ez2=ez(ixep1,iye  ,ize  )
                ez3=ez(ixe  ,iyep1,ize  )
                ez4=ez(ixep1,iyep1,ize  )
 
 
                fox1=fox(ixe  ,iye  ,ize  )
                fox2=fox(ixep1,iye  ,ize  )
                fox3=fox(ixe  ,iyep1,ize  )
                fox4=fox(ixep1,iyep1,ize  )
 
 
                foy1=foy(ixe  ,iye  ,ize  )
                foy2=foy(ixep1,iye  ,ize  )
                foy3=foy(ixe  ,iyep1,ize  )
                foy4=foy(ixep1,iyep1,ize  )
 
 
                foz1=foz(ixe  ,iye  ,ize  )
                foz2=foz(ixep1,iye  ,ize  )
                foz3=foz(ixe  ,iyep1,ize  )
                foz4=foz(ixep1,iyep1,ize  )
 
 
                exa=w1e* ex1+w2e* ex2+w3e* ex3+w4e* ex4&
                   +w1e*fox1+w2e*fox2+w3e*fox3+w4e*fox4
                eya=w1e* ey1+w2e* ey2+w3e* ey3+w4e* ey4&
                   +w1e*foy1+w2e*foy2+w3e*foy3+w4e*foy4
                eza=w1e* ez1+w2e* ez2+w3e* ez3+w4e* ez4&
                   +w1e*foz1+w2e*foz2+w3e*foz3+w4e*foz4
 
 
                ff=2./(1.+hh*hh*(bxa**2+bya**2+bza**2))
                vex=vx(l)+exa*hh
                vey=vy(l)+eya*hh
                vez=vz(l)+eza*hh
                p2xs=vex+(vey*bza-vez*bya)*hh
                p2ys=vey+(vez*bxa-vex*bza)*hh
                p2zs=vez+(vex*bya-vey*bxa)*hh

                vx(l)=vex+ff*(p2ys*bza-p2zs*bya)*hh+exa*hh
                vy(l)=vey+ff*(p2zs*bxa-p2xs*bza)*hh+eya*hh
                vz(l)=vez+ff*(p2xs*bya-p2ys*bxa)*hh+eza*hh

                x_disp = dt*vx(l)
                y_disp = dt*vy(l)

                x(l)=x(l)+ x_disp
                y(l)=y(l)+ y_disp

                if ( abs(x_disp/meshX%dxn(ixep1)) > 1.0 .or.                                               &
                     abs(y_disp/meshY%dxn(iyep1)) > 1.0) Courant_Violation_p = Courant_Violation_p + 1
                x_disp_max_p = max(x_disp_max_p,abs(x_disp)/meshX%dxn(ixep1))
                y_disp_max_p = max(y_disp_max_p,abs(y_disp)/meshY%dxn(iyep1))

                NP=LINK(NP)
              ENDDO
            ENDDO
           ENDDO
          ENDDO
        ENDIF

        call MPI_ALLREDUCE(x_disp_max_p,x_disp_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(y_disp_max_p,y_disp_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
        if (myid == 0) then
          write(6,*) " maximum x-displacement/dx = ",x_disp_max
          write(6,*) " maximum y-displacement/dy = ",y_disp_max
        endif

        call MPI_ALLREDUCE(Courant_Violation_p,Courant_Violation,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

        if (Courant_Violation /= 0) then
           if (myid == 0) write(6,*) "Particle displacements exceed cell size ",Courant_Violation," times"
           call MPI_FINALIZE(IERR)
           STOP
        endif



        call date_and_time(values=time_end_array(:,13))
        call accumulate_time_difference(time_begin_array(1,13),time_end_array(1,13),time_elapsed(13))
 
 
        call date_and_time(values=time_begin_array(:,14))
 
 
        ICOUNT=0                  !ICOUNT records how many times the particle
  11    ICOUNT=ICOUNT+1           !  exchange routine has been run. If there               
        irepeatp=0                !  are no fast particles, it will be called 
        do ipe=0,numprocs-1       !  once.
          nsendp(ipe)=0
          nrecvp(ipe)=0
        enddo
        ipleft (is)=0
        iprite (is)=0
        ipsendleft(is)=0
        ipsendrite(is)=0
 
 
        ipsend(is)=0
        ipsendtop(is)=0
        ipsendbot(is)=0
        ipsendlefttop(is)=0
        ipsendleftbot(is)=0
        ipsendritetop(is)=0
        ipsendritebot(is)=0
 

        nescape(is)=0
        nescape_yz(is)=0
        nescape_zy(is)=0
        nescape_xy(is)=0
        nescape_yx(is)=0
        nescape_zx(is)=0
        nescape_xz(is)=0
        DO IZE = KB-1,KE
         DO IYE = JB-1,JE
          DO IXE = 1,NX1
            iptemp(ixe,iye,ize,is)=0
!
!=======================================================================
!
!     mark particles that need to be sent to other processors
!
            NP=IPHEAD(IXE,IYE,IZE,IS)
            DO WHILE (NP.NE.0)
              xpart=x(np)
              ypart=y(np)
              zpart=z(np)
              l=np
              rx0=xpart-x_dipole
              ry0=ypart-y_dipole
              r_particle=sqrt((xpart-x_dipole)**2+(ypart-y_dipole)**2)
              sqrr=r_particle+eps2
              if ((vx(l)**2+vy(l)**2+vz(l)**2)*wpiwci**2 >= 400.) then
                iphead(ixe,iye,ize,is)=link(np)
                link(np)=ipstore
                ipstore=np
                np=iphead(ixe,iye,ize,is)
                goto 15
              endif
!
!
              if (r_particle < dipole_sphere_radius) then
 
                if (absorbing_dipole) then
                   iphead(ixe,iye,ize,is)=link(np)
                   link(np)=ipstore
                   ipstore=np
                   np=iphead(ixe,iye,ize,is)
                else
                  rrat = (2.0 *dipole_sphere_radius -sqrr) /sqrr
                  rrat=min(rrat,outer_radius/sqrr)
                  x(l) = rx0 *rrat + x_dipole
                  y(l) = ry0 *rrat + y_dipole
                  vx(np)=-vx(np)
                  vy(np)=-vy(np)
                  x (np)=x(np)+vx(np)*dt
                  y (np)=y(np)+vy(np)*dt
                  iphead(ixe,iye,ize,is)=link(np)
                  link(np)=iptemp(ixe,iye,ize,is)
                  iptemp(ixe,iye,ize,is)=np
                  np=iphead(ixe,iye,ize,is)
                endif
                goto 15

              endif
 
! first take care of particle leakage through Y and Z sides of box (reflect)
 
!              if (y(l) <  epsilon) then

!                  y(l) = 2.*epsilon - y(l)
!                  vy(l)= 2.*vybar(is)-vy(l)
!                  ypart= y(l)

!              endif
 

!              if (y(l) >  ymax-epsilon) then

!                  y(l) = 2.*(ymax-epsilon) - y(l)
!                  vy(l)= 2.*vybar(is)-vy(l)
!                  ypart= y(l)

!              endif
 

                  if       (y(l) <   epsilon     ) then
                    
                    y(l) = myranf()*epsilon

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vy_tmp=vpar(is)*u_array_injection(iv,jv)
                    vy(l)=+vy_tmp
                    vmag=sqrt(-log(1.-.999999*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*cos(th)
                    vz(l)=          vmag*vper(is)*sin(th)
                    x (l)=   myranf()*xmax
                    z (l)=zb+myranf()*(ze-zb)

                  endif

                  if       (y(l) >   ymax-epsilon   ) then

                    y(l) = ymax-myranf()*epsilon

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vy_tmp=vpar(is)*u_array_injection(iv,jv)
                    vy(l)=-vy_tmp
                    vmag=sqrt(-log(1.-.999999*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*cos(th)
                    vz(l)=          vmag*vper(is)*sin(th)
                    x (l)=     myranf()*xmax
                    z (l)=zb  +myranf()*(ze-zb)

                  endif

 
!  identify particles that have crossed outside the physical domain 
!  of simulation
 
              if (    xpart.lt.0..or.xpart.gt.xmax&
                  .or.ypart.lt.0..or.ypart.gt.ymax&       
                  .or.zpart.lt.0..or.zpart.gt.zmax) then
 
!  if particle np is outside physcial domain, set link(np) to previous
!  head of "empty" list, and set ipstore to np. Thus, np becomes new
!  head of "empty" list. 
 
                if (xpart.lt.0.) then
                  nescape_yz(is)=nescape_yz(is)+1
                else if (xpart.gt.xmax) then
                  nescape_zy(is)=nescape_zy(is)+1
                else if (ypart.lt.0.) then
                  nescape_xz(is)=nescape_xz(is)+1
                else if (ypart.gt.ymax) then
                  nescape_zx(is)=nescape_zx(is)+1
                else if (zpart.lt.0.) then
                  nescape_xy(is)=nescape_xy(is)+1
                else 
                  nescape_yx(is)=nescape_yx(is)+1
                endif
 
 
                npleavingp=npleavingp+1
                nescape(is) = nescape(is) + 1
                iphead(ixe,iye,ize,is)=link(np)
                link(np)=ipstore
                ipstore=np
                np=iphead(ixe,iye,ize,is)
 
 
              else
 
 
                if ((zpart.le.ze.and.zpart.ge.zb).and.(ypart.le.ye.and.ypart.ge.yb)) then
                  iphead(ixe,iye,ize,is)=link(np)
                  link(np)=iptemp(ixe,iye,ize,is)
                  iptemp(ixe,iye,ize,is)=np
                  np=iphead(ixe,iye,ize,is)
                ELSE
                  if (ypart.le.ye.and.ypart.ge.yb) then
                    iye_cc=jb
                  else
                    if (ypart.gt.ye) then
                      iye_cc=je+1 
                    else
                      iye_cc=jb-1 
                    endif
                  endif
                  if (zpart.le.ze.and.zpart.ge.zb) then
                    ize_cc=kb
                  else
                    if (zpart.gt.ze) then
                      ize_cc=ke+1 
                    else
                      ize_cc=kb-1 
                    endif
                  endif
                  nsendp(idmap_yz(iye_cc,ize_cc))=nsendp(idmap_yz(iye_cc,ize_cc))+1
                  iphead(ixe,iye,ize,is)=link(np)
                  link(np)=ipsend(is)
                  ipsend(is)=np
                  np=iphead(ixe,iye,ize,is)
                ENDIF
              ENDIF
 15           CONTINUE
            ENDDO
 
!  set iphead(is) to be the head of the list for species is, and reset
!  iptemp(is) to zero
 
            iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
            iptemp(ixe,iye,ize,is)=0

          ENDDO
        ENDDO
      ENDDO

!        if (dt /= 0.) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF
!
!************************************************************************
!
!     exchange data among processes and compute to see how many
!     particles each process has to send to, and receive from, other
!     processes
!
        call MPI_SENDRECV(nsendp(NBRLEFT   ),1,MPI_INTEGER8,NBRLEFT   ,0,&
                          nrecvp(NBRRITE   ),1,MPI_INTEGER8,NBRRITE   ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRRITE   ),1,MPI_INTEGER8,NBRRITE   ,0,&
                          nrecvp(NBRLEFT   ),1,MPI_INTEGER8,NBRLEFT   ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRTOP    ),1,MPI_INTEGER8,NBRTOP    ,0,&
                          nrecvp(NBRBOT    ),1,MPI_INTEGER8,NBRBOT    ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRBOT    ),1,MPI_INTEGER8,NBRBOT    ,0,&
                          nrecvp(NBRTOP    ),1,MPI_INTEGER8,NBRTOP    ,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRLEFTTOP),1,MPI_INTEGER8,NBRLEFTTOP,0,&
                          nrecvp(NBRRITEBOT),1,MPI_INTEGER8,NBRRITEBOT,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRRITEBOT),1,MPI_INTEGER8,NBRRITEBOT,0,&
                          nrecvp(NBRLEFTTOP),1,MPI_INTEGER8,NBRLEFTTOP,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRRITETOP),1,MPI_INTEGER8,NBRRITETOP,0,&
                          nrecvp(NBRLEFTBOT),1,MPI_INTEGER8,NBRLEFTBOT,0,&
                          mpi_comm_world,status,ierr)
        call MPI_SENDRECV(nsendp(NBRLEFTBOT),1,MPI_INTEGER8,NBRLEFTBOT,0,&
                          nrecvp(NBRRITETOP),1,MPI_INTEGER8,NBRRITETOP,0,&
                          mpi_comm_world,status,ierr)
!
!***********************************************************************
!
          nsendtotp=sum(nsendp)
          nrecvtotp=sum(nrecvp)
          call MPI_ALLREDUCE(nsendtotp,nsendtot,1,MPI_INTEGER8,MPI_SUM,&
                             MPI_COMM_WORLD,IERR)
          call MPI_ALLREDUCE(nrecvtotp,nrecvtot,1,MPI_INTEGER8,MPI_SUM,&
                             MPI_COMM_WORLD,IERR)
          if (myid.eq.0) then
            write(6,*) " FINISHED COMPILING LISTS "
            write(6,*) " # OF PARTICLES TO BE SENT     = ",NSENDTOT
            write(6,*) " # OF PARTICLES TO BE RECEIVED = ",NRECVTOT
          endif
          if (NSENDTOT.NE.NRECVTOT) THEN
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF

!        
!       Check to see if each processor has enough particle storage
!       to handle incoming particles
!        
        IF (NPTOTP+NRECVTOTP < NPLMAX) THEN
          EXIT_CODE_P = 0
        ELSE 
          EXIT_CODE_P = 1 
          write(6,*) " PROCESSOR # ",myid," RAN OUT OF PARTICLE STORAGE"
        ENDIF 
        call MPI_ALLREDUCE(EXIT_CODE_P,EXIT_CODE,1,MPI_INTEGER8,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
        IF (EXIT_CODE /= 0) THEN 
           if (myid == 0) write(6,*) "3DHybrid is stopped"
          CALL MPI_FINALIZE(IERR)
          STOP
        ENDIF
!
!***********************************************************************
!
!     exchange particles with other processors
!
        nsendactualp=0
        nrecvactualp=0
        do irepeat=1,4
          if (isendid(irepeat).eq.1) then
            NP=IPSEND(IS)
            DO WHILE (NP.NE.0)
              nsendactualp=nsendactualp+1

!             Uniform mesh - Same as is in ver 5.0
!              ixe=hxi*x(np)   +1.5000000000000001d+00
!              iye=hyi*y(np)   +0.5000000000000001d+00
!              ize=hzi*z(np)   +0.5000000000000001d+00

!             Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
              rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
              rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
              ixe=rxe
              iye=rye
              ize=rze
              IZE=1
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0
 
              ypart=y(np)
              zpart=z(np)
              if (ypart.le.ye.and.ypart.ge.yb) then
                iye_cc=jb 
              else
                if (ypart.gt.ye) then
                  iye_cc=je+1 
                else
                  iye_cc=jb-1 
                endif
              endif
              if (zpart.le.ze.and.zpart.ge.zb) then
                ize_cc=kb 
              else
                if (zpart.gt.ze) then
                  ize_cc=ke+1 
                else
                  ize_cc=kb-1 
                endif
              endif
              pdata(1)=x(np)
              pdata(2)=y(np)
              pdata(3)=z(np)
              pdata(4)=vx(np)
              pdata(5)=vy(np)
              pdata(6)=vz(np)
              pdata(7)=qp(np)
              call MPI_SEND(pdata,7,MPI_DOUBLE_PRECISION,&
                             idmap_yz(iye_cc,ize_cc),IT,MPI_COMM_WORLD,IERR)
              ipsend(is)=link(np)
              link(np)=ipstore
              ipstore=np
              np=ipsend(is)
            ENDDO
          else
            nprecvtmp=0
            do itmp=1,4
              ipe=irecvid(itmp,irepeat)
              if (ipe.ne.-1) then
                nprecvtmp=nprecvtmp+nrecvp(irecvid(itmp,irepeat))
                nrecvp(irecvid(itmp,irepeat))=0
              endif
            enddo
            do ii=1,nprecvtmp
              nrecvactualp=nrecvactualp+1
              nprecv=ipstore
              call MPI_RECV(pdata,7,MPI_DOUBLE_PRECISION,&
                            MPI_ANY_SOURCE,IT  ,&
                            MPI_COMM_WORLD,STATUS2,IERR)

              if (ipstore == 0) then
                Storage_Error_p = 1
              else

                x(nprecv)=pdata(1)
                y(nprecv)=pdata(2)
                z(nprecv)=pdata(3)
                vx(nprecv)=pdata(4)
                vy(nprecv)=pdata(5)
                vz(nprecv)=pdata(6)
                qp(nprecv)=pdata(7)

!               Nonuniform mesh - using MESH_UNMAP
                rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000d+00
                rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000d+00
                rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000d+00
                ixe=rxe 
                iye=rye 
                ize=rze 
                IZE=1
                iye=iye-1             ! integer index in y direction starts at 0
                ize=ize-1             ! integer index in z direction starts at 0
             

                ipstore=link(nprecv)

                if ((ixe > nx+1  .or. ixe < 1 ) .or. (iye > je+1    .or. iye < jb-1)) then
                  Field_Diverge_p = 1
                  ixe = min(max(ixe,1 ),nx+1)
                  iye = min(max(iye,jb-1),je+1  )
                endif

                link(nprecv)=iphead(ixe,iye,ize,is)
                iphead(ixe,iye,ize,is)=nprecv

              endif

            enddo
          endif
        enddo

        call MPI_ALLREDUCE(Storage_Error_p,Storage_Error,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
        if (Storage_Error /= 0) then
           if (myid == 0) then
             write(6,*)" "
             write(6,*)" "
             write(6,*) "Particle storage allocated is exceeded."
             write(6,*) "3DHybrid is stopped"
             write(6,*)" "
             write(6,*)" "
           endif
           call MPI_FINALIZE(IERR)
           STOP
        endif

        call MPI_ALLREDUCE(Field_Diverge_p,Field_Diverge,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
        if (Field_Diverge /= 0) then
           if (myid == 0) then
             write(6,*)" "
             write(6,*)" "
             write(6,*) "Field Solver Diverges"
             write(6,*) "3DHybrid is stopped"
             write(6,*)" "
             write(6,*)" "
           endif
           call MPI_FINALIZE(IERR)
           STOP
        endif

!
!=======================================================================
!
          call MPI_ALLREDUCE(nsendactualp,nsendactual,1,MPI_INTEGER8,&
                             MPI_SUM,MPI_COMM_WORLD,IERR)
          call MPI_ALLREDUCE(nrecvactualp,nrecvactual,1,MPI_INTEGER8,&
                             MPI_SUM,MPI_COMM_WORLD,IERR)
          if (myid.eq.0) then
            write(6,*) " FINISHED EXCHANGING PARTICLES "
            write(6,*) " # OF PARTICLES       SENT     = ",NSENDACTUAL
            write(6,*) " # OF PARTICLES       RECEIVED = ",NRECVACTUAL
          endif

!        if (dt /= 0.) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF
!
!=======================================================================
!
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1
                NP=IPHEAD(ixe,iye,ize,is)
                DO WHILE (NP.NE.0)
                  NPTOTP=NPTOTP+1
                  NP=LINK(NP)
                ENDDO
             enddo
           enddo
         enddo
         call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,&
                            MPI_COMM_WORLD,IERR)
         IF (MYID.EQ.0) THEN
           WRITE(6,*) " IS = ",IS
           WRITE(6,*) " # OF PARTICLES AFTER  PARMOV = ",NPTOT
         ENDIF
 999     CONTINUE
 
 
        call date_and_time(values=time_end_array(:,14))
        call accumulate_time_difference(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))
 
 
        call date_and_time(values=time_begin_array(:,15))
!
!*****************************************************************************
!
        if (testorbt) goto 10
!
!  Go through particle list and collect moments
!
        DO IIZ=KB-1,KE
          DO IIY=JB-1,JE
            DO IIX=1,NX1
              np_ijk=0
              np=iphead(iix,iiy,iiz,is)
              do while (np.ne.0)
                np_ijk=np_ijk+1
                nptotp=nptotp+1           !count particles
                npart(is) = npart(is) + 1 !count particles in each species
                L=NP

!               ver8.0
!
                q_p=qp(np)
                if (q_p <= 0.) write(6,*) " q_p =",q_p,iix,iiy,iiz
!

!               Uniform mesh - Same as in ver 5.0
!                rx=hxi*x(l)+1.5000000000000001d+00
!                ry=hyi*y(l)+0.5000000000000001d+00
!                rz=hzi*z(l)+0.5000000000000001d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

         
!              Nonuniform mesh - using MESH_UNMAP
                rx=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
                ry=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
                rz=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
                ix=rx 
                iy=ry 
                iz=rz 
                fx=rx-ix 
                fy=ry-iy 
                fz=rz-iz 
                iy=iy-1             ! integer index in y direction starts at 0
                iz=iz-1             ! integer index in z direction starts at 0
 
                ixp1 = ix+1
                iyp1 = iy+1
                izp1 = iz+1
 

                w1=q_p*(1.-fx)*(1.-fy)
                w2=q_p*fx     *(1.-fy)
                w3=q_p*(1.-fx)*fy     
                w4=q_p*fx     *fy     
 
                dns(ix  ,iy  ,iz  ,is)=dns(ix  ,iy  ,iz  ,is)+w1
                dns(ix+1,iy  ,iz  ,is)=dns(ix+1,iy  ,iz  ,is)+w2
                dns(ix  ,iy+1,iz  ,is)=dns(ix  ,iy+1,iz  ,is)+w3
                dns(ix+1,iy+1,iz  ,is)=dns(ix+1,iy+1,iz  ,is)+w4
 
                vxs(ix  ,iy  ,iz  ,is)=vxs(ix  ,iy  ,iz  ,is)+w1*vx(l)
                vxs(ix+1,iy  ,iz  ,is)=vxs(ix+1,iy  ,iz  ,is)+w2*vx(l)
                vxs(ix  ,iy+1,iz  ,is)=vxs(ix  ,iy+1,iz  ,is)+w3*vx(l)
                vxs(ix+1,iy+1,iz  ,is)=vxs(ix+1,iy+1,iz  ,is)+w4*vx(l)
 
                vys(ix  ,iy  ,iz  ,is)=vys(ix  ,iy  ,iz  ,is)+w1*vy(l)
                vys(ix+1,iy  ,iz  ,is)=vys(ix+1,iy  ,iz  ,is)+w2*vy(l)
                vys(ix  ,iy+1,iz  ,is)=vys(ix  ,iy+1,iz  ,is)+w3*vy(l)
                vys(ix+1,iy+1,iz  ,is)=vys(ix+1,iy+1,iz  ,is)+w4*vy(l)
 
                vzs(ix  ,iy  ,iz  ,is)=vzs(ix  ,iy  ,iz  ,is)+w1*vz(l)
                vzs(ix+1,iy  ,iz  ,is)=vzs(ix+1,iy  ,iz  ,is)+w2*vz(l)
                vzs(ix  ,iy+1,iz  ,is)=vzs(ix  ,iy+1,iz  ,is)+w3*vz(l)
                vzs(ix+1,iy+1,iz  ,is)=vzs(ix+1,iy+1,iz  ,is)+w4*vz(l)
 
 
                np=link(np)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
 10     kspc=is
 
 
        call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape(is),nescape_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape_yz(is),nescape_yz_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape_zy(is),nescape_zy_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape_xy(is),nescape_xy_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape_yx(is),nescape_yx_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape_zx(is),nescape_zx_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nescape_xz(is),nescape_xz_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
        deltime2 = deltime2 + real(clock_time1-clock_time)
 
 
        call XREAL_2D(DNS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL_2D(VXS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL_2D(VYS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL_2D(VZS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREALBCC_2D(DNS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC_2D(VXS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC_2D(VYS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC_2D(VZS(1,jb-1,kb-1,is),1,NX,NY,NZ)
!
! Up to this point dns, vxs, vys are still the accumulated charge and current, not charge density and current density.
! Apply Yuri's smoothing scheme
!
        if (smoothing .and. smooth_coef /= 0.) then
          if (myid == 0) write(6,*) " smooth_coef = ",smooth_coef
          DO  IIZ = KB,KE
            DO IIY=JB,JE
              DO IIX=2,NX1
                v_pc = meshX%dxc(iix  )*meshY%dxc(iiy+2)*meshZ%dxc(iiz+1)
                v_cc = meshX%dxc(iix  )*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1)
                v_mc = meshX%dxc(iix  )*meshY%dxc(iiy  )*meshZ%dxc(iiz+1)
                v_pp = meshX%dxc(iix  )*meshY%dxc(iiy+2)*meshZ%dxc(iiz+2)
                v_cp = meshX%dxc(iix  )*meshY%dxc(iiy+1)*meshZ%dxc(iiz+2)
                v_mp = meshX%dxc(iix  )*meshY%dxc(iiy  )*meshZ%dxc(iiz+2)
                v_pm = meshX%dxc(iix  )*meshY%dxc(iiy+2)*meshZ%dxc(iiz  )
                v_cm = meshX%dxc(iix  )*meshY%dxc(iiy+1)*meshZ%dxc(iiz  )
                v_mm = meshX%dxc(iix  )*meshY%dxc(iiy  )*meshZ%dxc(iiz  )
                norm_fac=coef_oncenter*v_cc+coef_offcenter*(v_pp+v_pc+v_pm+v_cp+v_cm+v_mp+v_mc+v_mm)
                tmp(iix,iiy,iiz) =  coef_oncenter  *    dns(iix  ,iiy  ,iiz  ,is)        &
                                  + coef_offcenter *  ( dns(iix  ,iiy  ,iiz+1,is)        &
                                                       +dns(iix  ,iiy  ,iiz-1,is)        &
                                                       +dns(iix  ,iiy+1,iiz+1,is)        &
                                                       +dns(iix  ,iiy+1,iiz  ,is)        &
                                                       +dns(iix  ,iiy+1,iiz-1,is)        &
                                                       +dns(iix  ,iiy-1,iiz+1,is)        &
                                                       +dns(iix  ,iiy-1,iiz  ,is)        &
                                                       +dns(iix  ,iiy-1,iiz-1,is)         )
                dns(iix,iiy,iiz,is) = tmp(iix,iiy,iiz) / norm_fac

                tmp(iix,iiy,iiz) =  coef_oncenter  *    vxs(iix  ,iiy  ,iiz  ,is)        &
                                  + coef_offcenter *  ( vxs(iix  ,iiy  ,iiz+1,is)        &
                                                       +vxs(iix  ,iiy  ,iiz-1,is)        &
                                                       +vxs(iix  ,iiy+1,iiz+1,is)        &
                                                       +vxs(iix  ,iiy+1,iiz  ,is)        &
                                                       +vxs(iix  ,iiy+1,iiz-1,is)        &
                                                       +vxs(iix  ,iiy-1,iiz+1,is)        &
                                                       +vxs(iix  ,iiy-1,iiz  ,is)        &
                                                       +vxs(iix  ,iiy-1,iiz-1,is)         )
                vxs(iix,iiy,iiz,is) = tmp(iix,iiy,iiz) / norm_fac

                tmp(iix,iiy,iiz) =  coef_oncenter  *    vys(iix  ,iiy  ,iiz  ,is)        &
                                  + coef_offcenter *  ( vys(iix  ,iiy  ,iiz+1,is)        &
                                                       +vys(iix  ,iiy  ,iiz-1,is)        &
                                                       +vys(iix  ,iiy+1,iiz+1,is)        &
                                                       +vys(iix  ,iiy+1,iiz  ,is)        &
                                                       +vys(iix  ,iiy+1,iiz-1,is)        &
                                                       +vys(iix  ,iiy-1,iiz+1,is)        &
                                                       +vys(iix  ,iiy-1,iiz  ,is)        &
                                                       +vys(iix  ,iiy-1,iiz-1,is)         )
                vys(iix,iiy,iiz,is) = tmp(iix,iiy,iiz) / norm_fac

                tmp(iix,iiy,iiz) =  coef_oncenter  *    vzs(iix  ,iiy  ,iiz  ,is)        &
                                  + coef_offcenter *  ( vzs(iix  ,iiy  ,iiz+1,is)        &
                                                       +vzs(iix  ,iiy  ,iiz-1,is)        &
                                                       +vzs(iix  ,iiy+1,iiz+1,is)        &
                                                       +vzs(iix  ,iiy+1,iiz  ,is)        &
                                                       +vzs(iix  ,iiy+1,iiz-1,is)        &
                                                       +vzs(iix  ,iiy-1,iiz+1,is)        &
                                                       +vzs(iix  ,iiy-1,iiz  ,is)        &
                                                       +vzs(iix  ,iiy-1,iiz-1,is)         )
                vzs(iix,iiy,iiz,is) = tmp(iix,iiy,iiz) / norm_fac
              ENDDO
            ENDDO
          ENDDO
          call XREALBCC_2D(DNS(1,jb-1,kb-1,is),1,NX,NY,NZ)
          call XREALBCC_2D(VXS(1,jb-1,kb-1,is),1,NX,NY,NZ)
          call XREALBCC_2D(VYS(1,jb-1,kb-1,is),1,NX,NY,NZ)
          call XREALBCC_2D(VZS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        else
          DO  IIZ = KB-1,KE+1
            DO IIY=JB-1,JE+1
              DO IIX=1,NX2
                dns(iix,iiy,iiz,is) = dns(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
                vxs(iix,iiy,iiz,is) = vxs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
                vys(iix,iiy,iiz,is) = vys(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
                vzs(iix,iiy,iiz,is) = vzs(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              ENDDO
            ENDDO
          ENDDO
        endif

 
        call date_and_time(values=time_end_array(:,15))
        call accumulate_time_difference(time_begin_array(1,15),time_end_array(1,15),time_elapsed(15))
      ENDDO

!        if (dt /= 0.) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF
 
!    end of main particle loop
 
      call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(npleavingp,npleaving,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

      if(prntinfo) then
 
 
      call date_and_time(values=time_end)
      clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 
 
      endif
      if(myid.eq.0.and.prntinfo) then
        do is=1,nspec
          if(is.eq.1) then
            write(6,*)
            write(6,*) " it = ",it
            write(6,*) "  species #    ninj     nescape     ntot  "
          endif
          write(6,1000) is,ninj_global(is),nescape_global(is)&
                       ,npart_global(is)
 
 
          write(6,1000) is,ninj_global(is),nescape_yz_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_zy_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_xy_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_yx_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_xz_global(is)&
                       ,npart_global(is)
          write(6,1000) is,ninj_global(is),nescape_zx_global(is)&
                       ,npart_global(is)
 
 
 1002     format(3i6,5x,i8,2x,i8,2x,i10)
 1000     format(5x,i2,5x,i8,2x,i8,2x,i10)
        enddo
        if (nspec >= 2) then
          write(6,1005) ninj_global(1)+ninj_global(2)     &
             ,nescape_global(1)+nescape_global(2)   &
             ,npart_global(1)+npart_global(2)
        else
          write(6,1005) ninj_global(1)     &
             ,nescape_global(1)    &
             ,npart_global(1)
        endif
 1005   format(5x,'sum',4x,i8,2x,i8,2x,i10)
 1006   format(2i6,4x,'sum',4x,i8,2x,i8,2x,i10)
      endif
      do is=1,nspec
        ninj(is) = 0
        ninj_global(is) = 0
      enddo
 
        call date_and_time(values=time_end_array(:,19))
        call accumulate_time_difference(time_begin_array(1,19),time_end_array(1,19),time_elapsed(19))
 

!        if (dt /= 0.) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF
      return
      end
