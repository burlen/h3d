!
!***********************************************************************
!
! ADDED BY HXV  1/09/03
!
      subroutine particle_newinject_linked_list(is)
!
!=======================================================================
!
      use parameter_mod
      use MESH2D
      use functions_f90
      use erf_mod
      integer*8, intent(in):: is
      integer*8:: ninjl,ninjr,n_inject_left_global,n_inject_rite_global      &
                 ,iv,jv,ixe,iye,ize,n_inject_l_est,n_inject_r_est            &
                 ,n_inject_left,n_inject_rite,ixe_v,iye_v,ize_v
      double precision:: vx_tmp,argum,vtherm,fluxran,residue,direction       &
                        ,fluxran_old,twopi,myranf,d_v,d_ranf,q_p,dqleft      &
                        ,dqrite,dtxi,dtyi,dtzi,rxe,rye,rze,particle_weight   &
                        ,xxl,xxr,xxl_global,xxr_global,y_p_logical,z_p_logical
      logical:: start_from_beginning


      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

!
!=======================================================================
!
      start_from_beginning=.true.
      it_inject=1
      sqrtpi=sqrt(pi)
      twopi=2.*pi
      d_v=1./16.
      d_ranf=1./1001.
      particle_weight = hx*hy*hz*dfac(is)*frac(is)
!
!=======================================================================
!
!     inject the particles into the simulation domain from left
!     and right boundaries
!
          dqleft=+0.5*dt*it_inject                                &
                 *( vpar(is)*exp(-(vxbar(is)/vpar(is))**2)/sqrtpi &
                   +vxbar(is)*(1.+erf( dble(vxbar(is)/vpar(is)))) &
                  )
          dqrite=+0.5*dt*it_inject                                &
                 *( vpar(is)*exp(-(vxbar(is)/vpar(is))**2)/sqrtpi &
                   -vxbar(is)*(1.+erf(-dble(vxbar(is)/vpar(is)))) &
                  )
!          dqrite=0.0
          qleft(is)=qleft(is)+dqleft                   
          qrite(is)=qrite(is)+dqrite

          if (myid == 0 .and. is == 1) then
            write(6,*) " qleft = ",qleft(is)
            write(6,*) " qrite = ",qrite(is)
           endif

          xxl            = 0.
          xxr            = 0.
          do j=jb,je
            do k=kb,ke
              xxl           = xxl                                                                        &
                                  +(qleft(is)/meshX%dxc(2   ))/dfac(is)
              xxr           = xxr                                                                        &
                                  +(qrite(is)*meshX%dxc(nx-1))/dfac(is)
            enddo
          enddo


          if (myid == 0 .and. is == 1) then
            write(6,*) " xxl   = ",xxl          
            write(6,*) " xxr   = ",xxr      
           endif


!         Uniform mesh - Same as is in version 5.0
!          q_p=fxsho*xmax/(npx(is)*npy(is)*npz(is))

!         Nonuniform mesh
          q_p=fxsho*xmax/(volume_fraction*npx(is)*npy(is)*npz(is)*npes)

          if (uniform_loading_in_logical_grid) then
            call MPI_ALLREDUCE(xxl           ,xxl_global   ,1,MPI_DOUBLE_PRECISION,&
                               MPI_SUM,MPI_COMM_WORLD,IERR)
            call MPI_ALLREDUCE(xxr           ,xxr_global   ,1,MPI_DOUBLE_PRECISION,&
                               MPI_SUM,MPI_COMM_WORLD,IERR)
            ninjl = xxl
            ninjr = xxr
            n_inject=ninjl+ninjr
            n_inject_left_global = xxl_global
            n_inject_rite_global = xxr_global
            if (myid == 0) then
              write(6,*) " # of particle injected from left  boundary = ",n_inject_left_global
              write(6,*) " # of particle injected from right boundary = ",n_inject_rite_global
              write(6,*) " dqleft = ",dqleft
              write(6,*) " qleft  = ", qleft(is)
              write(6,*) " dt     = ",dt
              write(6,*) " dqrite = ",dqrite
              write(6,*) " qrite  = ", qrite(is)
            endif
          else
            ninjl=qleft(is)/q_p
            ninjr=qrite(is)/q_p
            qinjl=ninjl*q_p
            qinjr=ninjr*q_p
            n_inject=ninjl+ninjr
            call MPI_ALLREDUCE(ninjl,n_inject_left_global,1,MPI_INTEGER8,&
                               MPI_SUM,MPI_COMM_WORLD,IERR)
            call MPI_ALLREDUCE(ninjr,n_inject_rite_global,1,MPI_INTEGER8,&
                             MPI_SUM,MPI_COMM_WORLD,IERR)
            if (myid == 0) then
              write(6,*) " # of particle injected from left  boundary = ",n_inject_left_global
              write(6,*) " # of particle injected from right boundary = ",n_inject_rite_global
              write(6,*) " dqleft = ",dqleft
              write(6,*) " qleft  = ", qleft(is)
              write(6,*) " dt     = ",dt
              write(6,*) " q_p    = ",q_p
            endif
          endif
!
!
          if (uniform_loading_in_logical_grid) then
            qinjl=0.
            qinjr=0.
          endif

          jv=int((vxbar(is)/vpar(is))/d_v)
	  if (myid == 0 .and. jv > 1000)               &
          write(6,*) " ERROR in injection.f90 - Increase the size of u_array_injection."
          jv=max(0,min(1000,jv))
          DO WHILE (ninjl.ne.0)
            np=ipstore
            if (np >= nplmax .or. np <= 0) then
              write(6,*) "myid = ",myid
              write(6,*) " np = ",np
              write(6,*) " nplmax = ",nplmax
            endif
            fluxran=myranf()
            iv=int((fluxran/d_ranf)+0.5)
            iv=max(1,min(1001,iv))
            vx_tmp=vpar(is)*u_array_injection(iv,jv)
            vx(np)=+vx_tmp
            vmag=sqrt(-log(1.-.999999*myranf()))
            th=twopi*myranf()
            vy(np)=vmag*vper(is)*cos(th)
            vz(np)=vmag*vper(is)*sin(th)
            x(np)=+vx(np)*dt*it_inject*myranf()

            if (uniform_loading_in_logical_grid) then
              Y_P_LOGICAL  = YB_LOGICAL+(YE_LOGICAL-YB_LOGICAL)*myranf()
              Z_P_LOGICAL  = ZB_LOGICAL+(ZE_LOGICAL-ZB_LOGICAL)*myranf()
              iye_v=dtyi*y_p_logical+1.00000000000d+00
              ize_v=dtzi*z_p_logical+1.00000000000d+00
              qp(np) = qp_cell(2,iye_v,ize_v,is)
              y (np) = MESH_MAP(meshY,y_p_logical)
              z (np) = MESH_MAP(meshZ,z_p_logical)
            else
              y (np)=yb+myranf()*(ye-yb)
              z (np)=zb+myranf()*(ze-zb)
              qp(np)=particle_weight
            endif

!           Uniform mesh - Same as in version 5.0
!            ixe=hxi*x(np)+1.5000001d+00
!            iye=hyi*y(np)+0.5000001d+00
!            ize=hzi*z(np)+0.5000001d+00

!           Nonuniform mesh - using MESH_UNMAP
            rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
            rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
            rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
            ixe=rxe
            iye=rye
            ize=rze
            iye=iye-1             ! integer index in y direction starts at 0
            ize=ize-1             ! integer index in z direction starts at 0

            if (uniform_loading_in_logical_grid) then
              qinjl = qinjl + qp(np)
            endif
 
            ipstore=link(np)
            link(np)=iphead(ixe,iye,ize,is)
            iphead(ixe,iye,ize,is)=np
            ninjl=ninjl-1
          ENDDO




          jv=-int((vxbar(is)/vpar(is))/d_v)
          DO WHILE (ninjr.ne.0)
            np=ipstore
            fluxran=myranf()
            iv=int((fluxran/d_ranf)+0.5)
            iv=max(1,min(1001,iv))
            vx_tmp=vpar(is)*u_array_injection(iv,jv)
            vx(np)=-vx_tmp
            vmag=sqrt(-log(1.-.999999*myranf()))
            th=twopi*myranf()
            vy(np)=vmag*vper(is)*cos(th)
            vz(np)=vmag*vper(is)*sin(th)
            x(np)=xmax+vx(np)*dt*it_inject*myranf()


            if (uniform_loading_in_logical_grid) then
              Y_P_LOGICAL  = YB_LOGICAL+(YE_LOGICAL-YB_LOGICAL)*myranf()
              Z_P_LOGICAL  = ZB_LOGICAL+(ZE_LOGICAL-ZB_LOGICAL)*myranf()
              iye_v=dtyi*y_p_logical+1.00000000000d+00
              ize_v=dtzi*z_p_logical+1.00000000000d+00
              qp(np) = qp_cell(nx-1,iye_v,ize_v,is)
              y (np) = MESH_MAP(meshY,y_p_logical)
              z (np) = MESH_MAP(meshZ,z_p_logical)
            else
              y (np)=yb+myranf()*(ye-yb)
              z (np)=zb+myranf()*(ze-zb)
              qp(np)=particle_weight
            endif

!           Uniform mesh - Same as in version 5.0
!            ixe=hxi*x(np)+1.50
!            iye=hyi*y(np)+0.50
!            ize=hzi*z(np)+0.50

!           Nonuniform mesh - using MESH_UNMAP
            rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
            rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
            rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
            ixe=rxe
            iye=rye
            ize=rze
            iye=iye-1             ! integer index in y direction starts at 0
            ize=ize-1             ! integer index in z direction starts at 0

            if (uniform_loading_in_logical_grid) then
              qinjr = qinjr + qp(np)
            endif
 
 
            ipstore=link(np)
            link(np)=iphead(ixe,iye,ize,is)
            iphead(ixe,iye,ize,is)=np
            ninjr=ninjr-1
          ENDDO
!
!
          if (uniform_loading_in_logical_grid) then
!            write(6,*) " dqleft,qinjl = ",dqleft,qinjl/((ye-yb)*(ze-zb))
            qleft(is)=qleft(is)-qinjl/((ye-yb)*(ze-zb))
            qrite(is)=qrite(is)-qinjr/((ye-yb)*(ze-zb))
          else
            qleft(is) = qleft(is) - qinjl
            qrite(is) = qrite(is) - qinjr
          endif
!
!=======================================================================
!
	  return
	  end subroutine
!
!***********************************************************************
!
! ADDED BY HXV  1/09/03
!
      subroutine pdf_injection
!
!=======================================================================
!
!     solve the particle equations of motion using leapfrog
!
!     This subroutine is called by HERCULES
!
!     This subroutine calls LISTMKR
!
!=======================================================================
!
      use parameter_mod
      implicit none
      integer*8:: i,j,k,irec
      double precision:: d_fluxran,fluxran,vth,u_injection                 &
                        ,argum,vtherm,residue,d_v,direction
      double precision:: arguma(-1000:1000)
      double precision::fluxrana(1001),rata(-1000:1000),ssa(-1000:1000),sqrtpi
      logical:: start_from_beginning
!
!***********************************************************************
!
      direction=1.
      d_fluxran=1./1001.
      vth=1.
      d_v=1./16.
      j=0
      argum = j*d_v
      vtherm=vth
      start_from_beginning=.true.
    
      sqrtpi =dsqrt(pi)
      do k=-1000,1000
        arguma(k) = k*d_v
      enddo
      do k=-1000,1000
        rata(k) = arguma(k)/vtherm
      enddo
      do k=-1000,1000
        ssa(k) = vtherm*dexp(-rata(k)*rata(k))+arguma(k)*sqrtpi* &
                       (1.0d0+erf(dble(rata(k)))) 
      enddo

      do i=1,1001
        fluxrana(i)=d_fluxran*(dble(i)-0.5d0)
      enddo

      do i=1,1001 
        call zeroflux(u_array_injection(i,j),arguma(0),vtherm,fluxrana(i),residue,direction,start_from_beginning,ssa(0))
      enddo
!
!
      direction=1.
      start_from_beginning=.false.
      do j=1,1000
        do i=1,1001 
          u_array_injection(i,j)=u_array_injection(i,j-1)
        enddo
        do i=1,1001
          call zeroflux(u_array_injection(i,j),arguma(j),vtherm,fluxrana(i),residue,direction,start_from_beginning,ssa(j))
        enddo
      enddo
!
!
      direction=1.
      start_from_beginning=.false.
      do j=1,1000
        do i=1,1001 
          u_array_injection(i,-j)=u_array_injection(i,-j+1)
        enddo
        do i=1,1001
          call zeroflux(u_array_injection(i,-j),arguma(-j),vtherm,fluxrana(i),residue,direction,start_from_beginning,ssa(-j))
!------ CHECK SSA INDEX -----
!---------------------------- 
        enddo
      enddo
!
!=======================================================================
!
      return
      end
!
!*********************************************************************
!
!
! ADDED BY HXV  1/09/03
!
      subroutine zeroflux(x,argum,vtherm,xran,f,direction,start_from_beginning &
                          ,ssa)
!
!=====================================================================
!
!     this subroutine calculates a value x such that flux(x)=0
!
!=====================================================================
!
      use flux_mod
      implicit none
      integer*8:: ibisect,inewton,i,iran
      double precision:: x,argum,vtherm,xran,f,x0,f0,x1,f1,xleft,xrite &
      ,fleft,frite,dfdx,direction,dv,ssa,ssal
      logical:: start_from_beginning
!
!=====================================================================
!
!     find a reasonable initial guess for the Newton-Raphson iteration
!     by bisection
!
      ssal=ssa
      if (start_from_beginning) then
        x0=max(0.0d+00,argum-10.*vtherm*direction)
        f0=flux(x0,argum,vtherm,xran,ssal)
        ibisect=0
   5    ibisect=ibisect+1
        x1=x0+0.01*vtherm*direction
        f1=flux(x1,argum,vtherm,xran,ssal)
        if (f0*f1.gt.0..and.x1.lt.(argum+10.*vtherm*direction)) then
        x0=x1
        f0=f1
        goto 5
        endif
        x=0.5*(x1+x0)
      endif
!
!=====================================================================
!
!     Newton-Raphson iteration
!
      inewton=0
 10   inewton=inewton+1
      f     =flux(x,argum,vtherm,xran,ssal)
      xleft=(1.-1.d-03)*x
      fleft=flux(xleft,argum,vtherm,xran,ssal)
      xrite=(1.+1.d-03)*x
      frite=flux(xrite,argum,vtherm,xran,ssal)
      dfdx=(frite-fleft)/(2.d-03*x)
      x=x-f/dfdx
      if ((inewton.le.100.and.abs(f).ge.1.d-07).or.inewton.le.5) goto 10
!
!=====================================================================
!
      end subroutine zeroflux
