!
!#######################################################################
!
      subroutine field_2d
!
!=======================================================================
!
      use parameter_mod
      use MESH2D
!
!=======================================================================
!
      call date_and_time(values=time_begin_array(:,21))


      call date_and_time(values=time_begin_array(:,9))
      call pressgrad_2d(1)
      call date_and_time(values=time_end_array(:,9))
      call accumulate_time_difference(time_begin_array(1,9),time_end_array(1,9),time_elapsed(9))
 
 
      call date_and_time(values=time_begin_array(:,10))
      call bcalc_2d
      call date_and_time(values=time_end_array(:,10))
      call accumulate_time_difference(time_begin_array(1,10),time_end_array(1,10),time_elapsed(10))
 
 
      call date_and_time(values=time_begin_array(:,9))
      call pressgrad_2d(0)
      call date_and_time(values=time_end_array(:,9))
      call accumulate_time_difference(time_begin_array(1,9),time_end_array(1,9),time_elapsed(9))
 
 
      call date_and_time(values=time_begin_array(:,11))
      call ecalc_2d( 0)
      call date_and_time(values=time_end_array(:,11))
      call accumulate_time_difference(time_begin_array(1,11),time_end_array(1,11),time_elapsed(11))
 
 
      call date_and_time(values=time_begin_array(:,12))
      call focalc_2d
      call date_and_time(values=time_end_array(:,12))
      call accumulate_time_difference(time_begin_array(1,12),time_end_array(1,12),time_elapsed(12))
 
 
      call date_and_time(values=time_end_array(:,21))
      call accumulate_time_difference(time_begin_array(1,21),time_end_array(1,21),time_elapsed(21))
 
 
      return
      end
!
!################################################################################
!
      subroutine pressgrad_2d(iflag)

      use parameter_mod
      use MESH2D

      do k=kb,ke
        do j = jb,je
          do i=2,nx1
            dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))&
                +(1.-iflag)*den(i,j,k)
            a=1/dena

!           Uniform mesh - Same as is in version 5.0
!            dxa=a/(2.*hx)
!            dya=a/(2.*hy)
!            dza=a/(2.*hz)

 
!           Nonuniform mesh
            dxa=a/((meshX%dxn(i  )+meshX%dxn(i+1)))
            dya=a/((meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
            dza=a/((meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0
 


            dpedx(i,j,k)=(  (pe(i+1,j-1,k  )+2.*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4.    &
                          - (pe(i-1,j-1,k  )+2.*pe(i-1,j,k  )+pe(i-1,j+1,k  ))/4.    &
                         )*dxa
            dpedy(i,j,k)=(  (pe(i-1,j+1,k  )+2.*pe(i,j+1,k  )+pe(i+1,j+1,k  ))/4.    &
                          - (pe(i-1,j-1,k  )+2.*pe(i,j-1,k  )+pe(i+1,j-1,k  ))/4.    &
                         )*dya
            dpedz(i,j,k)=0.0
          enddo
        enddo
      enddo
 
      return
      end
!
!################################################################################
!
      subroutine ecalc_2d( iflag)

      use parameter_mod
      use MESH2D
      implicit none
      integer*4:: iflag, k, i, j
      real(kind=CUSTOM_REAL) :: bx1,bx2,bx3,bx4,by1,by2,by3,by4,bz1,bz2,bz3,bz4
      real(kind=CUSTOM_REAL) :: vixa,viya,viza,curlbx_scalar,curlby_scalar,curlbz_scalar
      real(kind=CUSTOM_REAL) :: bxav,byav,bzav,xj,yj,zj,tenx,tenz,bxx,byy,bzz
      real(kind=CUSTOM_REAL) :: dena,a,dxa,dya,dza,dbxdy,dbxdz,dbydx,dbydz,dbzdx,dbzdy
      real(kind=CUSTOM_REAL) :: teny,btot,curr_tot,tjdotb,xmtp1m,dexdy,dexdz,deydx,deydz
      real(kind=CUSTOM_REAL) :: dezdx,dezdy
 
      do k=kb,ke
        do j = jb,je
          do i=2,nx1
            bx1=bx(i+1,j+1,k)  +bdipole_x(i+1,j+1,k)
            bx2=bx(i  ,j+1,k)  +bdipole_x(i  ,j+1,k)
            bx3=bx(i  ,j  ,k)  +bdipole_x(i  ,j  ,k)
            bx4=bx(i+1,j  ,k)  +bdipole_x(i+1,j  ,k)
            by1=by(i+1,j+1,k)  +bdipole_y(i+1,j+1,k)
            by2=by(i  ,j+1,k)  +bdipole_y(i  ,j+1,k)
            by3=by(i  ,j  ,k)  +bdipole_y(i  ,j  ,k)
            by4=by(i+1,j  ,k)  +bdipole_y(i+1,j  ,k)
            bz1=bz(i+1,j+1,k)  +bdipole_z(i+1,j+1,k)
            bz2=bz(i  ,j+1,k)  +bdipole_z(i  ,j+1,k)
            bz3=bz(i  ,j  ,k)  +bdipole_z(i  ,j  ,k)
            bz4=bz(i+1,j  ,k)  +bdipole_z(i+1,j  ,k)
            vixa=(1.-iflag)*(1.5*vix(i,j,k)-0.5*vixo(i,j,k))&
                +iflag*vix(i,j,k)
            viya=(1.-iflag)*(1.5*viy(i,j,k)-0.5*viyo(i,j,k))&
                +iflag*viy(i,j,k)
            viza=(1.-iflag)*(1.5*viz(i,j,k)-0.5*vizo(i,j,k))&
                +iflag*viz(i,j,k)
            dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))&
                +(1.-iflag)*den(i,j,k)
            a=1/dena


!           Uniform mesh - Same as is in version 5.0
!            dxa=a/(2.*hx)
!            dya=a/(2.*hy)
!            dza=a/(2.*hz)

 
!          Nonuniform mesh
            dxa=a/(2.*meshX%dxc(i))
            dya=a/(2.*meshY%dxc(j+1))                  ! integer index in y direction starts at 0
            dza=a/(2.*meshZ%dxc(k+1))                  ! integer index in z direction starts at 0
 

            dbxdy=+bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
                  -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
            dbxdz= 0.0
            dbydx=+by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                  -by(i  ,j  ,k  )-by(i  ,j+1,k  )
            dbydz= 0.0
            dbzdx=+bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                  -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
            dbzdy=+bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                  -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
            curlbx_scalar=dya*dbzdy-dza*dbydz
            curlby_scalar=dza*dbxdz-dxa*dbzdx
            curlbz_scalar=dxa*dbydx-dya*dbxdy
            bxav=0.25*(bx1+bx2+bx3+bx4)
            byav=0.25*(by1+by2+by3+by4)
            bzav=0.25*(bz1+bz2+bz3+bz4)
!
! 6/25/2006 New eta_par option: tensor eta
!
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar

            if (eta_par.eq.0) then
		    tenx=eta(i,j,k)*xj
		    teny=eta(i,j,k)*yj
		    tenz=eta(i,j,k)*zj
	     else
		    bxx = bxav
		    byy = byav
		    bzz = bzav
		    btot = sqrt(bxx**2 + byy**2 + bzz**2)
                    xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
                    xtmp2m = 1.E-12_CUSTOM_REAL	!LAURA
	            curr_tot = max(xtmp2m,xtmp1m)
		 
	            if (eta_par.eq.1) then
		       tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
                       tenx = tjdotb*bxx/btot
		       teny = tjdotb*byy/btot
		       tenz = tjdotb*bzz/btot
		    else if (eta_par.eq.2) then
                       xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
                       xtmp2m = 1.E-12_CUSTOM_REAL	!LAURA
		       curr_tot = max(xtmp2m,xmtp1m)
                       tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
		       tenx = min(resis,tenx)
		       tenx = tenx*xj
                       teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
		       teny = min(resis,teny)
		       teny = teny*yj
                       tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
		       tenz = min(resis,tenz)
		       tenz = tenz*zj
		    endif
	     endif

             ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                       -dpedx(i,j,k)+tenx/a     
             ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                       -dpedy(i,j,k)+teny/a     
             ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                       -dpedz(i,j,k)+tenz/a                     
 
 
          enddo
        enddo
      enddo
 
 
      ex=ex*dipole_sphere_ex
      ey=ey*dipole_sphere_ey
      ez=ez*dipole_sphere_ez
!
!******************************************
!
!     boundary conditions
!
      call date_and_time(values=time_begin_array(:,18))
      call XREALBCC_PACK_E_2D(EX,EY,EZ,1,NX,NY,NZ)
      call date_and_time(values=time_end_array(:,18))
      call accumulate_time_difference(time_begin_array(1,18),time_end_array(1,18),time_elapsed(18))

      do k=kb-1,ke+1
        do j = jb-1,je+1
          ex(nx2,j,k)=-ex(nx1,j,k)
          ey(nx2,j,k)=+ey(nx1,j,k)
          ez(nx2,j,k)=+ez(nx1,j,k)
        enddo
      enddo
!
!***********************************************************
!
!        calculate curl E
!
      do k=kb,ke+1
         do j = jb,je+1
          do i=2,nx2
            dexdy= ex(i  ,j  ,k)+ex(i-1,j  ,k)&
                  -ex(i  ,j-1,k)-ex(i-1,j-1,k)
            dexdz= 0.0
            deydx= ey(i  ,j  ,k)+ey(i  ,j-1,k)&
                  -ey(i-1,j  ,k)-ey(i-1,j-1,k)
            deydz= 0.0
            dezdx= ez(i  ,j  ,k)+ez(i  ,j-1,k)&
                  -ez(i-1,j  ,k)-ez(i-1,j-1,k)
            dezdy= ez(i  ,j  ,k)+ez(i-1,j  ,k)&
                  -ez(i  ,j-1,k)-ez(i-1,j-1,k)

!           Uniform mesh - Same as is in version 5.0
!            curlex(i,j,k)=dezdy/(2.*hy)-deydz/(2.*hz)
!            curley(i,j,k)=dexdz/(2.*hz)-dezdx/(2.*hx)
!            curlez(i,j,k)=deydx/(2.*hx)-dexdy/(2.*hy)

 
!          Nonuniform mesh
            curlex(i,j,k)=dezdy/(2.*meshY%dxn(j+1))-deydz/(2.*meshZ%dxn(k+1))        ! integer index in y and z directions start  at 0
            curley(i,j,k)=dexdz/(2.*meshZ%dxn(k+1))-dezdx/(2.*meshX%dxn(i  ))        ! integer index in z       direction  starts at 0
            curlez(i,j,k)=deydx/(2.*meshX%dxn(i  ))-dexdy/(2.*meshY%dxn(j+1))        ! integer index in y       direction  starts at 0
 

          enddo
        enddo
      enddo
 

      return
      end
!
!################################################################################
!
      subroutine bcalc_2d

      use parameter_mod
      use MESH2D
      real(kind=CUSTOM_REAL):: tempx1(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
                        ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1)&
                        ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)
!
!=======================================================================
!
       call date_and_time(values=time_begin_array(:,22))


       dts=dt/real(iterb)
       dts2=dts/2.
       dts6=dts/6.
!
!***********************************
!   subcycle into iterb interations
!***********************************

      do 10 ii=1,iterb
 
        bxs=bx
        bys=by
        bzs=bz
!
!*******************
!   R-K first part
!*******************
!
         call date_and_time(values=time_begin_array(:,16))
         call ecalc_2d( 1 )
         call date_and_time(values=time_end_array(:,16))
         call accumulate_time_difference(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))
 
!***********************
!   B = B(n)+dt*K1/2
!***********************

         bx=bxs-dts2*curlex
         by=bys-dts2*curley
         bz=bzs-dts2*curlez
 
!******************
!   temp1 = K1
!******************

         tempx1=curlex
         tempy1=curley
         tempz1=curlez


!***************
!   R-K part 2
!***************

         call date_and_time(values=time_begin_array(:,16))
         call ecalc_2d( 1 )
         call date_and_time(values=time_end_array(:,16))
         call accumulate_time_difference(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))
 

!*********************
!  B = B(n)+dt*K2/2
!*********************

         bx=bxs-dts2*curlex
         by=bys-dts2*curley
         bz=bzs-dts2*curlez
!
!********************
!  temp2 = K2
!********************

         tempx1=tempx1+2.*curlex
         tempy1=tempy1+2.*curley
         tempz1=tempz1+2.*curlez


!*****************
!  R-K  part 3
!*****************

         call date_and_time(values=time_begin_array(:,16))
         call ecalc_2d( 1 )
         call date_and_time(values=time_end_array(:,16))
         call accumulate_time_difference(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))


!*********************
!  B = B(n)+dt*K3
!*********************

         bx=bxs-dts*curlex
         by=bys-dts*curley
         bz=bzs-dts*curlez

 
!*********************
!  temp3 = K3
!*********************

               tempx1=tempx1+2.*curlex
               tempy1=tempy1+2.*curley
               tempz1=tempz1+2.*curlez


!***************
!   R-K  part 4
!***************
 
         call date_and_time(values=time_begin_array(:,16))
         call ecalc_2d( 1 )
         call date_and_time(values=time_end_array(:,16)) 
         call accumulate_time_difference(time_begin_array(1,16),time_end_array(1,16),time_elapsed(16))
 

!*************************************
!  B = B(n) + dt*(K1+2K2+2K3+K4)/6
!*************************************

         bx=bxs-dts6*(tempx1+curlex)
         by=bys-dts6*(tempy1+curley)
         bz=bzs-dts6*(tempz1+curlez)

!************************
! restore boundary values
!************************

         bx(1,:   ,:   )=bxs(1,:   ,:   )
         bx(:,jb-1,:   )=bxs(:,jb-1,:   )
         bx(:,:   ,kb-1)=bxs(:,:   ,kb-1)
         by(1,:   ,:   )=bys(1,:   ,:   )
         by(:,jb-1,:   )=bys(:,jb-1,:   )
         by(:,:   ,kb-1)=bys(:,:   ,kb-1)
         bz(1,:   ,:   )=bzs(1,:   ,:   )
         bz(:,jb-1,:   )=bzs(:,jb-1,:   )
         bz(:,:   ,kb-1)=bzs(:,:   ,kb-1)


!************************
!  end of iteration loop
!************************

  10     continue
 
         CALL XREALBCC_PACK_B_2D(BX,BY,BZ,1,NX,NY,NZ)
 

      do k=kb-1,ke+1
        do j = jb-1,je+1
          bx(1  ,j,k)=bx(2  ,j,k)
          by(1  ,j,k)=by(2  ,j,k)
          bz(1  ,j,k)=bz(2  ,j,k)
        enddo
      enddo
 

       call date_and_time(values=time_end_array(:,22))
       call accumulate_time_difference(time_begin_array(1,22),time_end_array(1,22),time_elapsed(22))
 
 
       return
       end
!
!################################################################################
!
      subroutine focalc_2d

      use parameter_mod
      use MESH2D
      real(kind=CUSTOM_REAL):: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb,curr_tot	!LAURA

      do k=kb,ke 
        do j = jb,je
          do i=2,nx1
            dbxdy= bx(i+1,j+1,k)+bx(i  ,j+1,k)&
                  -bx(i+1,j  ,k)-bx(i  ,j  ,k)
            dbxdz= 0.0
            dbydx= by(i+1,j+1,k)+by(i+1,j  ,k)&
                  -by(i  ,j+1,k)-by(i  ,j  ,k)
            dbydz= 0.0
            dbzdx= bz(i+1,j+1,k)+bz(i+1,j  ,k)&
                  -bz(i  ,j+1,k)-bz(i  ,j  ,k)
            dbzdy= bz(i+1,j+1,k)+bz(i  ,j+1,k)&
                  -bz(i+1,j  ,k)-bz(i  ,j  ,k)


!           Uniform mesh - Same as is in version 5.0
!            curlbx_scalar=dbzdy/(2.*hy)-dbydz/(2.*hz)
!            curlby_scalar=dbxdz/(2.*hz)-dbzdx/(2.*hx)
!            curlbz_scalar=dbydx/(2.*hx)-dbxdy/(2.*hy)

 
!          Nonuniform mesh
            curlbx_scalar=dbzdy/(2.*meshY%dxc(j+1))-dbydz/(2.*meshZ%dxc(k+1))
            curlby_scalar=dbxdz/(2.*meshZ%dxc(k+1))-dbzdx/(2.*meshX%dxc(i  ))
            curlbz_scalar=dbydx/(2.*meshX%dxc(i  ))-dbxdy/(2.*meshY%dxc(j+1))
 

!
! 6/25/2006 New eta_par option: tensor eta
!
            bx1=bx(i+1,j+1,k)  +bdipole_x(i+1,j+1,k)
            bx2=bx(i  ,j+1,k)  +bdipole_x(i  ,j+1,k)
            bx3=bx(i  ,j  ,k)  +bdipole_x(i  ,j  ,k)
            bx4=bx(i+1,j  ,k)  +bdipole_x(i+1,j  ,k)
            by1=by(i+1,j+1,k)  +bdipole_y(i+1,j+1,k)
            by2=by(i  ,j+1,k)  +bdipole_y(i  ,j+1,k)
            by3=by(i  ,j  ,k)  +bdipole_y(i  ,j  ,k)
            by4=by(i+1,j  ,k)  +bdipole_y(i+1,j  ,k)
            bz1=bz(i+1,j+1,k)  +bdipole_z(i+1,j+1,k)
            bz2=bz(i  ,j+1,k)  +bdipole_z(i  ,j+1,k)
            bz3=bz(i  ,j  ,k)  +bdipole_z(i  ,j  ,k)
            bz4=bz(i+1,j  ,k)  +bdipole_z(i+1,j  ,k)
            bxav=.25*(bx1+bx2+bx3+bx4)
            byav=.25*(by1+by2+by3+by4)
            bzav=.25*(bz1+bz2+bz3+bz4)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar
	    bxx = bxav
	    byy = byav
	    bzz = bzav
            xtmp1m = sqrt(bxx**2 + byy**2 + bzz**2)
	    btot = max(1.E-12_CUSTOM_REAL,xtmp1m)	!LAURA
            xtmp1m = sqrt(xj**2 + yj**2 + zj**2) 
	    curr_tot = max(1.E-12_CUSTOM_REAL,xtmp1m)	!LAURA
            if (eta_par.eq.0) then
		    tenx=eta(i,j,k)*xj
		    teny=eta(i,j,k)*yj
		    tenz=eta(i,j,k)*zj
	    else
		    bxx = bxav
		    byy = byav
		    bzz = bzav
                    xtmp1m = sqrt(bxx**2 + byy**2 + bzz**2)
		    btot = max(1.E-12_CUSTOM_REAL,xtmp1m)	!LAURA
                    xtmp1m = sqrt(xj**2 + yj**2 + zj**2)
	            curr_tot = max(1.E-12_CUSTOM_REAL,xtmp1m)	!LAURA
		 
		    if (eta_par.eq.1) then
		       tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
                       tenx = tjdotb*bxx/btot
		       teny = tjdotb*byy/btot
		       tenz = tjdotb*bzz/btot
		    else if (eta_par.eq.2) then
                       tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
		       tenx = min(resis,tenx)
		       tenx = tenx*xj
                       teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
		       teny = min(resis,teny)
		       teny = teny*yj
                       tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
		       tenz = min(resis,tenz)
		       tenz = tenz*zj
		    endif
 
                    eta_times_b_dot_j(i,j,k) = min(resis,eta(i,j,k)*abs((bxav*xj + byav*yj + bzav*zj))/(curr_tot*btot))
	    endif


            fox(i,j,k)=-tenx
            foy(i,j,k)=-teny
            foz(i,j,k)=-tenz
 
 
          enddo
        enddo
      enddo
!
!******************************************
!
!     boundary conditions
!
!  first update internal ghost cells so that all processors
!  have latest information. Care must be exercised so that
!  BCs are set wrt proecessor that corresponds to that location. 
!  To that end, keep z loops set to limits of kb and ke.
!
       call XREALBCC_PACK_E_2D(fox,foy,foz,1,NX,NY,NZ)


      do k=kb-1,ke+1 
        do j = jb-1,je+1
          fox(1  ,j,k)=+fox(2  ,j,k)
          foy(1  ,j,k)=+foy(2  ,j,k)
          foz(1  ,j,k)=+foz(2  ,j,k)
          fox(nx2,j,k)=+fox(nx1,j,k)
          foy(nx2,j,k)=+foy(nx1,j,k)
          foz(nx2,j,k)=+foz(nx1,j,k)
        enddo
      enddo
 
      do k=kb-1,ke+1 
        do i=1,nx2
          if (jb.eq.1) then
            fox(i,jb-1,k)=fox(i,jb,k)
            foy(i,jb-1,k)=foy(i,jb,k)
            foz(i,jb-1,k)=foz(i,jb,k)
          endif
          if (je.eq.ny) then
            fox(i,je+1,k)=fox(i,je,k)
            foy(i,je+1,k)=foy(i,je,k)
            foz(i,je+1,k)=foz(i,je,k)
          endif
        enddo
      enddo

       if (kb.eq.1) then
         do j = jb-1,je+1
           do i=1,nx2
             fox(i,j,kb-1)=fox(i,j,kb)
             foy(i,j,kb-1)=foy(i,j,kb)
             foz(i,j,kb-1)=foz(i,j,kb)
           enddo
         enddo
      endif

      if (ke.eq.nz) then
        do j = jb-1,je+1
          do i=1,nx2
            fox(i,j,ke+1)=fox(i,j,ke)
            foy(i,j,ke+1)=foy(i,j,ke)
            foz(i,j,ke+1)=foz(i,j,ke)
          enddo
        enddo
      endif
 
      return
      end
!
!#######################################################################
!
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
      real(kind=CUSTOM_REAL):: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,epsilon	!LAURA
      real(kind=CUSTOM_REAL):: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi	!LAURA
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: &	!LAURA
      bx_av,by_av,bz_av
      real(kind=CUSTOM_REAL):: v_limit,eps2,rx0,ry0,rrat,sqrr,outer_radius,q_p	!LAURA
      real(kind=CUSTOM_REAL):: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min,x_disp,y_disp,z_disp          &	!LAURA
                        ,y_disp_max_p,x_disp_max_p,y_disp_max,x_disp_max
      integer*8:: Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p
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

 
      epsilon=1.E-05_CUSTOM_REAL	!LAURA
      outer_radius=dipole_sphere_radius+min(hx,hy)
      eps2=1.E-25_CUSTOM_REAL	!LAURA


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
!        v_limit=1.E+10_CUSTOM_REAL	!LAURA
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
!                rxb=hxi*x(l)+1.9999999999999999E+00_CUSTOM_REAL	!LAURA
!                ryb=hyi*y(l)+0.9999999999999999E+00_CUSTOM_REAL	!LAURA
!                rzb=hzi*z(l)+0.9999999999999999E+00_CUSTOM_REAL	!LAURA
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
                rxb=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                ryb=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                rzb=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                ixb=rxb
                iyb=ryb
                izb=rzb

                IZB=1
                fxb=rxb-ixb
                fyb=ryb-iyb
                fzb=rzb-izb
                iyb=iyb-1             ! integer index in y direction starts at 0
                izb=izb-1             ! integer index in z direction starts at 0

                if (fxb > 0.5) then
                  fxb = fxb -0.5
                  ixb = ixb + 1
                else
                  fxb=fxb + 0.5
                endif
                if (fyb > 0.5) then
                  fyb = fyb -0.5
                  iyb = iyb + 1
                else
                  fyb=fyb + 0.5
                endif
                iyb=min(iyb,je)


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
!                rxe=hxi*x(l)+1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                rye=hyi*y(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                rze=hzi*z(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
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
                rxe=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                rye=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                rze=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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

        call MPI_ALLREDUCE(x_disp_max_p,x_disp_max,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,IERR)	!LAURA
        call MPI_ALLREDUCE(y_disp_max_p,y_disp_max,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,IERR)	!LAURA
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
 
              if (y(l) <  epsilon) then

                  y(l) = 2.*epsilon - y(l)
                  vy(l)= 2.*vybar(is)-vy(l)
                  ypart= y(l)

              endif
 

              if (y(l) >  ymax-epsilon) then

                  y(l) = 2.*(ymax-epsilon) - y(l)
                  vy(l)= 2.*vybar(is)-vy(l)
                  ypart= y(l)

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
!              ixe=hxi*x(np)   +1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!              iye=hyi*y(np)   +0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!              ize=hzi*z(np)   +0.5000000000000001E+00_CUSTOM_REAL	!LAURA

!             Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
              call MPI_SEND(pdata,7,CUSTOM_MPI_TYPE,&	!LAURA
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
              call MPI_RECV(pdata,7,CUSTOM_MPI_TYPE,&	!LAURA
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
                rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
!                rx=hxi*x(l)+1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                ry=hyi*y(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                rz=hzi*z(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                ix=rx
!                iy=ry
!                iz=rz
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

         
!              Nonuniform mesh - using MESH_UNMAP
                rx=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                ry=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                rz=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
!
!#######################################################################
!
      subroutine xreal_2d(a,nx1m,ny1m,nz1m)
 
 
      use parameter_mod
      integer*8 i,j,nx1m,ny1m,nz1m,k
      real(kind=CUSTOM_REAL) a(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
     ,tmp(nxmax,jb-1:je+1,kb-1:ke+1)
 
      a(2   ,:,:   )=a(2   ,:,:   )+a(1   ,:,:   )
      a(nx1m+1,:,:   )=a(nx1m+1,:,:   )+a(nx1m+2,:,:   )
      a(:   ,:,kb  )=a(:   ,:,kb  )+a(:   ,:,kb-1)+a(:,:,kb+1)
      call MPI_SENDRECV(a    (1    ,je+1,kb-1),1,stridery,nbrrite,0, &
                        tmp  (1    ,jb-1,kb-1),1,stridery,nbrleft,0, &
                        mpi_comm_world,status,ierr)
      if (jb.eq.1) tmp(:,jb-1,:)=a(:,jb-1,:)
      a(:,jb,:)=a(:,jb,:)+tmp  (:,jb-1,:)
      call MPI_SENDRECV(a    (1    ,jb-1,kb-1),1,stridery,nbrleft,1, &
                        tmp  (1    ,je+1,kb-1),1,stridery,nbrrite,1, &
                        mpi_comm_world,status,ierr)
      if (je.eq.ny1m) tmp(:,je+1,:)=a(:,je+1,:)
      a(:,je,:)=a(:,je,:)+tmp  (:,je+1,:)
 
 
      a(1   ,:,:   )=a(2   ,:,: )
      a(nx1m+2,:,:   )=a(nx1m+1,:,: )
      a(:   ,:,kb-1)=a(:   ,:,kb)
      a(:   ,:,kb+1)=a(:   ,:,kb)
 
 
      return
      end
!
!#######################################################################
!
      subroutine xrealbcc_2d(a, ibnd, nx1m, ny1m,nz1m)
 
 
      use parameter_mod
      integer*8      i,j,nx1m,ny1m,nz1m
      integer   ibnd
      real(kind=CUSTOM_REAL) a(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
     ,tmp(nxmax,jb-1:je+1,kb-1:ke+1)
 
      tmp=a
 
 
      a(1   ,:,:   )=a(2   ,:,: )
      a(nx1m+2,:,:   )=a(nx1m+1,:,: )
      a(:   ,:,kb-1)=a(:   ,:,kb)
      a(:   ,:,kb+1)=a(:   ,:,kb)
 
 
      call MPI_SENDRECV(a(1    ,je  ,kb-1),1,stridery,nbrrite,0,&
                        a(1    ,jb-1,kb-1),1,stridery,nbrleft,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(a(1    ,jb  ,kb-1),1,stridery,nbrleft,1,&
                        a(1    ,je+1,kb-1),1,stridery,nbrrite,1,&
                        mpi_comm_world,status,ierr)
      if (jb.eq.1) then
        if(ibnd.eq.1)  then
          a(:,jb-1,:)=a(:,jb,:)
        else
          a(:,jb-1,:)=tmp(:,jb-1,:)
        endif
      endif
      if (je.eq.ny1m) then
        if(ibnd.eq.1)  then
          a(:,je+1,:)=a(:,je,:)
        else
          a(:,je+1,:)=tmp(:,je+1,:)
        endif
      endif
 
 
      return
      end
!
!#######################################################################
!
      subroutine xrealbcc_pack_b_2d(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)
 
 
      use parameter_mod
      integer*8      i,j,nx1m,ny1m,nz1m,k
      integer   ibnd
      real(kind=CUSTOM_REAL) a_x(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
                      ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                      ,packed_data_xy_recv(nxmax,jb-1:je+1,3)
 
          a_x(:,:,kb-1)=a_x(:,:,kb)
          a_y(:,:,kb-1)=a_y(:,:,kb)
          a_z(:,:,kb-1)=a_z(:,:,kb)
          a_x(:,:,ke+1)=a_x(:,:,ke)
          a_y(:,:,ke+1)=a_y(:,:,ke)
          a_z(:,:,ke+1)=a_z(:,:,ke)
 
 
      do k=kb-1,ke+1
        do i=1,nxmax
          packed_data_xz_send(i,k,1)=a_x(i,je,k)
          packed_data_xz_send(i,k,2)=a_y(i,je,k)
          packed_data_xz_send(i,k,3)=a_z(i,je,k)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
                        mpi_comm_world,status,ierr)
      if (jb.eq.1) then
        if(ibnd.eq.1)  then
          a_x(:,jb-1,:)=a_x(:,jb,:)
          a_y(:,jb-1,:)=a_y(:,jb,:)
          a_z(:,jb-1,:)=a_z(:,jb,:)
        endif
      else
        do k=kb-1,ke+1
          do i=1,nxmax
            a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
            a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
            a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
          enddo
        enddo
      endif
 
 
      do k=kb-1,ke+1
        do i=1,nxmax
          packed_data_xz_send(i,k,1)=a_x(i,jb,k)
          packed_data_xz_send(i,k,2)=a_y(i,jb,k)
          packed_data_xz_send(i,k,3)=a_z(i,jb,k)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
                        mpi_comm_world,status,ierr)
      if (jb.eq.ny1m) then
        if(ibnd.eq.1)  then
          a_x(:,je+1,:)=a_x(:,je,:)
          a_y(:,je+1,:)=a_y(:,je,:)
          a_z(:,je+1,:)=a_z(:,je,:)
        endif
      else
        do k=kb-1,ke+1
          do i=1,nxmax
            a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
            a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
            a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
          enddo
        enddo
      endif
 
 
      return
      end
!
!***********************************************************************
!
      subroutine caltemp2_global_2d
 
 
      use parameter_mod
      use MESH2D
      real(kind=CUSTOM_REAL):: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi	!LAURA
      integer*8 ix,iy,iz,ixp1,iyp1,izp1

 
 
      call date_and_time(values=time_begin_array(:,23))
 
 
      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt
 


      tpar  = 0.
      tperp = 0.
      rfrac = 0.
 
      if (nspec >= 2) rfrac = frac(2)/frac(1)
 
 
      call date_and_time(values=time_begin_array(:,26))
      DO IS=1,NSPEC
        wmult=wspec(is)
        h=dt*qspec(is)/wmult
        hh=.5*h
        dpedx = 0.

        DO IIZE = KB-1,KE
         DO IIYE = JB-1,JE
          DO IIXE = 1, NX1
            NP=IPHEAD(IIXE,IIYE,IIZE,IS)
            DO WHILE (NP.NE.0)
              L=NP

!             Uniform mesh - Same as is in version 5.0
!              rx=hxi*x(l)+1.5000000000000001
!              ry=hyi*y(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!              rz=hzi*z(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!              ix=rx
!              iy=ry
!              iz=rz
!              IZ=1
!              fx=rx-ix
!              fy=ry-iy
!              fz=rz-iz

 
!             Nonuniform mesh - using MESH_UNMAP
              rx=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              ry=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              rz=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              ix=rx
              iy=ry
              iz=rz
              IZ=1
              fx=rx-ix
              fy=ry-iy
              fz=rz-iz
              iy=iy-1             ! integer index in y direction starts at 0
              iz=iz-1             ! integer index in z direction starts at 0
 


              w1=(1.-fx)*(1.-fy)
              w2=    fx *(1.-fy)
              w3=(1.-fx)*    fy
              w4=    fx*     fy

              dns1= dns(ix,iy  ,iz  ,1)*w1+dns(ix+1,iy  ,iz  ,1)*w2&
              +     dns(ix,iy+1,iz  ,1)*w3+dns(ix+1,iy+1,iz  ,1)*w4

              dns2= 0.

              dnst = dns1 + dns2
      
              vxavg1=vxs(ix,iy  ,iz  ,1)*w1+vxs(ix+1,iy  ,iz  ,1)*w2&
              +      vxs(ix,iy+1,iz  ,1)*w3+vxs(ix+1,iy+1,iz  ,1)*w4

              vxavg2= 0.

              vxavg = (dns1*vxavg1 + dns2*vxavg2)/dnst

              vyavg1=vys(ix,iy  ,iz  ,1)*w1+vys(ix+1,iy  ,iz  ,1)*w2&
              +      vys(ix,iy+1,iz  ,1)*w3+vys(ix+1,iy+1,iz  ,1)*w4

              vyavg2=0.

              vyavg = (dns1*vyavg1 + dns2*vyavg2)/dnst

              vzavg1=vzs(ix,iy  ,iz  ,1)*w1+vzs(ix+1,iy  ,iz  ,1)*w2&
              +      vzs(ix,iy+1,iz  ,1)*w3+vzs(ix+1,iy+1,iz  ,1)*w4

              vzavg2=0.

              vzavg = (dns1*vzavg1 + dns2*vzavg2)/dnst

              vxa=vx(l)-vxavg
              vya=vy(l)-vyavg
              vza=vz(l)-vzavg

              bxa  =bx       (ix,iy  ,iz  )*w1+bx       (ix+1,iy  ,iz  )*w2&
              +     bx       (ix,iy+1,iz  )*w3+bx(       ix+1,iy+1,iz  )*w4&
              +     bdipole_x(ix,iy  ,iz  )*w1+bdipole_x(ix+1,iy  ,iz  )*w2&
              +     bdipole_x(ix,iy+1,iz  )*w3+bdipole_x(ix+1,iy+1,iz  )*w4

              bya  =by(       ix,iy  ,iz  )*w1+by(       ix+1,iy  ,iz  )*w2&
              +     by(       ix,iy+1,iz  )*w3+by       (ix+1,iy+1,iz  )*w4&
              +     bdipole_y(ix,iy  ,iz  )*w1+bdipole_y(ix+1,iy  ,iz  )*w2&
              +     bdipole_y(ix,iy+1,iz  )*w3+bdipole_y(ix+1,iy+1,iz  )*w4

              bza  =bz       (ix,iy  ,iz  )*w1+bz       (ix+1,iy  ,iz  )*w2&
              +     bz       (ix,iy+1,iz  )*w3+bz       (ix+1,iy+1,iz  )*w4&
              +     bdipole_z(ix,iy  ,iz  )*w1+bdipole_z(ix+1,iy  ,iz  )*w2&
              +     bdipole_z(ix,iy+1,iz  )*w3+bdipole_z(ix+1,iy+1,iz  )*w4

              btota=sqrt(bxa**2+bya**2+bza**2)
              if(btota.lt.1.e-20) btota=1.e-20
              wpar=(vxa*bxa+vya*bya+vza*bza)/btota
              wperp2=vxa**2+vya**2+vza**2-wpar**2

              tpar (ix  ,iy  ,iz,is)=tpar (ix  ,iy  ,iz,is)+qp(np)*w1*wpar*wpar
              tpar (ix+1,iy  ,iz,is)=tpar (ix+1,iy  ,iz,is)+qp(np)*w2*wpar*wpar 
              tpar (ix  ,iy+1,iz,is)=tpar (ix  ,iy+1,iz,is)+qp(np)*w3*wpar*wpar 
              tpar (ix+1,iy+1,iz,is)=tpar (ix+1,iy+1,iz,is)+qp(np)*w4*wpar*wpar 
              tperp(ix  ,iy  ,iz,is)=tperp(ix  ,iy  ,iz,is)+qp(np)*w1*wperp2 
              tperp(ix+1,iy  ,iz,is)=tperp(ix+1,iy  ,iz,is)+qp(np)*w2*wperp2 
              tperp(ix  ,iy+1,iz,is)=tperp(ix  ,iy+1,iz,is)+qp(np)*w3*wperp2 
              tperp(ix+1,iy+1,iz,is)=tperp(ix+1,iy+1,iz,is)+qp(np)*w4*wperp2
              dpedx(ix  ,iy  ,iz)=dpedx(ix  ,iy  ,iz)+qp(np)*w1
              dpedx(ix+1,iy  ,iz)=dpedx(ix+1,iy  ,iz)+qp(np)*w2 
              dpedx(ix  ,iy+1,iz)=dpedx(ix  ,iy+1,iz)+qp(np)*w3 
              dpedx(ix+1,iy+1,iz)=dpedx(ix+1,iy+1,iz)+qp(np)*w4
 
 
              np=link(np)
            ENDDO
          ENDDO
         ENDDO
        ENDDO

        DO IZ = KB-1,KE
          DO IY = JB-1,JE
            DO IX = 1, NX1
              if (dpedx(ix,iy,iz) /= 0.) then
                tpar (ix,iy,iz,is) = tpar (ix,iy,iz,is)/(   tx0(is)*dpedx(ix,iy,iz))
                tperp(ix,iy,iz,is) = tperp(ix,iy,iz,is)/(2.*tx0(is)*dpedx(ix,iy,iz))
              endif
            ENDDO
          ENDDO
        ENDDO

      ENDDO
 
      call date_and_time(values=time_end_array(:,26))
      call accumulate_time_difference(time_begin_array(1,26),time_end_array(1,26),time_elapsed(26))
!
!

      do is=1,nspec
        call date_and_time(values=time_begin_array(:,24))
        call XREAL_2D(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL_2D(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
        call date_and_time(values=time_end_array(:,24))
        call accumulate_time_difference(time_begin_array(1,24),time_end_array(1,24),time_elapsed(24))
 
 
        call date_and_time(values=time_begin_array(:,25))
        call XREALBCC_2D(tpar (1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC_2D(tperp(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call date_and_time(values=time_end_array(:,25))
        call accumulate_time_difference(time_begin_array(1,25),time_end_array(1,25),time_elapsed(25))
      enddo


      call date_and_time(values=time_begin_array(:,26))
 
!      GOTO 10  
!
!      do is=1,nspec
!        do k=kb-1,ke+1
!         do j = jb-1,je+1
!            do i=1,nx2
!              if(is.eq.1) then
!                dns1=dns(i,j,k,1)/(dfac(1)*frac(1))
!                dns2=0.
!                denum=dns1+rfrac*dns2
!              else
!                denum=dns(i,j,k,is)/(dfac(is)*frac(is))
!              endif
!              if(den(i,j,k).le.E_CUSTOM_REALenmin)  then	!LAURA
!                tpar(i,j,k,is)=1.e-5
!                tperp(i,j,k,is)=1.e-5
!              else
!                denum=denum*tx0(is)
!                tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
!                tperp(i,j,k,is)=0.5*tperp(i,j,k,is)*wspec(is)/denum
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
!
! 10   CONTINUE

      call date_and_time(values=time_end_array(:,26))
      call accumulate_time_difference(time_begin_array(1,26),time_end_array(1,26),time_elapsed(26))

 
      call date_and_time(values=time_end_array(:,23))
      call accumulate_time_difference(time_begin_array(1,23),time_end_array(1,23),time_elapsed(23))
 

      return
      end
!
!#######################################################################
!
      subroutine xrealbcc_pack_e_2d(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)
 
 
      use parameter_mod
      integer*8      i,j,nx1m,ny1m,nz1m,k
      integer   ibnd
      real(kind=CUSTOM_REAL) a_x(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
                      ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                      ,packed_data_xy_recv(nxmax,jb-1:je+1,3)
 
          a_x(:,:,kb-1)=a_x(:,:,kb)
          a_y(:,:,kb-1)=a_y(:,:,kb)
          a_z(:,:,kb-1)=a_z(:,:,kb)
          a_x(:,:,ke+1)=a_x(:,:,ke)
          a_y(:,:,ke+1)=a_y(:,:,ke)
          a_z(:,:,ke+1)=a_z(:,:,ke)
 
 
      do k=kb-1,ke+1
        do i=1,nxmax
          packed_data_xz_send(i,k,1)=a_x(i,je,k)
          packed_data_xz_send(i,k,2)=a_y(i,je,k)
          packed_data_xz_send(i,k,3)=a_z(i,je,k)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
                        mpi_comm_world,status,ierr)
      if (jb.eq.1) then
        if(ibnd.eq.1)  then
          a_x(:,jb-1,:)= a_x(:,jb,:)
          a_y(:,jb-1,:)=-a_y(:,jb,:)
          a_z(:,jb-1,:)= a_z(:,jb,:)
        endif
      else
        do k=kb-1,ke+1
          do i=1,nxmax
            a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
            a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
            a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
          enddo
        enddo
      endif
 
 
      do k=kb-1,ke+1
        do i=1,nxmax
          packed_data_xz_send(i,k,1)=a_x(i,jb,k)
          packed_data_xz_send(i,k,2)=a_y(i,jb,k)
          packed_data_xz_send(i,k,3)=a_z(i,jb,k)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
                        mpi_comm_world,status,ierr)
      if (jb.eq.ny1m) then
        if(ibnd.eq.1)  then
          a_x(:,je+1,:)= a_x(:,je,:)
          a_y(:,je+1,:)=-a_y(:,je,:)
          a_z(:,je+1,:)= a_z(:,je,:)
        endif
      else
        do k=kb-1,ke+1
          do i=1,nxmax
            a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
            a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
            a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
          enddo
        enddo
      endif
 
 
      return
      end
!
!#######################################################################
!
      subroutine nsmth_2d (a,nx2m,ny2m,nz2m)
 

      use parameter_mod
      integer*8 :: nx2m,ny2m,nz2m
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: temp, a	!LAURA
 
!************************
!  smoothing routine--assumes aperiodic in x
!************************
      call XREALBCC_2D(a,0,NX,NY,NZ)
      temp=a

      do k=kb-1,ke+1
        do j = jb,je
          do i=2,nx1
            a(i,j,k)=temp(i,j,k)/4.&
                     +( temp(i-1,j  ,k)+temp(i+1,j  ,k)+temp(i  ,j+1,k)   &
                     +temp(i  ,j-1,k))/8.&
                     +( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k)   & 
                     +temp(i-1,j-1,k))/16.
          enddo
        enddo
      enddo
 
      if (jb.eq.1) then
        do k=kb,ke
          do i=2,nx1
            a(i,jb-1,k)=a(i,jb,k)
          enddo
        enddo
      endif

      if (je.eq.ny) then
        do k=kb,ke
          do i=2,nx1
            a(i,je+1,k)=a(i,je,k)
          enddo
        enddo
      endif
 
 
      do k=kb-1,ke+1
        do j = jb-1,je+1
          a(1  ,j,k)=a(2  ,j,k)
          a(nx2m,j,k)=a(nx1,j,k)
        enddo
      enddo

      call XREALBCC_2D(a,0,NX,NY,NZ)

      return
      end

