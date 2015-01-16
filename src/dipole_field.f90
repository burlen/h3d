!
!***********************************************************************
!
      subroutine dipole_field
 
 
      use parameter_mod
      use functions_f90
      use erf_mod
      use MESH2D
 
! Mirror dipole
 
      real(kind=CUSTOM_REAL) :: x_dipole_mirror(2),y_dipole_mirror(2),z_dipole_mirror(2)
      real(kind=CUSTOM_REAL) :: x_dipole, y_dipole , z_dipole ! position of the dipole
      real(kind=CUSTOM_REAL) :: r_dipole , th_dipole , phi 
      real(kind=CUSTOM_REAL) :: b_r , b_t, b_eta
      real(kind=CUSTOM_REAL) :: x_pos , y_pos , z_pos
      INTEGER*8 :: i , j , k, i_array
      
 
!        Uniform mesh - Same as is in version 5.0
!         x_dipole=hx*(i_dipole+0.5)
!         y_dipole=hy*(j_dipole+0.5)
!         z_dipole=hz*(k_dipole+0.5)
 
!        Nonuniform mesh
         x_dipole=meshX%xc(i_dipole)
         y_dipole=meshY%xc(j_dipole)
         z_dipole=meshZ%xc(k_dipole)

         x_dipole_mirror(1) =        x_dipole
         y_dipole_mirror(1) =        y_dipole
         z_dipole_mirror(1) =        z_dipole
         x_dipole_mirror(2) =       -x_dipole
         y_dipole_mirror(2) =        y_dipole
         z_dipole_mirror(2) =        z_dipole
 
 
! Set up bdipole_x

	  do k = kb-1,ke+1

!           Uniform mesh - Same as is in version 5.0
!	    z_pos = real(k-1)*hz+0.5*hz	!LAURA
 
!           Nonuniform mesh
	    z_pos = meshZ%xc(k+1)         ! integer index in z direction starts  at 0

		do j = jb-1,je+1

!                 Uniform mesh - Same as is in version 5.0
!		  y_pos = real(j-1)*hy+0.5*hy	!LAURA

!                 Nonuniform mesh
		  y_pos = meshY%xc(j+1)   ! integer index in y direction starts  at 0
 
		  do i = 1,nx2

!                   Uniform mesh - Same as is in version 5.0
!		    x_pos = real(i-2)*hx+0.5*hx	!LAURA

!                   Nonuniform mesh
		    x_pos = meshX%xc(i)

                    bdipole_x(i,j,k)=0._CUSTOM_REAL
                    bdipole_y(i,j,k)=0._CUSTOM_REAL
                    bdipole_z(i,j,k)=0._CUSTOM_REAL
                        do i_array=1,2
		          r_dipole = dsqrt( (x_pos-x_dipole_mirror(i_array))**2_CUSTOM_REAL &
                                         + (y_pos-y_dipole_mirror(i_array))**2_CUSTOM_REAL &
                                         + (z_pos-z_dipole_mirror(i_array))**2_CUSTOM_REAL )
		   	  if ((r_dipole /= 0.)) then
		            th_dipole = dacos( (y_pos-y_dipole_mirror(i_array))/r_dipole )
			    if(x_dipole_mirror(i_array) /= x_pos) then
		              phi = datan (dabs(z_pos-z_dipole_mirror(i_array)) &
                                         /dabs(x_pos-x_dipole_mirror(i_array)))
			      if ((phi==0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi
   			      if (((z_pos-z_dipole_mirror(i_array))>0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi-phi
			      if (((z_pos-z_dipole_mirror(i_array))<0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi+phi
			      if (((z_pos-z_dipole_mirror(i_array))<0.).and.((x_pos-x_dipole_mirror(i_array))>0.)) phi=2._CUSTOM_REAL*pi-phi
			    else
			      if((z_dipole_mirror(i_array)-z_pos) > 0.) phi = pi / 2._CUSTOM_REAL
			      if((z_dipole_mirror(i_array)-z_pos) < 0.) phi = 3._CUSTOM_REAL*pi/2._CUSTOM_REAL
			    endif
		            b_r = -2._CUSTOM_REAL*dipole_moment*dcos(th_dipole)/r_dipole**3_CUSTOM_REAL
		            b_t = - dipole_moment*dsin(th_dipole)/r_dipole**3_CUSTOM_REAL
		            b_eta = b_r*dsin(th_dipole)+b_t*dcos(th_dipole)
                              bdipole_x(i,j,k) = bdipole_x(i,j,k)+b_eta*dcos(phi)
                              bdipole_z(i,j,k) = bdipole_z(i,j,k)+b_eta*dsin(phi)
                              bdipole_y(i,j,k) = bdipole_y(i,j,k)+b_r*dcos(th_dipole)-b_t*dsin(th_dipole)
			  else
			    if (myid == 0) write(6,*) 'r_dipole == 0 '
			    bdipole_y(i,j,k) = bdipole_y(i,j,k)-2._CUSTOM_REAL*dipole_moment/hz**3_CUSTOM_REAL
			  endif
                        enddo
		  enddo
		enddo
          enddo
          call xrealbcc(bdipole_x,1,nx,ny,nz)
          call xrealbcc(bdipole_y,1,nx,ny,nz)
          call xrealbcc(bdipole_z,1,nx,ny,nz)
 
 
          bdipole_z=bdipole_z+bz_IMF/wpiwci
 
 
      dipole_sphere_ex=1._CUSTOM_REAL
      dipole_sphere_ey=1._CUSTOM_REAL
      dipole_sphere_ez=1._CUSTOM_REAL

      do k = kb-1,ke+1

!       Uniform mesh - Same as is in version 5.0
!       z_pos = real(k-1)*hz+0.5*hz	!LAURA
 
!       Nonuniform mesh
        z_pos = meshZ%xc(k+1)         ! integer index in z direction starts  at 0

	do j = jb-1,je+1

!         Uniform mesh - Same as is in version 5.0
!	  y_pos = real(j-1)*hy+0.5*hy	!LAURA

!         Nonuniform mesh
	  y_pos = meshY%xc(j+1)   ! integer index in y direction starts  at 0
 
	  do i = 1,nx2

!           Uniform mesh - Same as is in version 5.0
!	    x_pos = real(i-2)*hx+0.5*hx	!LAURA

!           Nonuniform mesh
	    x_pos = meshX%xc(i)

            r_dipole = dsqrt((x_pos-x_dipole)**2_CUSTOM_REAL+(y_pos-y_dipole)**2_CUSTOM_REAL+(z_pos-z_dipole)**2_CUSTOM_REAL)
            if (r_dipole < dipole_sphere_radius) then
              dipole_sphere_ex(i,j,k)=0._CUSTOM_REAL
              dipole_sphere_ey(i,j,k)=0._CUSTOM_REAL
              dipole_sphere_ez(i,j,k)=0._CUSTOM_REAL
            endif
          enddo
        enddo
      enddo
      call xrealbcc(dipole_sphere_ex,1,nx,ny,nz)
      call xrealbcc(dipole_sphere_ey,1,nx,ny,nz)
      call xrealbcc(dipole_sphere_ez,1,nx,ny,nz)
 
      return
      end
