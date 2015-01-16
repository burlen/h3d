!
!***********************************************************************
!
      subroutine dipole_field
 
 
      use parameter_mod
      use functions_f90
      use erf_mod
      use MESH2D
 
! Mirror dipole
 
      REAL :: x_dipole_mirror(2),y_dipole_mirror(2),z_dipole_mirror(2)
      REAL :: x_dipole, y_dipole , z_dipole ! position of the dipole
	  REAL :: r_dipole , th_dipole , phi 
	  REAL :: b_r , b_th, b_eta
	  REAL :: x_pos , y_pos , z_pos
	  INTEGER*8 :: i , j , k
      
 
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
!	    z_pos = float(k-1)*hz+0.5*hz
 
!           Nonuniform mesh
	    z_pos = meshZ%xc(k+1)         ! integer index in z direction starts  at 0

		do j = jb-1,je+1

!                 Uniform mesh - Same as is in version 5.0
!		  y_pos = float(j-1)*hy+0.5*hy

!                 Nonuniform mesh
		  y_pos = meshY%xc(j+1)   ! integer index in y direction starts  at 0
 
		  do i = 1,nx2

!                   Uniform mesh - Same as is in version 5.0
!		    x_pos = float(i-2)*hx+0.5*hx

!                   Nonuniform mesh
		    x_pos = meshX%xc(i)

                    bdipole_x(i,j,k)=0.
                    bdipole_y(i,j,k)=0.
                    bdipole_z(i,j,k)=0.
                      if (image_dipole) then
                        do i_array=1,2
		          r_dipole = sqrt( (x_pos-x_dipole_mirror(i_array))**2 &
                                         + (y_pos-y_dipole_mirror(i_array))**2 &
                                         + (z_pos-z_dipole_mirror(i_array))**2 )
		   	  if ((r_dipole /= 0.)) then
		            th_dipole = acos( (y_pos-y_dipole_mirror(i_array))/r_dipole )
			    if(x_dipole_mirror(i_array) /= x_pos) then
		              phi = atan (abs(z_pos-z_dipole_mirror(i_array)) &
                                         /abs(x_pos-x_dipole_mirror(i_array)))
			      if ((phi==0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi
   			      if (((z_pos-z_dipole_mirror(i_array))>0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi-phi
			      if (((z_pos-z_dipole_mirror(i_array))<0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi+phi
			      if (((z_pos-z_dipole_mirror(i_array))<0.).and.((x_pos-x_dipole_mirror(i_array))>0.)) phi=2.*pi-phi
			    else
			      if((z_dipole_mirror(i_array)-z_pos) > 0.) phi = pi / 2.
			      if((z_dipole_mirror(i_array)-z_pos) < 0.) phi = 3*pi/2.
			    endif
		            b_r = -2.*dipole_moment*cos(th_dipole)/r_dipole**3
		            b_t = - dipole_moment*sin(th_dipole)/r_dipole**3
		            b_eta = b_r*sin(th_dipole)+b_t*cos(th_dipole)
                              bdipole_x(i,j,k) = bdipole_x(i,j,k)+b_eta*cos(phi)
                              bdipole_z(i,j,k) = bdipole_z(i,j,k)+b_eta*sin(phi)
                              bdipole_y(i,j,k) = bdipole_y(i,j,k)+b_r*cos(th_dipole)-b_t*sin(th_dipole)
			  else
			    if (myid == 0) write(6,*) 'r_dipole == 0 '
			    bdipole_y(i,j,k) = bdipole_y(i,j,k)-2.*dipole_moment/hz**3
			  endif
                        enddo
                      else
                          i_array = 1
		          r_dipole = sqrt( (x_pos-x_dipole_mirror(i_array))**2 &
                                         + (y_pos-y_dipole_mirror(i_array))**2 &
                                         + (z_pos-z_dipole_mirror(i_array))**2 )
		   	  if ((r_dipole /= 0.)) then
		            th_dipole = acos( (y_pos-y_dipole_mirror(i_array))/r_dipole )
			    if(x_dipole_mirror(i_array) /= x_pos) then
		              phi = atan (abs(z_pos-z_dipole_mirror(i_array)) &
                                         /abs(x_pos-x_dipole_mirror(i_array)))
			      if ((phi==0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi
   			      if (((z_pos-z_dipole_mirror(i_array))>0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi-phi
			      if (((z_pos-z_dipole_mirror(i_array))<0.).and.((x_pos-x_dipole_mirror(i_array))<0.)) phi=pi+phi
			      if (((z_pos-z_dipole_mirror(i_array))<0.).and.((x_pos-x_dipole_mirror(i_array))>0.)) phi=2.*pi-phi
			    else
			      if((z_dipole_mirror(i_array)-z_pos) > 0.) phi = pi / 2.
			      if((z_dipole_mirror(i_array)-z_pos) < 0.) phi = 3*pi/2.
			    endif
		            b_r = -2.*dipole_moment*cos(th_dipole)/r_dipole**3
		            b_t = - dipole_moment*sin(th_dipole)/r_dipole**3
		            b_eta = b_r*sin(th_dipole)+b_t*cos(th_dipole)
                              bdipole_x(i,j,k) = bdipole_x(i,j,k)+b_eta*cos(phi)
                              bdipole_z(i,j,k) = bdipole_z(i,j,k)+b_eta*sin(phi)
                              bdipole_y(i,j,k) = bdipole_y(i,j,k)+b_r*cos(th_dipole)-b_t*sin(th_dipole)
			  else
			    if (myid == 0) write(6,*) 'r_dipole == 0 '
			    bdipole_y(i,j,k) = bdipole_y(i,j,k)-2.*dipole_moment/hz**3
			  endif
                      endif
		  enddo
		enddo
          enddo
          call xrealbcc(bdipole_x,1,nx,ny,nz)
          call xrealbcc(bdipole_y,1,nx,ny,nz)
          call xrealbcc(bdipole_z,1,nx,ny,nz)
 
 
          bdipole_z=bdipole_z+bz_IMF/wpiwci
 
 
      dipole_sphere_ex=1.
      dipole_sphere_ey=1.
      dipole_sphere_ez=1.

      do k = kb-1,ke+1

!       Uniform mesh - Same as is in version 5.0
!       z_pos = float(k-1)*hz+0.5*hz
 
!       Nonuniform mesh
        z_pos = meshZ%xc(k+1)         ! integer index in z direction starts  at 0

	do j = jb-1,je+1

!         Uniform mesh - Same as is in version 5.0
!	  y_pos = float(j-1)*hy+0.5*hy

!         Nonuniform mesh
	  y_pos = meshY%xc(j+1)   ! integer index in y direction starts  at 0
 
	  do i = 1,nx2

!           Uniform mesh - Same as is in version 5.0
!	    x_pos = float(i-2)*hx+0.5*hx

!           Nonuniform mesh
	    x_pos = meshX%xc(i)

            r_dipole = sqrt((x_pos-x_dipole)**2+(y_pos-y_dipole)**2+(z_pos-z_dipole)**2)
            if (r_dipole < dipole_sphere_radius) then
              dipole_sphere_ex(i,j,k)=0.
              dipole_sphere_ey(i,j,k)=0.
              dipole_sphere_ez(i,j,k)=0.
            endif
          enddo
        enddo
      enddo
      call xrealbcc(dipole_sphere_ex,1,nx,ny,nz)
      call xrealbcc(dipole_sphere_ey,1,nx,ny,nz)
      call xrealbcc(dipole_sphere_ez,1,nx,ny,nz)
 
      return
      end
