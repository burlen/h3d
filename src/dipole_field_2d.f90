!
!***********************************************************************
!
      subroutine dipole_field_2d
 
 
      use parameter_mod
      use functions_f90
      use erf_mod
!
!=======================================================================
!     
! Mirror dipole
!
      REAL :: x_dipole_mirror(2),y_dipole_mirror(2),z_dipole_mirror(2)
 
 
      REAL :: x_dipole, y_dipole , z_dipole ! position of the dipole
	  REAL :: r_dipole , th_dipole , phi 
	  REAL :: b_r , b_th, b_eta,b_x,b_y
	  REAL :: x_pos , y_pos , z_pos,x_rel,y_rel
	  INTEGER*8 :: i , j , k
 
      
         if (myid == 0) write(6,*) " Begin dipole_field_2d"
 
         x_dipole=hx*(i_dipole+0.5)
         y_dipole=hy*(j_dipole+0.5)
         z_dipole=hz*(k_dipole+0.5)
 
 
         x_dipole_mirror(1) =        x_dipole
         y_dipole_mirror(1) =        y_dipole
         z_dipole_mirror(1) =        z_dipole
         x_dipole_mirror(2) =       -x_dipole
         y_dipole_mirror(2) =        y_dipole
         z_dipole_mirror(2) =        z_dipole
 
 
         if (myid == 0) write(6,*) " x_dipole = ",x_dipole
         if (myid == 0) write(6,*) " y_dipole = ",y_dipole
         if (myid == 0) write(6,*) " z_dipole = ",z_dipole
 
! Set up bdipole_x
 
          if (myid == 0) then
            write(6,*) " Begin set up bdipole_2d "
            write(6,*) " dipole_moment = ",dipole_moment
          endif
	  do k = kb-1,ke+1
	    z_pos = float(k-1)*hz
		do j = jb-1,je+1
		  y_pos = float(j-1)*hy
		  do i = 1,nx2
		    x_pos = float(i-2)*hx
                    bdipole_x(i,j,k)=0.
                    bdipole_y(i,j,k)=0.
                    bdipole_z(i,j,k)=0.
                      if (image_dipole) then
                        do i_array=1,2
		          x_rel    =        x_pos-x_dipole_mirror(i_array)
		          y_rel    =        y_pos-y_dipole_mirror(i_array)
		          r_dipole = sqrt(x_rel**2+y_rel**2)
                          b_x = -2.*dipole_moment*x_rel*y_rel/r_dipole**4
                          b_y =     dipole_moment*(x_rel**2-y_rel**2)/r_dipole**4
                          bdipole_x(i,j,k) = bdipole_x(i,j,k)+b_x
                          bdipole_y(i,j,k) = bdipole_y(i,j,k)+b_y
                        enddo
                      else
		          x_rel    =        x_pos-x_dipole_mirror(1)
		          y_rel    =        y_pos-y_dipole_mirror(1)
		          r_dipole = sqrt(x_rel**2+y_rel**2)
                          b_x = -2.*dipole_moment*x_rel*y_rel/r_dipole**4
                          b_y =     dipole_moment*(x_rel**2-y_rel**2)/r_dipole**4
                          bdipole_x(i,j,k) = bdipole_x(i,j,k)+b_x
                          bdipole_y(i,j,k) = bdipole_y(i,j,k)+b_y
                      endif
		  enddo
		enddo
          enddo
          call xrealbcc_2d(bdipole_x,1,nx,ny,nz)
          call xrealbcc_2d(bdipole_y,1,nx,ny,nz)
          call xrealbcc_2d(bdipole_z,1,nx,ny,nz)
 
 
          bdipole_z=bdipole_z+bz_IMF/wpiwci
 
 
      dipole_sphere_ex=1.
      dipole_sphere_ey=1.
      dipole_sphere_ez=1.
      do k = kb-1,ke+1
        z_pos = float(k-1)*hz+0.5*hz
	do j = jb-1,je+1
	  y_pos = float(j-1)*hy+0.5*hy
	  do i = 1,nx2
	    x_pos = float(i-2)*hx+0.5*hx
            r_dipole = sqrt((x_pos-x_dipole)**2+(y_pos-y_dipole)**2)
            if (r_dipole < dipole_sphere_radius) then
              dipole_sphere_ex(i,j,k)=0.
              dipole_sphere_ey(i,j,k)=0.
              dipole_sphere_ez(i,j,k)=0.
            endif
          enddo
        enddo
      enddo
      call xrealbcc_2d(dipole_sphere_ex,1,nx,ny,nz)
      call xrealbcc_2d(dipole_sphere_ey,1,nx,ny,nz)
      call xrealbcc_2d(dipole_sphere_ez,1,nx,ny,nz)


      return
      end
