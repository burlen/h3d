!***************************************************************************
!                                                                          *
!                                 VERSION 6.0                              *
!                           YURI'S NONUNIFORM MESH                         *
!                           3D IMPLEMENTATION ONLY                         *
!                      UNIFORM LOADING IN PHYSICAL SPACE                   *
!               UNIFORM LOADING IN LOGICAL SPACE NOT YET IMPLEMENTED       *
!                                                                          *
!***************************************************************************

      program hybrid
 
      use parameter_mod
      use functions_f90
      use MESH2D
      integer*8:: time_begin(8),time_end(8),time_per_cycle,input_error,is
      integer*8 clock_time_re1, itstart, itfinish
      integer*8, pointer :: itest
      double precision, dimension(:,:,:), allocatable:: uniform_mesh,nonuniform_mesh_global
      character (len=240)::fileNameX
      integer:: iwrite
      integer*8 file_unit_101
      external get_environment_variable
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  Declaration
!
!----------------------------------------------------------------------
      integer :: isrtx,   & ! Extent of the grid per processor
                 iendx,   & !
                 isrty,   & !
                 iendy,   & !
                 isrtz,   & !
                 iendz
      integer :: insitu_now ! Perform insitu for this timestep?
!----------------------------------------------------------------------
! end Declaration
!----------------------------------------------------------------------


!
!***********************************************************************
!
       namelist /datum/ denmin,resis,iterb,testorbt,norbskip                &
       ,restart,nspec,nx,xmax,ny,ymax,nz,zmax,npx,npy,npz                   &
       ,dt,nprint,nwrtdata,qspec,wspec,nskipx                               &
       ,nskipy,nskipz,bxc,byc,bzc,frac,vxbar,vybar,vzbar,anisot             &
       ,btspec,bete,fxsho,nxcel,phib,netax,netay,netaz,ave1,ave2            &
       ,wpiwci,theta_max,restrt_write,Yee,b_in_xy,global,dtwci,dipole_moment&
       ,dipole_sphere_radius,i_dipole,j_dipole,k_dipole,absorbing_dipole    &
       ,harris,rcorr,ishape,teti,nwrtrestart,ieta,etamin,etamax             &
       ,R_MP,R_obstacle_to_MP,eta_par,bz_IMF,gama,QUOTA                     &
       ,maximum_simulation_time,n_subcycles,buffer_zone                     &
       ,xaa,xbb,nax,nbx,yaa,ybb,nay,nby,zaa,zbb,naz,nbz,setup_mesh          &
       ,profile_power,uniform_loading_in_logical_grid,moat_zone             &
       ,MPI_IO_format,dipole_ramp_time,image_dipole,smoothing               &
       ,time_reverse_B,reverse_B_ramp_time,smooth_coef,post_process         &
       ,xbox_l,xbox_r,ybox_l,ybox_r,zbox_l,zbox_r,t_begin,t_end,nwrtparticle &
       ,insitu

       parameter(file_unit_101 = 101)
       time_elapsed=0.;time_begin_array=0;time_end_array=0
       buffer_zone=0.
       notime=1
       bxc=0.;byc=0.;bzc=0.
       insitu = .true.
!
!***********************************************************************
!
!     MPI initialization
!
      call MPI_INIT(IERR)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

!     time stamp
!      call date_and_time(values=wall_clock_begin)
      initial_time=MPI_Wtime()
!
!***********************************************************************
!
      if (myid == 0) write(6,*) " 3DHYB is starting"
!
!***********************************************************************
!
!     get the i/o data directory name from the environment variable
!     DATA_DIRECTORY
!
      if (myid==0) then
        call get_environment_variable1(data_directory,len(data_directory))
        data_directory=trim(adjustl(data_directory))//'/'
        call get_environment_variable2(restart_directory,len(restart_directory))
        restart_directory=trim(adjustl(restart_directory))//'/'
        write(6,*)
        write(6,*) " I/O DATA_DIRECTORY = ",trim(adjustl(data_directory))
        write(6,*) " RESTART_DIRECTORY  = ",trim(adjustl(restart_directory))
        write(6,*)
      endif
      call MPI_BCAST(data_directory,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(restart_directory,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      restart_index_suffix(1)='.1'
      restart_index_suffix(2)='.2'

      my_short_int=myid
      call integer_to_character(myid_char,len(myid_char),my_short_int)
      if (myid_char == '') myid_char='0'
!
!***********************************************************************
!
!  set default values
!
      iwt=0;nskipx=1;nskipy=1;nskipz=1;testorbt=.false.;pi=acos(-1.d+00);frac=1.d+00;t_stopped=0.
!
!***********************************************************************
!
      if (myid == 0) then
!         open ( 5,file='finput.dat',form='formatted',status='old')
         open ( 5,file=trim(adjustl(data_directory))//'finput.dat',form='formatted',status='old')
         read(5,nml=datum,iostat=input_error)
         write(6, datum)
      endif
      if (myid == 0) then
        write(6,*) " VXBAR ==== ",vxbar
      endif
 
! hxv - 12/02/2008 -Automatic restart
      inquire(file=trim(adjustl(restart_directory))//'restart_index.dat',exist=restart)

      call MPI_BCAST(denmin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(resis                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(iterb                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(testorbt               ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(MPI_IO_format          ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(norbskip               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(restart                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nspec                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nx                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ny                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nz                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(xmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ymax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(zmax                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(npx                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(npy                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(npz                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dt                     ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nprint                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nwrtdata               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(qspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(wspec                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nskipx                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nskipy                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nskipz                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(bxc                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(byc                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(bzc                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(frac                   ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(vxbar                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(vybar                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(vzbar                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(anisot                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(gama                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(btspec                 ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(bete                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(fxsho                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nxcel                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(phib                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(netax                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(netay                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(netaz                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ave1                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ave2                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(wpiwci                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(theta_max              ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(restrt_write           ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(b_in_xy                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(global                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dtwci                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dipole_moment          ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dipole_sphere_radius   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(i_dipole               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(j_dipole               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(k_dipole               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(absorbing_dipole       ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(rcorr                  ,5     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ishape                 ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(teti                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nwrtrestart            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ieta                   ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(etamin                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(etamax                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(R_MP                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(R_obstacle_to_MP       ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(eta_par                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(bz_IMF                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(b_in_xy                ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(global                 ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(absorbing_dipole       ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(R_MP                   ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(R_obstacle_to_MP       ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nwrtrestart            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(QUOTA                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(maximum_simulation_time,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(n_subcycles            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(buffer_zone            ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(xaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(xbb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nax                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nbx                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(yaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ybb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nay                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nby                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(zaa                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(zbb                    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(naz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nbz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(setup_mesh             ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(profile_power          ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(uniform_loading_in_logical_grid,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(moat_zone              ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dipole_ramp_time       ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(image_dipole           ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(smoothing              ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(reverse_B_ramp_time    ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(time_reverse_B         ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(smooth_coef            ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(post_process           ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(xbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(xbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ybox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ybox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(zbox_l                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(zbox_r                 ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(t_begin                ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(t_end                  ,1     ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nwrtparticle           ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  Broadcast insitu logical
!
!----------------------------------------------------------------------
      call MPI_BCAST(insitu,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
!----------------------------------------------------------------------
! end Broadcast
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  Initialize
!
!    coprocessorinitializewithpython(python insitu script filename,
!                                    length of python insitu script filename)
!
!----------------------------------------------------------------------
      if (insitu) then
         call coprocessorinitializewithpython("insitu.py", 9)
      endif
!----------------------------------------------------------------------
! end Initialize
!----------------------------------------------------------------------
!
!

!
! 
      theta = theta_max
!
!
       dipole_ramp_time    = dipole_ramp_time*wpiwci
       time_reverse_B      = time_reverse_B*wpiwci
       reverse_B_ramp_time = reverse_B_ramp_time*wpiwci
       vxbar=vxbar/wpiwci
       vybar=vybar/wpiwci
       vzbar=vzbar/wpiwci
       dt=dtwci*wpiwci
       if (nz == 1) then
         dipole_moment_max=sqrt(2.)*vxbar(1)*(R_MP**2)
       else
         dipole_moment_max=sqrt(2.)*vxbar(1)*(R_MP**3)
       endif
       dipole_sphere_radius=R_MP*R_obstacle_to_MP

       if (myid==0) then
         write(6,*) " dipole_moment_max    = ",dipole_moment_max
         write(6,*) " dipole_sphere_radius = ",dipole_sphere_radius
       endif

! hxv - 12/02/2008 -Automatic restart
       if (restart .and. myid == 0) then
         write(6,*) " "
         write(6,*) " "
         write(6,*) " RUN IS RESTARTED FROM "//trim(adjustl(restart_directory))
         write(6,*) " "
         write(6,*) " "
       endif
!
!      field subcycling
!
       n_subcycles=max(n_subcycles,1)
       dt_field=dt/n_subcycles
!
!***********************************************************************
!
!     set MPI Cartesian geometry, define stride vector types, obtain new
!     ID for the processors, perform 2D decomposition of the
!     computational mesh, and find nearest neighbors (in y and z
!     directions)
!
      if (MYID.EQ.0) then
        write(6,*) " # of processors available= ",NUMPROCS
      endif
      if (nz == 1.and.ny == 1) then
        ndim=1
        dims(1)=1
        dims(2)=1
      else if (nz == 1) then
        ndim=1
        dims(2)=1
      else
        ndim=2
      endif
      PERIODS=.TRUE. 
      REORDER=.TRUE.
      call MPI_DIMS_CREATE(NUMPROCS,NDIM,DIMS,IERR)
      npy=npy/dims(1)
      npz=npz/dims(2)
      if (myid == 0) then
        do i=1,ndim
          write(6,*) " DIMENSION = ",I," DIMS = ",DIMS(I)
        enddo
      endif
      call MPI_CART_CREATE(MPI_COMM_WORLD,NDIM,DIMS,PERIODS,REORDER,COMM2D,IERR)
      call MPI_COMM_RANK(COMM2D,MYID,IERR)
      call MPI_CART_GET(COMM2D,NDIM,DIMS,PERIODS,COORDS,IERR)
      call MPE_DECOMP1D(NZ,DIMS(2),COORDS(2),KB,KE)
      call MPE_DECOMP1D(NY,DIMS(1),COORDS(1),JB,JE)
 
 
      nspecm = nspec
      nxmax  = nx+2
      nymax  = ny+2
      nzmax  = nz+2
      nylmax = je-jb+1
      nzlmax = ke-kb+1
      call set_parameters(NUMPROCS)
      myid_stop(myid)=0 
      do is = 1 , nspecm
         qleft(is)=0
         qrite(is)=0
      enddo

      if (MYID == 0.) then
        write(6,*) " LOCAL ARRAY SIZE IN Y-DIRECTION = ",JE-JB+1
        write(6,*) " LOCAL ARRAY SIZE IN Z-DIRECTION = ",KE-KB+1
      endif


      if (nzlmax < ke-kb+1) then
          print *,'myid = ',myid,' nzlmax lt ke-kb+1'
          print *,'myid = ',myid,' nzlmax,ke,kb= ',nzlmax,ke,kb
          myid_stop(myid) = 1
      endif
      if (nylmax < je-jb+1) then
          print *,'myid = ',myid,' nylmax lt je-jb+1'
          print *,'myid = ',myid,' nylmax,je,jb= ',nylmax,je,jb
          myid_stop(myid) = 1
      endif
      do i = 0, npes-1
!         IF (myid.eq.i) THEN
!            write(*,*)"Node no:",i,"myid_stop=",MYID_STOP(I)
!         ENDIF
         i_i = i
         CALL MPI_BCAST(MYID_STOP(i),1,MPI_INTEGER8,i_i,MPI_COMM_WORLD, IERR)
      enddo
      do i = 0,npes-1
         if (myid_stop(i).ne.0) then
            call MPI_FINALIZE(IERR)
            write(*,*)"TEST HERE"
            write(*,*)i, myid_stop(i)
            STOP
         endif
      enddo


!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  CreateVTKUniformGrid - At this point we can define the uniform grid
!
!    for 2D - x-dim is unchanged and y-dim is decomposed
!    for 3D - x-dim is unchanged and both y-dim and z-dim are 
!             decomposed
!
!----------------------------------------------------------------------
      if (insitu) then
         isrtx = -1
         iendx = nx
         isrty = jb-2
         iendy = jb+nylmax-1
         isrtz = kb-2
         iendz = kb+nzlmax-1
         call createvtkuniformgrid(isrtx,iendx,isrty,iendy,isrtz,iendz)
      endif
!----------------------------------------------------------------------
! end CreateVTKUniformGrid
!----------------------------------------------------------------------

!  Use CART_SHIFT to determine processor to immediate left
! (NBRLEFT) and right (NBRRITE) of processor MYID
!  Since code is aperiodic in z, need to manually set the
!  left boundary for processor 0 and right boundary for npes-1
 
 
      if (ndim == 2) then
        call MPI_CART_SHIFT(COMM2D,0,1,NBRLEFT,NBRRITE,IERR)
        call MPI_CART_SHIFT(COMM2D,1,1,NBRBOT ,NBRTOP ,IERR)
      else if (ndim == 1) then
        call MPI_CART_SHIFT(COMM2D,0,1,NBRLEFT,NBRRITE,IERR)
        NBRTOP=MYID
        NBRBOT=MYID
      else if (ndim == 0) then
        NBRLEFT=MYID
        NBRRITE=MYID
        NBRTOP =MYID
        NBRBOT =MYID
      endif
      call MPI_SENDRECV(NBRTOP    ,1,MPI_INTEGER ,NBRRITE,0,&
                        NBRLEFTTOP,1,MPI_INTEGER ,NBRLEFT,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(NBRTOP    ,1,MPI_INTEGER ,NBRLEFT,0,&
                        NBRRITETOP,1,MPI_INTEGER ,NBRRITE,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(NBRBOT    ,1,MPI_INTEGER ,NBRRITE,0,&
                        NBRLEFTBOT,1,MPI_INTEGER ,NBRLEFT,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(NBRBOT    ,1,MPI_INTEGER ,NBRLEFT,0,&
                        NBRRITEBOT,1,MPI_INTEGER ,NBRRITE,0,&
                        mpi_comm_world,status,ierr)
      nbrid(1)=NBRRITETOP
      nbrid(2)=NBRTOP
      nbrid(3)=NBRLEFTTOP
      nbrid(4)=NBRLEFT
      nbrid(5)=NBRLEFTBOT
      nbrid(6)=NBRBOT
      nbrid(7)=NBRRITEBOT
      nbrid(8)=NBRRITE
 
 
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        isendid(1)=1
      else
        isendid(1)=0
      endif
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,1)=nbrrite
        irecvid(2,1)=-1
        irecvid(3,1)=nbrleft
        irecvid(4,1)=-1
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,1)=-1
        irecvid(2,1)=nbrtop
        irecvid(3,1)=-1
        irecvid(4,1)=nbrbot
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,1)=nbrritetop
        irecvid(2,1)=nbrlefttop
        irecvid(3,1)=nbrleftbot
        irecvid(4,1)=nbrritebot
      endif
 
 
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        isendid(2)=1
      else
        isendid(2)=0
      endif
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,2)=nbrrite
        irecvid(2,2)=-1
        irecvid(3,2)=nbrleft
        irecvid(4,2)=-1
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,2)=-1
        irecvid(2,2)=nbrtop
        irecvid(3,2)=-1
        irecvid(4,2)=nbrbot
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,2)=nbrritetop
        irecvid(2,2)=nbrlefttop
        irecvid(3,2)=nbrleftbot
        irecvid(4,2)=nbrritebot
      endif
 
 
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        isendid(3)=1
      else
        isendid(3)=0
      endif
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,3)=nbrrite
        irecvid(2,3)=-1
        irecvid(3,3)=nbrleft
        irecvid(4,3)=-1
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,3)=-1
        irecvid(2,3)=nbrtop
        irecvid(3,3)=-1
        irecvid(4,3)=nbrbot
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,3)=nbrritetop
        irecvid(2,3)=nbrlefttop
        irecvid(3,3)=nbrleftbot
        irecvid(4,3)=nbrritebot
      endif
 
 
      if (mod(coords(1)+1,2) == 0.and.mod(coords(2)+1,2) == 0) then
        isendid(4)=1
      else
        isendid(4)=0
      endif
      if (mod(coords(1)  ,2) == 0.and.mod(coords(2)+1,2) == 0) then
        irecvid(1,4)=nbrrite
        irecvid(2,4)=-1
        irecvid(3,4)=nbrleft
        irecvid(4,4)=-1
      else if (mod(coords(1)+1,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,4)=-1
        irecvid(2,4)=nbrtop
        irecvid(3,4)=-1
        irecvid(4,4)=nbrbot
      else if (mod(coords(1)  ,2) == 0.and.mod(coords(2)  ,2) == 0) then
        irecvid(1,4)=nbrritetop
        irecvid(2,4)=nbrlefttop
        irecvid(3,4)=nbrleftbot
        irecvid(4,4)=nbrritebot
      endif
      nzl=nzlmax
      nyl=nylmax
      if (myid == 0) then
        jbglobal(myid)=jb
        jeglobal(myid)=je
        kbglobal(myid)=kb
        keglobal(myid)=ke
        do ipe=numprocs-1,1,-1
          call MPI_IRECV(jbglobal(ipe),1,MPI_INTEGER8,IPE,0,MPI_COMM_WORLD,req(1),IERR)
          call MPI_IRECV(jeglobal(ipe),1,MPI_INTEGER8,IPE,1,MPI_COMM_WORLD,req(2),IERR)
          call MPI_IRECV(kbglobal(ipe),1,MPI_INTEGER8,IPE,2,MPI_COMM_WORLD,req(3),IERR)
          call MPI_IRECV(keglobal(ipe),1,MPI_INTEGER8,IPE,3,MPI_COMM_WORLD,req(4),IERR)
          call MPI_WAITALL(4,req,status_array,IERR)
        enddo
      else
        call MPI_ISEND(jb           ,1,MPI_INTEGER8,0,0,MPI_COMM_WORLD,req(1),IERR)
        call MPI_ISEND(je           ,1,MPI_INTEGER8,0,1,MPI_COMM_WORLD,req(2),IERR)
        call MPI_ISEND(kb           ,1,MPI_INTEGER8,0,2,MPI_COMM_WORLD,req(3),IERR)
        call MPI_ISEND(ke           ,1,MPI_INTEGER8,0,3,MPI_COMM_WORLD,req(4),IERR)
        call MPI_WAITALL(4,req,status_array,IERR)
      endif
      call MPI_BCAST(JBGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(JEGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(KBGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(KEGLOBAL,NUMPROCS,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      if (myid.ne.0) then
        do k=kb,ke
          do j=1,je-jb+1
            jvec(j)=myid
          enddo
          i_length=je-jb+1
          call MPI_ISEND(jvec(1),i_length,MPI_INTEGER8,0,0,MPI_COMM_WORLD,req(1),IERR)
          call MPI_WAITALL(1,req,status_array,IERR)
        enddo
      else
        do k=kbglobal(myid),keglobal(myid)
          do j=jbglobal(myid),jeglobal(myid)
            idmap_yz(j,k)=myid
          enddo
        enddo
        do ipe=1,numprocs-1
          jbt=jbglobal(ipe)
          jet=jeglobal(ipe)
          kbt=kbglobal(ipe)
          ket=keglobal(ipe)
          do k=kbt,ket
            i_length=jet-jbt+1
            call MPI_IRECV(idmap_yz(jbt,k),i_length,MPI_INTEGER8,IPE,0,MPI_COMM_WORLD,req(1),IERR)
            call MPI_WAITALL(1,req,status_array,IERR)
          enddo
        enddo
      endif
      call MPI_BCAST(idmap_yz,size(idmap_yz),MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
      idmap(0)=idmap(1)
      call MPI_TYPE_VECTOR(nzl+2,nx+2,(nx+2)*(nyl+2),MPI_DOUBLE_PRECISION,stridery,IERR)
      call MPI_TYPE_COMMIT(stridery,IERR)
      call MPI_TYPE_VECTOR(nyl+2,nx+2,nx+2          ,MPI_DOUBLE_PRECISION,STRIDERZ,IERR)
      call MPI_TYPE_COMMIT(STRIDERZ,IERR)
 
 
      nptotp=0
      do is=1,nspec
        nptotp=nptotp+npx(is)*npy(is)*npz(is)
      enddo
      if (nptotp > nplmax) then
        if (myid == 0) then
          write(6,*) ' Increase nplmax in the input file '
          write(6,*) 'nptotp = ',nptotp
        endif
        myid_stop(myid) = 1
      endif
      do i = 0, npes-1
         i_i = i
         CALL MPI_BCAST(MYID_STOP(i),1,MPI_INTEGER8,i_i,MPI_COMM_WORLD, IERR)
      enddo
      do i = 0,npes-1
         if (myid_stop(i).ne.0) then
            call MPI_FINALIZE(IERR)
            write(*,*)"TEST HERE"
            STOP
         endif
      enddo
      if (.not.testorbt) norbskip=1
!
!=======================================================================
!
      call allocate_global_arrays
      call pdf_injection

!
!=======================================================================
!
!     Initialize nonuniform mesh
!
      call MESH_INIT(meshX,xaa,xbb,xmax,nax,nbx,nx) ! initialize x-mesh
      call MESH_INIT(meshY,yaa,ybb,ymax,nay,nby,ny) ! initialize y-mesh
      call MESH_INIT(meshZ,zaa,zbb,zmax,naz,nbz,nz) ! initialize z-mesh

      call MESH_INDEX(meshX,CELL,ixv_2_c_map)
      call MESH_INDEX(meshY,CELL,iyv_2_c_map)
      call MESH_INDEX(meshZ,CELL,izv_2_c_map)
      call MESH_INDEX(meshX,NODE,ixv_2_v_map)
      call MESH_INDEX(meshY,NODE,iyv_2_v_map)
      call MESH_INDEX(meshZ,NODE,izv_2_v_map)
      call MESH_INDEX(meshX,CELL,ixc_2_c_map,CELL)
      call MESH_INDEX(meshY,CELL,iyc_2_c_map,CELL)
      call MESH_INDEX(meshZ,CELL,izc_2_c_map,CELL)
      call MESH_INDEX(meshX,NODE,ixc_2_v_map,CELL)
      call MESH_INDEX(meshY,NODE,iyc_2_v_map,CELL)
      call MESH_INDEX(meshZ,NODE,izc_2_v_map,CELL)

      if (myid == 0) then
!
!        write(6,*) " nl_x = ",meshX%nl
!        do i=1,meshX%nl+3 
!          write(6,*) " i, x;",i,meshX%xn(i),meshX%xc(i)
!        enddo
!        write(6,*) " nl_y = ",meshY%nl
!        do i=1,meshY%nl+3 
!          write(6,*) " i, y;",i,meshY%xn(i),meshY%xc(i)
!        enddo
!        write(6,*) " nl_z = ",meshZ%nl
!        do i=1,meshZ%nl+3 
!          write(6,*) " i, z;",i,meshZ%xn(i),meshZ%xc(i)
!        enddo

        open(unit=11,file='mesh_vertices.dat',status='unknown',form='formatted')

        write(11,*) meshX%nl+1,meshY%nl+1,meshZ%nl+1
        do i=2,meshX%nl+2 
          write(11,*) meshX%xn(i)
        enddo
        do i=2,meshX%nl+2 
          write(11,*) meshX%dxc(i)
        enddo
!
        do i=2,meshY%nl+2 
          write(11,*) meshY%xn(i)
        enddo
        do i=2,meshY%nl+2 
          write(11,*) meshY%dxc(i)
        enddo
!
        do i=2,meshZ%nl+2 
          write(11,*) meshZ%xn(i)
        enddo
        do i=2,meshZ%nl+2 
          write(11,*) meshZ%dxc(i)
        enddo
!         
        close(unit=11)
      endif

!     STOP HERE IF SETUP_MESH=.TRUE.

      if (setup_mesh) then
        call MPI_FINALIZE(IERR)
        STOP
      endif

      allocate (uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1))
      allocate (nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1))
!
!=======================================================================
!
      call date_and_time(values=time_begin)
      clock_time_re1=(time_begin(5)*3600.+time_begin(6)*60.+time_begin(7)+time_begin(8)*0.001)
 
!      if (myid == 0) inquire(file='restart_index.dat',exist=restart)
      if (myid == 0) inquire(file=trim(adjustl(restart_directory))//'restart_index.dat',exist=restart)
      call MPI_BCAST(restart         ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
   
      if (restart) then
 
        call makelist

        if (myid == 0) then
!          open(unit=222,file='restart_index.dat' ,status='old')
          open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='old')
          read(222,*) restart_index,itfin
          close(222)
        endif

        call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
        nplmax6 = 6*nplmax
! hxv 01/10/2014
!
        if (myid == 0) then
          write(6,*) " "
          write(6,*) " RESTARTED FROM SET # ",restart_index
          write(6,*) " "
        endif

 
!       comment out for timing on LANL machine

        do iwrite = 0,npes_over_60 
         if (mod(myid,npes_over_60 + 1).eq.iwrite) then
           call restrtrw(-1.0,itstart)
           call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
         endif
!         call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        enddo

        if (restart_index == 1) then
          restart_index=2
        else
          restart_index=1
        endif

!       Uniform mesh - Same as is in version 5.0
        yb=(jb-1)*hy
        ye= je   *hy
        zb=(kb-1)*hz
        ze= ke   *hz

!       Nonuniform mesh
        zb=meshZ%xn(kb+1)
        ze=meshZ%xn(ke+2)
        do ipe=0,npes-1
          zbglobal(ipe)=meshZ%xn(kbglobal(ipe)+1)
          zeglobal(ipe)=meshZ%xn(keglobal(ipe)+2)
        enddo
        yb=meshY%xn(jb+1)
        ye=meshY%xn(je+2)
        do ipe=0,npes-1
          ybglobal(ipe)=meshY%xn(jbglobal(ipe)+1)
          yeglobal(ipe)=meshY%xn(jeglobal(ipe)+2)
        enddo


        volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)

        xb        = 0.
        xe        = xmax
        xb_logical=MESH_UNMAP(meshX,xb)
        xe_logical=MESH_UNMAP(meshX,xe)
        yb_logical=MESH_UNMAP(meshY,yb)
        ye_logical=MESH_UNMAP(meshY,ye)
        zb_logical=MESH_UNMAP(meshZ,zb)
        ze_logical=MESH_UNMAP(meshZ,ze)


        do is=1,nspec
          npm=npx(is)*npy(is)*npz(is)*npes
          dfac(is)=real(ny*nz*nx)/real(npm)
          do ixe=1,nx2
            do iye=jb-1,je+1
              do ize=kb-1,ke+1
!                 qp_cell(ixe,iye,ize,is) = (meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)/(hx*hy*hz))*dfac(is)*frac(is)
                 qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
              enddo
            enddo
          enddo
        enddo


        do i=1,nxmax
          xc_uniform(i) = hx*(i-1.5)
          xv_uniform(i) = hx*(i-2.0)
        enddo
        do j=1,nymax
          yc_uniform(j) = hy*(j-0.5)
          yv_uniform(j) = hy*(j-1.0)
        enddo
        do k=1,nzmax
          zc_uniform(k) = hz*(k-0.5)
          zv_uniform(k) = hz*(k-1.0)
        enddo

        if (myid ==0) then
          my_short_int=it
          call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
          cycle_ascii_new=trim(adjustl(cycle_ascii))
          write(6,*) " cycle = ",cycle_ascii_new
        endif

        call MPI_BCAST(cycle_ascii    ,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

        call opendiagfiles

      else

        time=0.0
        call makelist
        if (myid == 0) write(6,*) " harris = ",harris
        if (myid == 0) write(6,*) " global = ",global
        if (myid == 0) write(6,*) " Yee    = ",yee
        restart_index=1
        if (dipole_ramp_time /= 0.) then
          dipole_moment = 0.
        else
          dipole_moment = dipole_moment_max
        endif
        call init_global

      endif
  
 !    Write dipole field to file

        irecnum = 1
        lenrec=(nxmax-2)*recl_for_real
!        call MESH_INTERPOLATED_3D(bdipole_x,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bdipole_x
        rnorm = wpiwci
        if (MPI_IO_format) then
	  fileNameX = trim(adjustl(data_directory))//'bx_dipole.gda'
          call wrtfile         (uniform_mesh,rnorm,fileNameX,irecnum,ny,nz)
        else
          irecnum = 1
          lenrec=(nxmax-2)*recl_for_real
          open (file_unit_101,                                                                                  &
                file= trim(adjustl(data_directory))//'bx_dipole.gda', &
                form='unformatted',                                                                   &
                action='write',access='direct', status='unknown',recl=lenrec)
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit_101,irecnum,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(bdipole_y,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bdipole_y
        if (MPI_IO_format) then
	  fileNameX = trim(adjustl(data_directory))//'by_dipole.gda'
          call wrtfile         (uniform_mesh,rnorm,fileNameX,irecnum,ny,nz)
        else
          open (file_unit_101,                                                                                  &
                file= trim(adjustl(data_directory))//'by_dipole.gda', &
                form='unformatted',                                                                   &
                action='write',access='direct', status='unknown',recl=lenrec)
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit_101,irecnum,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(bdipole_z,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bdipole_z
        if (MPI_IO_format) then
	  fileNameX = trim(adjustl(data_directory))//'bz_dipole.gda'
          call wrtfile         (uniform_mesh,rnorm,fileNameX,irecnum,ny,nz)
        else
          open (file_unit_101,                                                                                  &
                file= trim(adjustl(data_directory))//'bz_dipole.gda', &
                form='unformatted',                                                                   &
                action='write',access='direct', status='unknown',recl=lenrec)
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit_101,irecnum,ny,nz)
        endif

      call date_and_time(values=time_end)
      clock_time_init=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
      if (myid == 0) then
        print *,'load time = ',real(clock_time_init-clock_time_re1)  
      endif
      clock_time_old = clock_time_init
!
! hxv 01/10/2014
!
      if (myid == 0) then
        write(6,*) " RESTARTED AT CYCLE # ",itfin
      endif
!
!=======================================================================
!
! Write particle data to file and exit if post_process=.true.
!
      if (post_process) then
        if (myid == 0) then
          my_short_int=itfin
          call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
        endif
        call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        write(6,*) " calling particle_in_volume_write with cycle_ascii = ",cycle_ascii
        call particle_in_volume_write
        goto 999
      endif
!
!=======================================================================
!
!     march forward in time
!
!=======================================================================
!
      itstart = itfin+1

!     change how itfinish is computed
      itfinish = (maximum_simulation_time-t_stopped)/dtwci+itstart-1

      if (myid == 0) write(6,*) 'itstart, itfinish = ',itstart,' ',itfinish
      if (myid == 0) write(6,*) 'nwrtdata = ',nwrtdata
      it = itstart
      time_elapsed=0.;time_begin_array=0;time_end_array=0

      do while(it <= itfinish)
!
!=======================================================================
!
        call get_cleanup_status(len(cleanup_status))

        if (myid == 0) then
          write(6,*) " "
          write(6,*) " "
          WRITE(6,*) " DT = ",DT
          WRITE(6,*) " T_STOPPED = ",T_STOPPED
        endif

!       time stamp
!        call date_and_time(values=wall_clock_end)
!        wall_clock_elapsed= (wall_clock_end(3)-wall_clock_begin(3))*24.                          &
!                           +(wall_clock_end(5)-wall_clock_begin(5))                              &
!                           +(wall_clock_end(6)-wall_clock_begin(6))/60.                          &
!                           +(wall_clock_end(7)-wall_clock_begin(7))/3600.                        &
!                           +(wall_clock_end(8)-wall_clock_begin(8))*0.001/3600.
!        call date_and_time(values=wall_clock_end)
!=======================================================================
        final_time=MPI_Wtime()
        wall_clock_elapsed= final_time-initial_time
        if (wall_clock_elapsed >= QUOTA*3600. .and. myid == 0) then
          cleanup_status = 'CLEANUP_STATUS=TRUE'
          write(6,*) " "
          write(6,*) " "
          write(6,*) " INSUFFICIENT RUNTIME ALLOWED: DUMPING RESTART DATA"
          write(6,*) " MY QUOTA = ",QUOTA
          write(6,*) " "
          write(6,*) " "
        endif
 
!       if (it == 3) then
!       write(6,*) " myid,",myid," is stopped"
!       call MPI_FINALIZE(IERR)
!       STOP
!       endif
        
        call MPI_BCAST(cleanup_status,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

        if (cleanup_status == 'CLEANUP_STATUS=EXIT') then
          goto 999
        else if (cleanup_status == 'CLEANUP_STATUS=TRUE') then
          if (myid == 0) then
            WRITE(6,*)
            WRITE(6,*) 'WRITING THE RESTART FILE'
            WRITE(6,*)
          endif

          itfin = it
!         comment out for timing on LANL machine

        do iwrite = 0,npes_over_60 
         if (mod(myid,npes_over_60 + 1).eq.iwrite) then
           call restrtrw(1.0,itstart)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        enddo

          if (myid == 0) then
          open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='unknown')
            write(222,*) restart_index,itfin
            close(222)
          endif
          goto 998
        endif
!
!=======================================================================
!
        call date_and_time(values=time_begin_array(:,1))
!
!=======================================================================
!
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
!
        if (notime == 0) then
          write(file_unit_time,"(i4,' begin    ',e20.8)") it,real(clock_time-clock_time_init)
        endif
        if (myid == 0.and.mod(it,10) == 0) then
          print *,"it = ",it
          print *,'system time (delta) = ',real(clock_time - clock_time_old)
          print *,'system time (total) = ',real(clock_time - clock_time_init)
        endif
        clock_time_old = clock_time
!
!=======================================================================
!
!     determine whether or not to print diagnostic information
!
        if (mod(it,nprint) == 0) then
          prntinfo=.true.
        else
          prntinfo=.false.
        endif
!
!=======================================================================
!
!     determine whether or not to write data into files
!
        if (mod(it,nwrtdata) == 0) then
          wrtdat=.true.
        else
          wrtdat=.false.
        endif
!
!=======================================================================
!
! HXV
!
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        if (notime == 0) then
          write(file_unit_time,"(i4,' trans    ',e20.8)") it,real(clock_time-clock_time_init)
        endif
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        if (notime == 0) then
          write(file_unit_time,"(i4,' injctpar ',e20.8)") it,real(clock_time-clock_time_init)
        endif

        call date_and_time(values=time_begin_array(:,3))
        if (it <=1 .or. harris) goto 2
        do is=1,nspec
           call particle_newinject_linked_list(is)
        enddo

!        if (it == 2 ) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF

 2      continue
        call date_and_time(values=time_end_array(:,3))
         call date_and_time(values=time_begin_array(:,2))
         if (ndim /= 1) then
           call etacalc       ! Dietmar's resistivity
         else
           call etacalc_2d    ! Dietmar's resistivity
         endif


         call trans

!        if (it == 1 ) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF

         call date_and_time(values=time_end_array(:,2))
 

        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 

        if (notime == 0) then
          write(file_unit_time,"(i4,' sortit   ',e20.8)") it,real(clock_time-clock_time_init)
        endif

 
        call date_and_time(values=time_begin_array(:,4))
        if (mod(it,10) == 0) call sortit    !  sort the particles
        call date_and_time(values=time_end_array(:,4))
 
 
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 

        if (notime == 0) then
          write(file_unit_time,"(i4,' field    ',e20.8)") it,real(clock_time-clock_time_init)
        endif

        call date_and_time(values=time_begin_array(:,5))
        if (.not.testorbt) then
          if (ndim /=1) then 
            call field
          else
            call field_2d
          endif
        endif
 
        call date_and_time(values=time_end)
        call date_and_time(values=time_end_array(:,5))
 
        
 998    if (    wrtdat.or.cleanup_status == 'CLEANUP_STATUS=TRUE'.or.it == itfinish       &
            .or.it == 1) then
          if (myid ==0) then
            my_short_int=it
            call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
            cycle_ascii_new=trim(adjustl(cycle_ascii))
            write(6,*) " cycle = ",cycle_ascii_new
          endif
          call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
          call MPI_BCAST(cycle_ascii_new,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

         if (myid == 0 .and. .not. MPI_IO_format) then
            call openfiles
         endif
         
          call date_and_time(values=time_begin_array(:,6))
          if (ndim /= 1) then
            call caltemp2_global
          else
            call caltemp2_global_2d
          endif
          call date_and_time(values=time_end_array(:,6))
          numvars = 18
          irecnum=1
          call dataout(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,tperp,                          &
                      p_xx,p_xy,p_xz,p_yy,p_yz,p_zz,                                          &
                      nxmax,nymax,nzmax,file_unit,myid,                                       &
	              numvars,irecnum,kb,ke,numprocs,wpiwci,jb,je,ny,nz,nylmax,nzlmax,nspecm, &
                      bdipole_x,bdipole_y,bdipole_z,eta,eta_times_b_dot_j,eta_par,            &
                      uniform_mesh,nonuniform_mesh_global,trim(adjustl(data_directory)), trim(adjustl(cycle_ascii)),      &
                      MPI_IO_format)
         if (myid == 0 .and. .not. MPI_IO_format) then
           do j=1,20
             close(file_unit(j))
           enddo
         endif
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  CoProcessing - At this point we coprocess
!
!----------------------------------------------------------------------
      if (insitu) then
!        if (mod(it,11).eq.0) then
         my_short_int=it
         call requestdatadescription(my_short_int-1,                  &
                                     DBLE(my_short_int-1),insitu_now)
         call addattributes(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,      &
                            tperp,eta,nxmax,nylmax,nzlmax,               &
                            kb,ke,jb,je,nspecm)
         call coprocess() 
!        endif
      endif
!----------------------------------------------------------------------
! end CoProcessing
!----------------------------------------------------------------------

!
!=======================================================================
!
! Write particle data to file and exit if post_process=.true.
!
      if (myid == 0) then
        write(6,*) " t_begin,time,t_end = ",t_begin*wpiwci,time,t_end*wpiwci
        write(6,*) "it,nwrtparticle = ",it,nwrtparticle
      endif
      if (t_begin*wpiwci <= time .and. time <= t_end*wpiwci .and. mod(it,nwrtparticle)==0) then
        if (myid == 0) then
          my_short_int=it
          call integer_to_character(cycle_ascii,len(cycle_ascii),my_short_int)
        endif
        call MPI_BCAST(cycle_ascii,160,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        write(6,*) " calling particle_in_volume_write with cycle_ascii = ",cycle_ascii
        call particle_in_volume_write
      endif
!
!=======================================================================
!
          if (cleanup_status == 'CLEANUP_STATUS=TRUE') goto 999
        endif
!
!=======================================================================
!
!       maximum_simulation_time
!       if (restrt_write == 1.and.mod(it,nwrtrestart)==0) then
        if (restrt_write == 1.and.(mod(it,nwrtrestart)==0 .or. it == itfinish)) then

          if (myid == 0) then
            WRITE(6,*)
            WRITE(6,*) 'WRITING THE RESTART FILE'
            WRITE(6,*)
          endif

          itfin = it

        do iwrite = 0,npes_over_60  
         if (mod(myid,npes_over_60 + 1).eq.iwrite) then
           call restrtrw(1.0,itstart)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        enddo

          if (myid == 0) then
          open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='unknown')
            write(222,*) restart_index,itfin
            close(222)
          endif

          if (restart_index == 1) then
            restart_index=2
          else
            restart_index=1
          endif

        endif

!        if (it == 1 ) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF
!
!=======================================================================
!
        time=time+dt
        it = it + 1
!
!=======================================================================
!
         theta_old = theta
         if (time_reverse_B >= 0. .and. time >= time_reverse_B) then
           if (myid == 0) write(6,*) "it = ",it," time = ",time," time_reverse_B = ",time_reverse_B
           if (reverse_B_ramp_time == 0.) then
             theta = -theta_max
           else
             theta = theta_max*(1.-2.*min(1.,(time-time_reverse_B)/reverse_B_ramp_time))
           endif
         endif
         if (theta /= theta_old) then
           	bx(1,:,:)=bx(1,:,:)-bxc
           	by(1,:,:)=by(1,:,:)-byc
           	bz(1,:,:)=bz(1,:,:)-bzc
           	pi=acos(-1.d+00)
           	pifac=180./pi
           	bxc = cos(theta/pifac)/wpiwci
           	byc = sin(theta/pifac)/wpiwci
           	bzc = 0.
           	exc= (vzbar(1))*byc-(vybar(1))*(bzc+bz_IMF/wpiwci)
           	eyc=+(vxbar(1))*(bzc+bz_IMF/wpiwci)-(vzbar(1))*bxc
           	ezc= (vybar(1))*bxc-(vxbar(1))*byc
           	bx(1,:,:)=bx(1,:,:)+bxc
           	by(1,:,:)=by(1,:,:)+byc
           	bz(1,:,:)=bz(1,:,:)+bzc
           	ex(1,:,:)=exc
           	ey(1,:,:)=eyc
           	ez(1,:,:)=ezc
         endif
!
!       Linear ramp of dipole moment
!
        if (dipole_ramp_time /= 0.) then
          dipole_moment = dipole_moment_max*min(1.,time/dipole_ramp_time)
        else
          dipole_moment = dipole_moment_max
        endif

        if (ndim /=1) then
          call dipole_field
        else
          call dipole_field_2d
        endif
        if (myid == 0) then 
          write(6,*) " dipole moment = ",dipole_moment
          write(6,*) " max(bdipole_x) = ",maxval(abs(bdipole_x))
          write(6,*) " max(bdipole_y) = ",maxval(abs(bdipole_y))
          write(6,*) " max(bdipole_z) = ",maxval(abs(bdipole_z))
        endif
!
!=======================================================================
!
        call date_and_time(values=time_end_array(:,1))
!
!=======================================================================
!
        do j=1,6
          call accumulate_time_difference(time_begin_array(1,j),time_end_array(1,j),time_elapsed(j))
        enddo
!
!=======================================================================
!

!        if (it == 1 ) then
!          CALL MPI_FINALIZE(IERR)
!          STOP
!        ENDIF

      enddo    ! IT do while loop
      
 999  if (notime == 0) close(file_unit_time)
!
!***********************************************************************
!
      if (MYID.EQ.0) then
        write(6,*) " "
        write(6,*) " "
        write(6,*) " *** RUN COMPLETED *** RUN COMPLETED *** RUN COMPLETED "
        write(6,*) " "
        write(6,*) " "
        write(6,*) " subroutine trans             (s)          =",time_elapsed(2)
        write(6,*) " subroutine injectpar         (s)          =",time_elapsed(3)
        write(6,*) " subroutine sort              (s)          =",time_elapsed(4)
        write(6,*) " subroutine field             (s)          =",time_elapsed(5)
        write(6,*) " subroutine caltemp2          (s)          =",time_elapsed(6)
        write(6,*) " total time                   (s)          =",time_elapsed(1)
        write(6,*) " "
        write(6,*) " "
        write(6,*) " In subroutine caltemp2,"
        write(6,*) "   subroutine xreal           (s)          =",time_elapsed(24)
        write(6,*) "   subroutine xrealbcc        (s)          =",time_elapsed(25)
        write(6,*) "   interpolation              (s)          =",time_elapsed(26)
        write(6,*) "   total caltemp2             (s)          =",time_elapsed(23)
        write(6,*) " "
        write(6,*) " "
        write(6,*) " In subroutine trans," 
        write(6,*) "   subroutine parmov          (s)          =",time_elapsed(7)
        write(6,*) "   subroutine energy          (s)          =",time_elapsed(8)
        write(6,*) "   total trans                (s)          =",time_elapsed(20)
        write(6,*) " "
        write(6,*) " "
        write(6,*) " In subroutine field,"
        write(6,*) "   subroutine pressgrad       (s)          =",time_elapsed(9)
        write(6,*) "   subroutine bcalc           (s)          =",time_elapsed(10)
        write(6,*) "   subroutine ecalc           (s)          =",time_elapsed(11)
        write(6,*) "   subroutine focalc          (s)          =",time_elapsed(12)
        write(6,*) "   total field                (s)          =",time_elapsed(21)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "   In subroutine parmov,"
        write(6,*) "     particle pushing         (s)          =",time_elapsed(13)
        write(6,*) "     particle exchange        (s)          =",time_elapsed(14)
        write(6,*) "     particle interpolation   (s)          =",time_elapsed(15)
        write(6,*) "     total parmov             (s)          =",time_elapsed(19)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "   In subroutine bcalc,"
        write(6,*) "     subroutine ecalc         (s)          =",time_elapsed(16)
        write(6,*) "     subroutine xrealbcc      (s)          =",time_elapsed(17)
        write(6,*) "     total bcalc              (s)          =",time_elapsed(22)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "     In subroutine ecalc (called from subroutine bcalc and field),"
        write(6,*) "       subroutine xrealbcc    (s)          =",time_elapsed(18)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "communication time/total time (%)          =" &
        ,100.*(time_elapsed(24)+time_elapsed(25)+time_elapsed(14)+time_elapsed(17) &
        +time_elapsed(18))/time_elapsed(1)
        write(6,*) " "
        write(6,*) " "
        write(6,*) "Further breakdown of communication time "
        write(6,*) "  particle exchage in subroutine parmov (%) =" &
        ,100.*time_elapsed(14)/time_elapsed(1)
        write(6,*) "  subroutine xrealbcc                   (%) =" &
        ,100.*(time_elapsed(25)+time_elapsed(17)+time_elapsed(18))/time_elapsed(1)
        write(6,*) "  subroutine xreal                      (%) =" &
        ,100.*time_elapsed(24)/time_elapsed(1)
      endif
!
!=======================================================================
!
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  Finalize
!
!----------------------------------------------------------------------
      call coprocessorfinalize()
!----------------------------------------------------------------------
! end Finalize
!----------------------------------------------------------------------
!
      call MPI_FINALIZE(IERR)
      stop
      end

!***********************************************************************
!
      subroutine addattributes(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,   &
                               tperp,eta,nxmax,nylmax,nzlmax,            &
                               kb,ke,jb,je,nspecm)
!======================================================================
!
! insitu Patrick O'Leary 2/21/2013
!
!  addattributes - adding raw scalars and vectors to insitu pipeline.
!
!======================================================================
      implicit none
      integer*8 nxmax,nylmax,nzlmax,kb,ke,jb,je,nspecm
      double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: &
          bx,by,bz,den,ex,ey,ez,vix,viy,viz,eta
      double precision,                                         &
          dimension(nxmax,jb-1:je+1,kb-1:ke+1,nspecm) ::              &
          tpar,tperp
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  AddAttributes - At this point we add raw scalars, otherwise we would 
!               have create new memory for the REAL*4, and scale by
!               other scalars and a real*4 normalizer.
!
!----------------------------------------------------------------------
        call adduniformgridvector('b',1,bx,by,bz,                     &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('den',3,den,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridvector('e',1,ex,ey,ez,                     &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridvector('vi',2,vix,viy,viz,                 &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('tpar',4,tpar(:,:,:,1),             &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('tperp',5,tperp(:,:,:,1),           &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('eta',3,eta,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
!----------------------------------------------------------------------
! end AddAttributes
!----------------------------------------------------------------------
        return
      end subroutine addattributes


!
!***********************************************************************
!
      subroutine dataout( bx, by, bz, den, ex, ey, ez, vix, viy, viz, tpar, tperp,                &
                          p_xx,p_xy,p_xz,p_yy,p_yz,p_zz,                                          &
                          nxmax,nymax,nzmax,file_unit,myid,                                       &
                          numvars,irecnum,kb,ke,numprocs,wpiwci,jb,je,ny,nz,nylmax,nzlmax,nspecm, &
                          bdipole_x,bdipole_y,bdipole_z,eta,eta_times_b_dot_j,eta_par,            &
                          uniform_mesh,nonuniform_mesh_global,data_directory, cycle_ascii,MPI_IO_format)

      integer*8 nxmax, nymax, nzmax, myid, numvars, irecnum, kb, ke, numprocs,jb,je,nylmax,nzlmax,nspecm
      integer*8, dimension(20) :: file_unit
      double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: bx, by, bz, den,                  &
      ex, ey, ez, vix, viy, viz,bdipole_x,bdipole_y,bdipole_z,eta,eta_times_b_dot_j
      double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1,nspecm) :: tpar, tperp
      double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1,nspecm) :: p_xx,p_xy,p_xz,p_yy,p_yz,p_zz
      integer*8 ny,nz,eta_par
      double precision rnorm, wpiwci
      integer*8 irecdel, ir1, ir2, irec_start , idebug,IERR
!      double precision:: uniform_mesh(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)
      double precision:: uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1)
      double precision:: nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1)
      character(len=240):: fileName
      character data_directory*(*)
      character cycle_ascii*(*)
      logical:: MPI_IO_format

 
 
!      irecdel = (nz+2)*(ny+2)
      irecdel = nz*ny
 
! Now determine the starting and ending record number 
! for a given variable and myid value.
 

      irec_start = irecnum
      rnorm = wpiwci
      bx=bx+bdipole_x
!      call MESH_INTERPOLATED_3D(bx,uniform_mesh,nonuniform_mesh_global)
      uniform_mesh=bx 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'bx_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(1),irec_start,ny,nz)
        endif
      bx=bx-bdipole_x
      idebug = 0
      if (idebug.ne.1) then
        rnorm = wpiwci
        by=by+bdipole_y
!        call MESH_INTERPOLATED_3D(by,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=by 
        if (MPI_IO_format) then
 	  fileName= trim(trim(adjustl(data_directory))//'by_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(2),irec_start,ny,nz)
        endif
        by=by-bdipole_y
        rnorm = wpiwci
        bz=bz+bdipole_z
!        call MESH_INTERPOLATED_3D(bz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bz
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'bz_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(3),irec_start,ny,nz)
        endif
        bz=bz-bdipole_z
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(den,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=den
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'den_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(4),irec_start,ny,nz)
        endif
        rnorm = wpiwci**2
!        call MESH_INTERPOLATED_3D(ex,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=ex 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'ex_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(5),irec_start,ny,nz)
        endif
        rnorm = wpiwci**2
!        call MESH_INTERPOLATED_3D(ey,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=ey 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'ey_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(6),irec_start,ny,nz)
        endif
        rnorm = wpiwci**2
!        call MESH_INTERPOLATED_3D(ez,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=ez 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'ez_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(7),irec_start,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(vix,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=vix
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'vix_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(8),irec_start,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(viy,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=viy
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'viy_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(9),irec_start,ny,nz)
        endif
        rnorm = wpiwci 
!        call MESH_INTERPOLATED_3D(viz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=viz
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'viz_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(10),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(tpar,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=tpar(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'tpar_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(11),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(tperp,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=tperp(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'tperp_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(12),irec_start,ny,nz)
        endif


        rnorm = 1.
!        call MESH_INTERPOLATED_3D(p_xx,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=p_xx(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'p_xx_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(14),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(p_xy,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=p_xy(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'p_xy_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(15),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(p_xz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=p_xz(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'p_xz_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(16),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(p_yy,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=p_yy(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'p_yy_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(17),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(p_yz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=p_yz(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'p_yz_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(18),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(p_zz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=p_zz(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'p_zz_'//trim(adjustl(cycle_ascii)))//'.gda'
          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(19),irec_start,ny,nz)
        endif




        if (eta_par == 0) then
          rnorm = 1.
!          call MESH_INTERPOLATED_3D(eta,uniform_mesh,nonuniform_mesh_global)
          uniform_mesh=eta
          if (MPI_IO_format) then
            fileName= trim(trim(adjustl(data_directory))//'eta_'//trim(adjustl(cycle_ascii)))//'.gda'
            call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
          else
            call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(13),irec_start,ny,nz)
          endif
        else
          rnorm = 1.
!          call MESH_INTERPOLATED_3D(eta_times_b_dot_j,uniform_mesh,nonuniform_mesh_global)
          uniform_mesh=eta_times_b_dot_j
          if (MPI_IO_format) then
            fileName= trim(trim(adjustl(data_directory))//'eta_par_'//trim(adjustl(cycle_ascii)))//'.gda'
            call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
          else
            call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(13),irec_start,ny,nz)
          endif
        endif
      endif
      irecnum = irecnum + irecdel

      return
      end subroutine dataout
!
!***********************************************************************
!
      subroutine makelist
 
      use parameter_mod
      integer*8:: ip
 
      ipstore=1
      ipleft =0
      iprite =0
      iprecv =0
      iphead =0
      iptemp =0
      do ip=1,nplmax-1
        link(ip)=ip+1
      enddo
      link(nplmax)=0
      return
      end
!
!***********************************************************************
!
      subroutine openfiles
 
      use parameter_mod
      integer*8  file_unit_ref
 
      if (.not.testorbt) then
!        lenrec=nxmax*recl_for_real
        lenrec=(nxmax-2)*recl_for_real
        file_unit_ref = 250
        do j=1,20
  	  file_unit(j) = file_unit_ref + j
        enddo
        open (file_unit(1),                                                                         &
!              file= 'bx_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'bx_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',                                                                   &
              action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(2),                                                                         &
!              file= 'by_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'by_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
              action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(3),                                                                         &
!             file= 'bz_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'bz_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(4),                                                                         &
!             file='den_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'den_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(5),                                                                         &
!            file= 'ex_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'ex_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
              action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(6),                                                                         &
!             file= 'ey_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'ey_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
              action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(7),                                                                         &
!             file= 'ez_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'ez_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(8),                                                                         &
!             file='vix_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'vix_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(9),                                                                          &
!             file='viy_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'viy_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(10),                                                                         &
!             file='viz_'//trim(adjustl(cycle_ascii))//'.gda', &
              file= trim(trim(adjustl(data_directory))//'viz_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(11),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'tpar_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(12),                                                                           &
!             file='tperp_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'tperp_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        if (eta_par == 0) then
          open (file_unit(13),                                                                         &
!               file='eta_'//trim(adjustl(cycle_ascii))//'.gda',&
                file= trim(trim(adjustl(data_directory))//'eta_'//trim(adjustl(cycle_ascii)))//'.gda', &
                form='unformatted',   &
	  	action='write',access='direct', status='unknown',recl=lenrec)
        else
          open (file_unit(13),                                                                             &
!               file='eta_par_'//trim(adjustl(cycle_ascii))//'.gda',&
                file= trim(trim(adjustl(data_directory))//'eta_par_'//trim(adjustl(cycle_ascii)))//'.gda', &
                form='unformatted',   &
  	  	action='write',access='direct', status='unknown',recl=lenrec)
        endif
        open (file_unit(14),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'p_xx_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(15),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'p_xy_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(16),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'p_xz_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(17),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'p_yy_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(18),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'p_yz_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
        open (file_unit(19),                                                                          &
!             file='tpar_'//trim(adjustl(cycle_ascii))//'.gda',&
              file= trim(trim(adjustl(data_directory))//'p_zz_'//trim(adjustl(cycle_ascii)))//'.gda', &
              form='unformatted',   &
	      action='write',access='direct', status='unknown',recl=lenrec)
      endif
      return
      end subroutine openfiles
!
!***********************************************************************
!
      subroutine opendiagfiles
 
      use parameter_mod
      character timeunit*3, file_name*20
      integer*8 file_unit_ref
 
      file_unit_time = myid + 500
      write(timeunit,"(i3)") file_unit_time

!
! The diagnostic timing file is defined as a standard COS
! blocked file. Note that this file is ONLY used during
! diagnostic runs
! 
      if (notime == 0) then
        file_name = "timing" // timeunit // ".txt" 
        open(UNIT=file_unit_time,FILE=file_name,status='unknown')
      endif
      return
      end subroutine opendiagfiles
!
!***********************************************************************
!
      subroutine restrtrw(rw,itstart)
!
!=======================================================================
!
      use parameter_mod
      integer*8 f_unit,itstart,np_count
      real rw
      double precision, dimension(:), allocatable:: particle_tmp_array
 
 
      if (rw == +1.0) then
 
        t_stopped = t_stopped + (it-itstart+1)*dtwci
        f_unit=215+myid
!        open(unit=f_unit,file='restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        open(unit=f_unit,file=trim(adjustl(restart_directory))//'restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        ,form='unformatted',status='unknown')
      
!        Old restart format - dump all allocated particle memory
!        write(f_unit) x,y,z,vx,vy,vz,qp,link,porder

!       New restart format - dump only non-trivial particles
        write(f_unit) nspec
        DO IS=1,NSPEC
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
          write(f_unit) nptotp

          allocate (particle_tmp_array(nptotp))

!         x
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = x(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

!         y
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = y(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

!         z
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = z(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

!         vx
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = vx(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

!         vy
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = vy(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

!         vz
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = vz(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

!         qp
          NPTOTP=0
          do ize=kb-1,ke
            do iye=jb-1,je
              do ixe=1,nx1   
                  NP=IPHEAD(ixe,iye,ize,is)
                  DO WHILE (NP.NE.0)
                    NPTOTP=NPTOTP+1
                    particle_tmp_array(nptotp) = qp(np)
                    NP=LINK(NP)
                  ENDDO
              enddo
            enddo
          enddo
          write(f_unit) particle_tmp_array

          deallocate (particle_tmp_array)

        ENDDO

        write(f_unit) ninj,ninj_global,nescape,nescape_global,npart, &
        npart_global,qleft,qrite,maximum_simulation_time

        write(f_unit) x0,x1,tx0,vpar,vper,vbal,bbal,rcorr,teti,ishape

        write(f_unit) btspec, qspec, wspec, frac, vxbar, vybar,      &
        vzbar, anisot, denmin, resis, wpiwci, bete, fxsho,ave1,      &
        ave2,phib, xmax,ymax,zmax,bxc,byc,bzc,gama,                  &
        npx, npy, npz,                                               &
        iterb,norbskip,restrt_write,nxcel,netax,netay,netaz,nspec,   &
        nx,ny,nz,nskipx,nskipy,                                      &
        nskipz, testorbt, restart,etamin,etamax,ieta,eta_par,bz_IMF

        write(f_unit) hx,hy,hz,hxi,hyi,hzi                           &
       ,pi,efld,bfld,time,te0                                        &
       ,prntinfo,wrtdat,itfin,iwt                                    &
       ,nx1,nx2,ny1,ny2,nz1,nz2,it                                   &
!       ,ipstore,nptot,npleaving,npentering,myid_stop                 &
       ,nptot,npleaving,npentering,myid_stop                 &
       ,iclock_speed,iopen,iseed, file_unit,file_unit_read           &
       ,file_unit_time,notime,file_unit_tmp                          &
       ,clock_time_init,clock_time_old,clock_time                    &
       ,clock_time1

        write(f_unit) dfac,nskip,ipleft,iprite,ipsendleft,ipsendrite &
        ,iprecv,ipsendtop,ipsendbot,ipsendlefttop,ipsendleftbot      &
        ,ipsendritetop,ipsendritebot,ipsend

        write(f_unit) bx,by,bz,den,pe,eta,ex,ey,ez,fox,foy,foz       &
        ,bdipole_x,bdipole_y,bdipole_z                               &
        ,dipole_sphere_ex,dipole_sphere_ey,dipole_sphere_ez          &
        ,eta_times_b_dot_j

        write(f_unit) vix, viy, viz, vixo, viyo, vizo

!        write(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp
        write(f_unit) tpar,tperp,dns,vxs,vys,vzs

        close(unit=f_unit)


      else if (rw == -1.0) then
 
        f_unit=215+myid
!        open(unit=f_unit,file='restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        open(unit=f_unit,file=trim(adjustl(restart_directory))//'restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        ,form='unformatted',status='unknown')
      
!        Old restart format - dump all allocated particle memory
!        read(f_unit) x,y,z,vx,vy,vz,qp,link,porder

!       New restart format - dump only non-trivial particles
        read(f_unit) nspec
        DO IS=1,NSPEC
          read(f_unit) NPTOTP
          allocate (particle_tmp_array(nptotp))

!         x
          read(f_unit) particle_tmp_array
          ixe=2
          iye=jb
          ize=kb
          do np_count=1,NPTOTP
            NP=IPSTORE
            x(np)=particle_tmp_array(np_count)
            ipstore=link(np)
            link(np)=iphead(ixe,iye,ize,is)
            iphead(ixe,iye,ize,is)=np
          enddo
          iptemp(ixe,iye,ize,is)=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            iphead(ixe,iye,ize,is)=link(np)
            link(np)=iptemp(ixe,iye,ize,is)
            iptemp(ixe,iye,ize,is)=np
            np=iphead(ixe,iye,ize,is)
          ENDDO
          iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
          iptemp(ixe,iye,ize,is)=0

!         y
          read(f_unit) particle_tmp_array
          NP_COUNT=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            NP_COUNT=NP_COUNT+1
            y(np)=particle_tmp_array(np_count)
            np=link(np)
          ENDDO

!         z
          read(f_unit) particle_tmp_array
          NP_COUNT=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            NP_COUNT=NP_COUNT+1
            z(np)=particle_tmp_array(np_count)
            np=link(np)
          ENDDO

!         vx
          read(f_unit) particle_tmp_array
          NP_COUNT=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            NP_COUNT=NP_COUNT+1
            vx(np)=particle_tmp_array(np_count)
            np=link(np)
          ENDDO

!         vy
          read(f_unit) particle_tmp_array
          NP_COUNT=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            NP_COUNT=NP_COUNT+1
            vy(np)=particle_tmp_array(np_count)
            np=link(np)
          ENDDO
   
!         vz
          read(f_unit) particle_tmp_array
          NP_COUNT=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            NP_COUNT=NP_COUNT+1
            vz(np)=particle_tmp_array(np_count)
            np=link(np)
          ENDDO
   
!         qp
          read(f_unit) particle_tmp_array
          NP_COUNT=0
          np=iphead(ixe,iye,ize,is)
          DO WHILE (NP.NE.0)
            NP_COUNT=NP_COUNT+1
            qp(np)=particle_tmp_array(np_count)
            np=link(np)
          ENDDO

          deallocate (particle_tmp_array)

        ENDDO

        read(f_unit) ninj,ninj_global,nescape,nescape_global,npart,  &
        npart_global,qleft,qrite,t_stopped

        read(f_unit) x0,x1,tx0,vpar,vper,vbal,bbal,rcorr,teti,ishape

        read(f_unit) btspec, qspec, wspec, frac, vxbar, vybar,       &
        vzbar, anisot, denmin, resis, wpiwci, bete, fxsho,ave1,      &
        ave2,phib, xmax,ymax,zmax,bxc,byc,bzc,gama,                  &
        npx, npy, npz,                                               &
        iterb,norbskip,restrt_write,nxcel,netax,netay,netaz,nspec,   &
        nx,ny,nz,nskipx,nskipy,                                      &
        nskipz, testorbt, restart,etamin,etamax,ieta,eta_par,bz_IMF

        read(f_unit) hx,hy,hz,hxi,hyi,hzi                            &
       ,pi,efld,bfld,time,te0                                        &
       ,prntinfo,wrtdat,itfin,iwt                     &
       ,nx1,nx2,ny1,ny2,nz1,nz2,it                                   &
!       ,ipstore,nptot,npleaving,npentering,myid_stop                 &
       ,nptot,npleaving,npentering,myid_stop                 &
       ,iclock_speed,iopen,iseed, file_unit,file_unit_read           &
       ,file_unit_time,notime,file_unit_tmp                          &
       ,clock_time_init,clock_time_old,clock_time                    &
       ,clock_time1

        read(f_unit) dfac,nskip,ipleft,iprite,ipsendleft,ipsendrite  &
        ,iprecv,ipsendtop,ipsendbot,ipsendlefttop,ipsendleftbot      &
        ,ipsendritetop,ipsendritebot,ipsend

        read(f_unit) bx,by,bz,den,pe,eta,ex,ey,ez,fox,foy,foz        &
        ,bdipole_x,bdipole_y,bdipole_z                               &
        ,dipole_sphere_ex,dipole_sphere_ey,dipole_sphere_ez          &
        ,eta_times_b_dot_j

        read(f_unit) vix, viy, viz, vixo, viyo, vizo

!        read(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp
        read(f_unit) tpar,tperp,dns,vxs,vys,vzs

        close(unit=f_unit)

! Reset electic field and fox,y,z 

         noresete = 1
         if (noresete == 0) then
           deno=den;vixo=vix;viyo=viy;vizo=viz
           call ecalc( 1 )
           call focalc
         endif

      endif
    
        call sortit
!
!***********************************************************************
!
      return
      end
!
!***********************************************************************
!
      subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
 
      integer   numprocs,myid
      integer*8 n
      integer*8 s, e
      integer*8 nlocal
      integer*8 deficit
 
      nlocal  = n / numprocs
      s       = myid * nlocal + 1
      deficit = mod(n,numprocs)
      s       = s + min(myid,deficit)
      if (myid  <  deficit) then
          nlocal = nlocal + 1
      endif
      e = s + nlocal - 1
      if (e  >  n .or. myid  ==  numprocs-1) e = n
 
 
      return
      end
!
!***********************************************************************
!
      subroutine count(n,nskip,ncount)

      ncount=0
      do i=1,n,nskip
        ncount=ncount+1
      enddo

      return
      end
!
!***********************************************************************
!
      subroutine trans
 
      use parameter_mod
      use MESH2D
      integer*8:: is,i,j,k,jbmin,jbmax,kbmin,kbmax
      integer*8:: time_begin(8),time_end(8)
      double precision:: dttmp,dns_tmp
 
      call date_and_time(values=time_begin_array(:,20))
 
      do is=1,nspec
        DO K=KB-1,KE+1
          do j=jb-1,je+1
            do i=1,nx2
              dns(i,j,k,is)=1.e-10
              vxs(i,j,k,is)=0.
              vys(i,j,k,is)=0.
              vzs(i,j,k,is)=0.
              if (is == 1) then
                deno(i,j,k)=den(i,j,k)
                vixo(i,j,k)=vix(i,j,k)
                viyo(i,j,k)=viy(i,j,k)
                vizo(i,j,k)=viz(i,j,k)
                den(i,j,k)=0.
                vix(i,j,k)=0.
                viy(i,j,k)=0.
                viz(i,j,k)=0.
              endif
            enddo
          enddo
        enddo
      enddo
 
      call date_and_time(values=time_end)
      clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 
      if (notime == 0) then
        write(file_unit_time,"(i4,' parmovin ',e20.8)") it,real(clock_time-clock_time_init)
      endif
      call date_and_time(values=time_begin_array(:,7))

  
      if (ndim /= 1) then
         call parmov
       else
         call parmov_2d
       endif

      call date_and_time(values=time_end_array(:,7))
      call accumulate_time_difference(time_begin_array(1,7),time_end_array(1,7),time_elapsed(7))
 
 
      call date_and_time(values=time_end)
      clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 
      if (notime == 0) then
        write(file_unit_time,"(i4,' prmovout ',e20.8)") it,real(clock_time-clock_time_init)
      endif

      if (testorbt) return

      do is=1,nspec
        do k=kb-1,ke+1
          do j=jb-1,je+1
            do i=1,nx2

!             Nonuniform mesh
              cell_volume_ratio = hx*hy*hz/(meshX%dxc(i)*meshY%dxc(j+1)*meshZ%dxc(k+1))

              dns_tmp=dns(i,j,k,is)

!             Uniform mesh - Same as is in version 5.0
!              if (dns_tmp <= denmin) dns_tmp=1.d+10

!             Nonuniform mesh
!              if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=1.d+10
              if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=denmin/cell_volume_ratio ! July 21, 2010

              vxs(i,j,k,is)=vxs(i,j,k,is)/dns_tmp
              vys(i,j,k,is)=vys(i,j,k,is)/dns_tmp
              vzs(i,j,k,is)=vzs(i,j,k,is)/dns_tmp

!             Uniform mesh - Same as is in version 5.0
!              dns(i,j,k,is)=dns(i,j,k,is)

!             Nonuniform mesh
!              dns(i,j,k,is)=dns(i,j,k,is)*cell_volume_ratio

            enddo
          enddo
        enddo
      enddo
 
      do is=1,nspec
        do k=kb-1,ke+1
          do j=jb-1,je+1
            do i=1,nx2 
            den(i,j,k)=den(i,j,k)+dns(i,j,k,is)*qspec(is) 
            vix(i,j,k)=vix(i,j,k)+qspec(is)*dns(i,j,k,is)*vxs(i,j,k,is) 
            viy(i,j,k)=viy(i,j,k)+qspec(is)*dns(i,j,k,is)*vys(i,j,k,is) 
            viz(i,j,k)=viz(i,j,k)+qspec(is)*dns(i,j,k,is)*vzs(i,j,k,is)
            enddo
          enddo
        enddo
      enddo

      do k=kb-1,ke+1
          do j=jb-1,je+1
          do i=1,nx2
            den(i,j,k)=max(denmin,den(i,j,k))
            pe(i,j,k) =te0*den(i,j,k)**gama
            vix(i,j,k)=vix(i,j,k)/den(i,j,k)
            viy(i,j,k)=viy(i,j,k)/den(i,j,k)
            viz(i,j,k)=viz(i,j,k)/den(i,j,k)
          enddo
        enddo
      enddo
 
!     Nonuniform mesh
      if (ndim /= 1) then
      call XREALBCC(PE ,1,NX,NY,NZ)
      call XREALBCC(DEN,1,NX,NY,NZ)
      call XREALBCC(VIX,1,NX,NY,NZ)
      call XREALBCC(VIY,1,NX,NY,NZ)
      call XREALBCC(VIZ,1,NX,NY,NZ)
      else
      call XREALBCC_2D(PE ,1,NX,NY,NZ)
      call XREALBCC_2D(DEN,1,NX,NY,NZ)
      call XREALBCC_2D(VIX,1,NX,NY,NZ)
      call XREALBCC_2D(VIY,1,NX,NY,NZ)
      call XREALBCC_2D(VIZ,1,NX,NY,NZ)
      endif

!     Smoothing
      if (smoothing .and. smooth_coef == 0.) then
        if (ndim /=1) then
          call nsmth(PE )
          call nsmth(DEN)
          call nsmth(VIX)
          call nsmth(VIY)
          call nsmth(VIZ)
          call XREALBCC(PE ,1,NX,NY,NZ)
          call XREALBCC(DEN,1,NX,NY,NZ)
          call XREALBCC(VIX,1,NX,NY,NZ)
          call XREALBCC(VIY,1,NX,NY,NZ)
          call XREALBCC(VIZ,1,NX,NY,NZ)
        else
          call nsmth_2d(PE ,NX2,NY2,NZ2)
          call nsmth_2d(DEN,NX2,NY2,NZ2)
          call nsmth_2d(VIX,NX2,NY2,NZ2)
          call nsmth_2d(VIY,NX2,NY2,NZ2)
          call nsmth_2d(VIZ,NX2,NY2,NZ2)
          call XREALBCC_2D(PE ,1,NX,NY,NZ)
          call XREALBCC_2D(DEN,1,NX,NY,NZ)
          call XREALBCC_2D(VIX,1,NX,NY,NZ)
          call XREALBCC_2D(VIY,1,NX,NY,NZ)
          call XREALBCC_2D(VIZ,1,NX,NY,NZ)
        endif
      endif

 
      if (it == 0) then
         deno=den;vixo=vix;viyo=viy;vizo=viz
      endif
 
      call date_and_time(values=time_begin_array(:,8))
 
      call energy
 
      call date_and_time(values=time_end_array(:,8))
      call accumulate_time_difference(time_begin_array(1,8),time_end_array(1,8),time_elapsed(8))
 
 
       kbmin = kb-1
       kbmax = ke+1
 
 
       jbmin = jb-1
       jbmax = je+1
 
 
      call date_and_time(values=time_end_array(:,20))
      call accumulate_time_difference(time_begin_array(1,20),time_end_array(1,20),time_elapsed(20))
 
      return
      end
!
!#######################################################################
!
      subroutine xreal(a,nx1m,ny1m,nz1m)
!
!  XREAL exchanges the buffer (ghost) cell information in the planes 
!        kb-1 and ke+1.
!  There are up to four exchanges that take place for each processor:
!  1 & 2)  Processors myid and nbrleft exchange contributions in cells kb-1 
!          (myid) and ke+1 (nbrleft)
!  3 & 4)  Processors myid and nbrrite exchange contributions in cells ke+1 
!          (myid) and kb-1 (nbrleft)
!
!  Fewer exchanges occur if myid is either at the left or right boundary of
!  the physical domain
!
!  Blocking is avoided by first sending to the right and then to the left
!

      use parameter_mod
      integer*8 i,j,nx1m,ny1m,nz1m,k
      double precision a(nxmax,jb-1:je+1,kb-1:ke+1),tmp(nxmax,jb-1:je+1,kb-1:ke+1)
      a(2   ,:,:)=a(2   ,:,:)+a(1   ,:,:)
      a(nx1m+1,:,:)=a(nx1m+1,:,:)+a(nx1m+2,:,:)
      call MPI_SENDRECV(a    (1    ,je+1,kb-1),1,stridery,nbrrite,0, &
                        tmp  (1    ,jb-1,kb-1),1,stridery,nbrleft,0, &
                        mpi_comm_world,status,ierr)
      if (jb == 1) tmp(:,jb-1,:)=a(:,jb-1,:)
      a(:,jb,:)=a(:,jb,:)+tmp  (:,jb-1,:)
      call MPI_SENDRECV(a    (1    ,jb-1,kb-1),1,stridery,nbrleft,1, &
                        tmp  (1    ,je+1,kb-1),1,stridery,nbrrite,1, &
                        mpi_comm_world,status,ierr)
      if (je == ny1m) tmp(:,je+1,:)=a(:,je+1,:)
      a(:,je,:)=a(:,je,:)+tmp  (:,je+1,:)
      call MPI_SENDRECV(a    (1    ,jb-1,ke+1),1,striderz,nbrtop ,0, &
                        tmp  (1    ,jb-1,kb-1),1,striderz,nbrbot ,0, &
                        mpi_comm_world,status,ierr)
      if (kb == 1) tmp(:,:,kb-1)=a(:,:,kb-1)
      a(:,:,kb)=a(:,:,kb)+tmp  (:,:,kb-1)
      call MPI_SENDRECV(a    (1    ,jb-1,kb-1),1,striderz,nbrbot ,1, &
                        tmp  (1    ,jb-1,ke+1),1,striderz,nbrtop ,1, &
                        mpi_comm_world,status,ierr)
      if (ke == nz1m) tmp(:,:,ke+1)=a(:,:,ke+1)
      a(:,:,ke)=a(:,:,ke)+tmp  (:,:,ke+1)
 
 
      a(1   ,:,:)=a(2   ,:,:)
      a(nx1m+2,:,:)=a(nx1m+1,:,:)
 
 
      return
      end
!
!#######################################################################
!
      subroutine xrealbcc(a, ibnd, nx1m, ny1m,nz1m)

!  XREALBCC updates the ghost cells in the planes kb-1 and ke+1 by obtaining
!  latest values for neighboring processors (nbrleft and nbrrite)
!
!  There are up to four exchanges that take place for each processor:
!  1) cell ke+1 in nbrleft is replaced by values in cell kb  in myid
!  2) cell kb-1 in nbrrite is replaced by values in cell ke  in myid
!  3) cell ke+1 in myid is replaced by values in cell kb  in nbrrite
!  4) cell kb-1 in myid is replaced by values in cell ke  in nbrleft
!
!  Fewer exchanges occur if myid is either at the left or right boundary of
!  the physical domain
!
!  Blocking is avoided by first sending to the right and then to the left
!
!  If ibnd = 1, then zero slope boundary conditions for k = 0 and nz1 are 
!             set at end of routine
!
      use parameter_mod
      integer*8 ibnd,i,j,nx1m,ny1m,nz1m
      double precision a(nxmax,jb-1:je+1,kb-1:ke+1),tmp(nxmax,jb-1:je+1,kb-1:ke+1)
 
      tmp=a
      call MPI_SENDRECV(a(1    ,jb-1,ke  ),1,striderz,nbrtop ,0,&
                        a(1    ,jb-1,kb-1),1,striderz,nbrbot ,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(a(1    ,jb-1,kb  ),1,striderz,nbrbot ,1,&
                        a(1    ,jb-1,ke+1),1,striderz,nbrtop ,1,&
                        mpi_comm_world,status,ierr)
      if (kb == 1) then
        if (ibnd == 1) then
          a(:,:,kb-1)=a(:,:,kb)
        else
          a(:,:,kb-1)=tmp(:,:,kb-1)
        endif
      endif
      if (ke == nz1m) then
        if (ibnd == 1) then
          a(:,:,ke+1)=a(:,:,ke)
        else
          a(:,:,ke+1)=tmp(:,:,ke+1)
        endif
      endif
      call MPI_SENDRECV(a(1    ,je  ,kb-1),1,stridery,nbrrite,0,&
                        a(1    ,jb-1,kb-1),1,stridery,nbrleft,0,&
                        mpi_comm_world,status,ierr)
      call MPI_SENDRECV(a(1    ,jb  ,kb-1),1,stridery,nbrleft,1,&
                        a(1    ,je+1,kb-1),1,stridery,nbrrite,1,&
                        mpi_comm_world,status,ierr)
      if (jb == 1) then
        if (ibnd == 1)  then
          a(:,jb-1,:)=a(:,jb,:)
        else
          a(:,jb-1,:)=tmp(:,jb-1,:)
        endif
      endif
      if (je == ny1m) then
        if (ibnd == 1)  then
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
      subroutine energy
 
      use parameter_mod
 
      efldp=0.
      bfldp=0.
      do k=kb,ke
        do j=jb,je
          do i=2,nx1
            efldp=efldp+ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2
            bfldp=bfldp+bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2
          enddo
        enddo
      enddo
      efldp=efldp*hx*hy*hz*0.5
      bfldp=bfldp*hx*hy*hz*0.5
      call MPI_ALLREDUCE(efldp,efld,1,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(bfldp,bfld,1,MPI_DOUBLE_PRECISION,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      return
      end
!
!#######################################################################
!
      subroutine parmov
 
      use parameter_mod
      use MESH2D
      integer*8 count_kbq,time_begin(8),time_end(8)
      integer*8 nptotp_kbq,npart_kbq(2),np_ijk,Storage_Error_p,Storage_Error
      data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
      data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
      data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
      integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                 ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
      double precision:: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,xpart,ypart,zpart,r_particle
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
      double precision, dimension(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bx_av,by_av,bz_av
      double precision:: v_limit,eps2,rx0,ry0,rz0,rrat,sqrr,outer_radius,myranf,twopi,fluxran,vxa,vyz,vza
      INTEGER*8:: L, EXIT_CODE_P, EXIT_CODE
      integer*8:: n_fast_removed,n_fast_removed_local,nptot_max,Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p
      double precision:: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min,x_disp,y_disp,z_disp          &
                        ,y_disp_max_p,x_disp_max_p,z_disp_max_p,y_disp_max,x_disp_max,z_disp_max
      double precision:: disp_max_p(3),disp_max(3),tx,ty,tz,v_x,v_y,v_z  
      INTEGER*4 :: nescapearr(8),nescapearr_global(8)
      INTEGER*4 :: ppacket(3),ppacketg(3),dpacket(4),dpacketg(4)
      INTEGER*8 :: epacket(2),epacketg(2),indx,loop
      INTEGER*8,dimension(:),allocatable :: nparr
      double precision, dimension(3,nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bxyz_av
      double precision:: TEX1,TEX2,TEX3,TEX4,TEX5,TEX6,TEX7,TEX8  
      double precision:: TEY1,TEY2,TEY3,TEY4,TEY5,TEY6,TEY7,TEY8  
      double precision:: TEZ1,TEZ2,TEZ3,TEZ4,TEZ5,TEZ6,TEZ7,TEZ8  
      double precision:: mX_xa,mX_ta,mX_ca1,mX_ca2,mX_xb,mX_dtdx,mX_tb,mX_cb1,mX_cb2
      double precision:: mY_xa,mY_ta,mY_ca1,mY_ca2,mY_xb,mY_dtdx,mY_tb,mY_cb1,mY_cb2
      double precision:: mZ_xa,mZ_ta,mZ_ca1,mZ_ca2,mZ_xb,mZ_dtdx,mZ_tb,mZ_cb1,mZ_cb2
      allocate (nparr(nplmax))
 
      call date_and_time(values=time_begin_array(:,19))

      Storage_Error_p = 0
      Field_Diverge_p = 0
      Courant_Violation_p = 0
      x_disp_max_p        = 0
      y_disp_max_p        = 0
      z_disp_max_p        = 0

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

 
      epsilon= buffer_zone
      d_ranf=1./1001.
      twopi=2.*acos(-1.)

      outer_radius=dipole_sphere_radius+min(hx,hy,hz)
      eps2=1.d-25
 
 
!     Uniform mesh - Same as is in version 5.0
!      x_dipole=hx*(i_dipole+0.5)
!      y_dipole=hy*(j_dipole+0.5)
!      z_dipole=hz*(k_dipole+0.5)
 
!      Nonuniform mesh
       x_dipole=meshX%xc(i_dipole)
       y_dipole=meshY%xc(j_dipole)
       z_dipole=meshZ%xc(k_dipole)


!     Uniform mesh - Same as in version 5.0
!      if (dt /=0.) then
!        if (ndim == 1) then
!         v_limit = min(hx/dt,hy/dt)
!        else
!         v_limit = min(hx/dt,hy/dt,hz/dt)
!      endif
!      else
!        v_limit=1.d+10
!      endif
!

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
 
      bx_av=0.;by_av=0.;bz_av=0.
      DO K = KB-1,KE
        DO J = JB-1,JE
          DO I = 1, NX1
            bx_av(i,j,k)=0.125*( bx(i  ,j  ,k  )+bdipole_x(i  ,j  ,k  )             &
                                +bx(i+1,j  ,k  )+bdipole_x(i+1,j  ,k  )             &
                                +bx(i  ,j+1,k  )+bdipole_x(i  ,j+1,k  )             &
                                +bx(i+1,j+1,k  )+bdipole_x(i+1,j+1,k  )             &
                                +bx(i  ,j  ,k+1)+bdipole_x(i  ,j  ,k+1)             &
                                +bx(i+1,j  ,k+1)+bdipole_x(i+1,j  ,k+1)             &
                                +bx(i  ,j+1,k+1)+bdipole_x(i  ,j+1,k+1)             &
                                +bx(i+1,j+1,k+1)+bdipole_x(i+1,j+1,k+1)             &
                               )
            by_av(i,j,k)=0.125*( by(i  ,j  ,k  )+bdipole_y(i  ,j  ,k  )             &
                                +by(i+1,j  ,k  )+bdipole_y(i+1,j  ,k  )             &
                                +by(i  ,j+1,k  )+bdipole_y(i  ,j+1,k  )             &
                                +by(i+1,j+1,k  )+bdipole_y(i+1,j+1,k  )             &
                                +by(i  ,j  ,k+1)+bdipole_y(i  ,j  ,k+1)             &
                                +by(i+1,j  ,k+1)+bdipole_y(i+1,j  ,k+1)             &
                                +by(i  ,j+1,k+1)+bdipole_y(i  ,j+1,k+1)             &
                                +by(i+1,j+1,k+1)+bdipole_y(i+1,j+1,k+1)             &
                               )
            bz_av(i,j,k)=0.125*( bz(i  ,j  ,k  )+bdipole_z(i  ,j  ,k  )             &
                                +bz(i+1,j  ,k  )+bdipole_z(i+1,j  ,k  )             &
                                +bz(i  ,j+1,k  )+bdipole_z(i  ,j+1,k  )             &
                                +bz(i+1,j+1,k  )+bdipole_z(i+1,j+1,k  )             &
                                +bz(i  ,j  ,k+1)+bdipole_z(i  ,j  ,k+1)             &
                                +bz(i+1,j  ,k+1)+bdipole_z(i+1,j  ,k+1)             &
                                +bz(i  ,j+1,k+1)+bdipole_z(i  ,j+1,k+1)             &
                                +bz(i+1,j+1,k+1)+bdipole_z(i+1,j+1,k+1)             &
                               )
          enddo
        enddo
      enddo
      CALL XREALBCC_PACK_B(BX_AV,BY_AV,BZ_AV,1,NX,NY,NZ)

      if (myid == 0) WRITE(6,*) " CALLING PARMOVE, NSPEC = ",NSPEC
      if (myid == 0) WRITE(6,*) " ABSORBING_DIPOLE       = ",absorbing_dipole
      if (myid == 0) WRITE(6,*) " DIPOLE_SPHERE_RADIUS   = ",dipole_sphere_radius
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
 
!  beginning of main particle loop
 
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



!        Diagnostic print statement
!
!        call MPI_ALLREDUCE(nptotp,nptot_max,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,IERR)
!        if (nptotp  == nptot_max) then
!          write(6,*) " PROCESSOR # ",MYID,", # OF PARTICLES = ",NPTOTP
!          write(6,*) " MAXIMUM # OF PARTICLES ALLOWED      = ",NPLMAX
!        endif

        call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
        IF (MYID.EQ.0) THEN
          WRITE(6,*) " IS = ",IS
          WRITE(6,*) " TOTAL # OF PARTICLES BEFORE PARMOV = ",NPTOT
        ENDIF

        wmult=wspec(is)
        h=dt*qspec(is)/wmult
        hh=.5*h
 
 
        call date_and_time(values=time_end)
        clock_time1=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        IF (DT.NE.0) THEN
          npart(is) = 0
          npart_kbq(is) = 0
          DO IIZE = KB-1,KE
            DO IIYE = JB-1,JE
              DO IIXE = 1, NX1
                NP=IPHEAD(IIXE,IIYE,IIZE,IS)
 
!  begin advance of particle position and velocity
!  If dt=0, skip
 
                DO WHILE (NP.NE.0)
                  L=NP
                  npart_kbq(is) = npart_kbq(is)+1
                  nptotp_kbq = nptotp_kbq + 1

!                 Uniform mesh - Same as in version 5.0
!                  rxe=hxi*x(l)+1.5000000000000001d+00
!                  rye=hyi*y(l)+0.5000000000000001d+00
!                  rze=hzi*z(l)+0.5000000000000001d+00
!                  ixe=rxe
!                  iye=rye
!                  ize=rze
!                  ixe=max(1   ,min(ixe,nx1))
!                  iye=max(jb-1,min(iye,je))
!                  ize=max(kb-1,min(ize,ke))
!                  ixep1 = ixe+1
!                  iyep1 = iye+1
!                  izep1 = ize+1
!                  fxe=rxe-ixe
!                  fye=rye-iye
!                  fze=rze-ize

!                 Nonuniform mesh - without using MESH_UNMAP
!                  rxe=hxi*x(l)+1.500000000000000d+00
!                  rye=hyi*y(l)+1.500000000000000d+00
!                  rze=hzi*z(l)+1.500000000000000d+00
!                  ixe=rxe
!                  iye=rye
!                  ize=rze
!                  ixe=ixc_2_c_map(ixe)
!                  iye=iyc_2_c_map(iye)
!                  ize=izc_2_c_map(ize)
!                  fxe=(x(l)-meshX%xc(ixe))/meshX%dxn(ixep1)
!                  fye=(y(l)-meshY%xc(iye))/meshY%dxn(iyep1)
!                  fze=(z(l)-meshZ%xc(ize))/meshZ%dxn(izep1)

!                 Nonuniform mesh - using MESH_UNMAP
                  rxe=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000d+00
                  rye=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000d+00
                  rze=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000d+00
                  ixe=rxe
                  iye=rye
                  ize=rze
                  fxe=rxe-ixe
                  fye=rye-iye
                  fze=rze-ize
                  iye=iye-1             ! integer index in y direction starts at 0
                  ize=ize-1             ! integer index in z direction starts at 0
 
                  ixep1 = ixe+1
                  iyep1 = iye+1
                  izep1 = ize+1

!                 Diagnosic test on particle cell index algorithm
!                  if (     fxe < 0. .or. fxe > 1.                &
!                      .or. fye < 0. .or. fye > 1.                &
!                      .or. fze < 0. .or. fze > 1.) then 
!                      write(6,*) " SCATTER LOOP"
!                      write(6,*) fxe,fye,fze
!                      write(6,*) " x;",x(l),meshX%xn(ixe),meshX%xc(ixe)
!                      write(6,*) " y;",y(l),meshY%xn(iye),meshY%xc(iye)
!                      write(6,*) " z;",z(l),meshZ%xn(ize),meshZ%xc(ize)
!                      write(6,*) " r; ",rxe,rye,rze
!                      write(6,*) " ixc_map; ",ixe,iye,ize
!                      write(6,*) " xc     ; ",meshX%xc(ixe),meshY%xc(iye),meshZ%xc(ize)
!                  endif
 

                  w1e=(1.-fxe)*(1.-fye)*(1.-fze)
                  w2e=fxe*(1.-fye)*(1.-fze)
                  w3e=(1.-fxe)*fye*(1.-fze)
                  w4e=fxe*fye*(1.-fze)
                  w5e=(1.-fxe)*(1.-fye)*fze
                  w6e=fxe*(1.-fye)*fze
                  w7e=(1.-fxe)*fye*fze
                  w8e=fxe*fye*fze
 
                  ex1=ex(ixe  ,iye  ,ize  )
                  ex2=ex(ixep1,iye  ,ize  )
                  ex3=ex(ixe  ,iyep1,ize  )
                  ex4=ex(ixep1,iyep1,ize  )
                  ex5=ex(ixe  ,iye  ,izep1)
                  ex6=ex(ixep1,iye  ,izep1)
                  ex7=ex(ixe  ,iyep1,izep1)
                  ex8=ex(ixep1,iyep1,izep1)
                  ey1=ey(ixe  ,iye  ,ize  )
                  ey2=ey(ixep1,iye  ,ize  )
                  ey3=ey(ixe  ,iyep1,ize  )
                  ey4=ey(ixep1,iyep1,ize  )
                  ey5=ey(ixe  ,iye  ,izep1)
                  ey6=ey(ixep1,iye  ,izep1)
                  ey7=ey(ixe  ,iyep1,izep1)
                  ey8=ey(ixep1,iyep1,izep1)
                  ez1=ez(ixe  ,iye  ,ize  )
                  ez2=ez(ixep1,iye  ,ize  )
                  ez3=ez(ixe  ,iyep1,ize  )
                  ez4=ez(ixep1,iyep1,ize  )
                  ez5=ez(ixe  ,iye  ,izep1)
                  ez6=ez(ixep1,iye  ,izep1)
                  ez7=ez(ixe  ,iyep1,izep1)
                  ez8=ez(ixep1,iyep1,izep1)
 
                  bx1=bx_av(ixe  ,iye  ,ize  )
                  bx2=bx_av(ixep1,iye  ,ize  )
                  bx3=bx_av(ixe  ,iyep1,ize  )
                  bx4=bx_av(ixep1,iyep1,ize  )
                  bx5=bx_av(ixe  ,iye  ,izep1)
                  bx6=bx_av(ixep1,iye  ,izep1)
                  bx7=bx_av(ixe  ,iyep1,izep1)
                  bx8=bx_av(ixep1,iyep1,izep1)
                  by1=by_av(ixe  ,iye  ,ize  )
                  by2=by_av(ixep1,iye  ,ize  )
                  by3=by_av(ixe  ,iyep1,ize  )
                  by4=by_av(ixep1,iyep1,ize  )
                  by5=by_av(ixe  ,iye  ,izep1)
                  by6=by_av(ixep1,iye  ,izep1)
                  by7=by_av(ixe  ,iyep1,izep1)
                  by8=by_av(ixep1,iyep1,izep1)
                  bz1=bz_av(ixe  ,iye  ,ize  )
                  bz2=bz_av(ixep1,iye  ,ize  )
                  bz3=bz_av(ixe  ,iyep1,ize  )
                  bz4=bz_av(ixep1,iyep1,ize  )
                  bz5=bz_av(ixe  ,iye  ,izep1)
                  bz6=bz_av(ixep1,iye  ,izep1)
                  bz7=bz_av(ixe  ,iyep1,izep1)
                  bz8=bz_av(ixep1,iyep1,izep1)
 
                  fox1=fox(ixe  ,iye  ,ize  )
                  fox2=fox(ixep1,iye  ,ize  )
                  fox3=fox(ixe  ,iyep1,ize  )
                  fox4=fox(ixep1,iyep1,ize  )
                  fox5=fox(ixe  ,iye  ,izep1)
                  fox6=fox(ixep1,iye  ,izep1)
                  fox7=fox(ixe  ,iyep1,izep1)
                  fox8=fox(ixep1,iyep1,izep1)
                  foy1=foy(ixe  ,iye  ,ize  )
                  foy2=foy(ixep1,iye  ,ize  )
                  foy3=foy(ixe  ,iyep1,ize  )
                  foy4=foy(ixep1,iyep1,ize  )
                  foy5=foy(ixe  ,iye  ,izep1)
                  foy6=foy(ixep1,iye  ,izep1)
                  foy7=foy(ixe  ,iyep1,izep1)
                  foy8=foy(ixep1,iyep1,izep1)
                  foz1=foz(ixe  ,iye  ,ize  )
                  foz2=foz(ixep1,iye  ,ize  )
                  foz3=foz(ixe  ,iyep1,ize  )
                  foz4=foz(ixep1,iyep1,ize  )
                  foz5=foz(ixe  ,iye  ,izep1)
                  foz6=foz(ixep1,iye  ,izep1)
                  foz7=foz(ixe  ,iyep1,izep1)
                  foz8=foz(ixep1,iyep1,izep1)

                  exa=w1e*ex1+w2e*ex2+w3e*ex3+w4e*ex4      &
                     +w5e*ex5+w6e*ex6+w7e*ex7+w8e*ex8      &
                     +w1e*fox1+w2e*fox2+w3e*fox3+w4e*fox4  &
                     +w5e*fox5+w6e*fox6+w7e*fox7+w8e*fox8
                  eya=w1e*ey1+w2e*ey2+w3e*ey3+w4e*ey4      &
                     +w5e*ey5+w6e*ey6+w7e*ey7+w8e*ey8      &
                     +w1e*foy1+w2e*foy2+w3e*foy3+w4e*foy4  &
                     +w5e*foy5+w6e*foy6+w7e*foy7+w8e*foy8
                  eza=w1e*ez1+w2e*ez2+w3e*ez3+w4e*ez4      &
                     +w5e*ez5+w6e*ez6+w7e*ez7+w8e*ez8      &
                     +w1e*foz1+w2e*foz2+w3e*foz3+w4e*foz4  &
                     +w5e*foz5+w6e*foz6+w7e*foz7+w8e*foz8
 
                  bxa=w1e*bx1+w2e*bx2+w3e*bx3+w4e*bx4      &
                     +w5e*bx5+w6e*bx6+w7e*bx7+w8e*bx8
                  bya=w1e*by1+w2e*by2+w3e*by3+w4e*by4      &
                     +w5e*by5+w6e*by6+w7e*by7+w8e*by8
                  bza=w1e*bz1+w2e*bz2+w3e*bz3+w4e*bz4      &
                     +w5e*bz5+w6e*bz6+w7e*bz7+w8e*bz8
 
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
                  z_disp = dt*vz(l)

                  x(l)=x(l)+ x_disp
                  y(l)=y(l)+ y_disp
                  z(l)=z(l)+ z_disp

                  if ( abs(x_disp/meshX%dxn(ixep1)) > 1.0 .or.                                               &
                       abs(y_disp/meshY%dxn(iyep1)) > 1.0 .or.                                               &
                       abs(z_disp/meshZ%dxn(izep1)) > 1.0) Courant_Violation_p = Courant_Violation_p + 1
                  x_disp_max_p = max(x_disp_max_p,abs(x_disp)/meshX%dxn(ixep1))
                  y_disp_max_p = max(y_disp_max_p,abs(y_disp)/meshY%dxn(iyep1))
                  z_disp_max_p = max(z_disp_max_p,abs(z_disp)/meshZ%dxn(izep1))

                  NP=LINK(NP)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        disp_max_p(1) = x_disp_max_p
        disp_max_p(2) = y_disp_max_p
        disp_max_p(3) = z_disp_max_p
        call MPI_ALLREDUCE(disp_max_p,disp_max,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)   !LAURA
        x_disp_max = disp_max(1)
        y_disp_max = disp_max(2)
        z_disp_max = disp_max(3)

        if (myid == 0) then
          write(6,*) " maximum x-displacement/dx = ",x_disp_max
          write(6,*) " maximum y-displacement/dy = ",y_disp_max
          write(6,*) " maximum z-displacement/dz = ",z_disp_max
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
        nsendp=0
        nrecvp=0
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

        n_fast_removed_local = 0 

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
!
!     Remove fast particles
!
                if (abs(vx(l)) >= v_limit .or. abs(vy(l)) >= v_limit .or. abs(vz(l)) >= v_limit) then
                  n_fast_removed_local = n_fast_removed_local + 1
                  iphead(ixe,iye,ize,is)=link(np)
                  link(np)=ipstore
                  ipstore=np
                  np=iphead(ixe,iye,ize,is)
                  goto 15
                endif
!
!

                if (ndim /= 1) then
                  r_particle=sqrt((xpart-x_dipole)**2+(ypart-y_dipole)**2+(zpart-z_dipole)**2)
                else
                  r_particle=sqrt((xpart-x_dipole)**2+(ypart-y_dipole)**2)
                endif
                rx0=xpart-x_dipole
                ry0=ypart-y_dipole
                rz0=zpart-z_dipole
                sqrr=r_particle+eps2

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
                      z(l) = rz0 *rrat + z_dipole
                      vx(np)=-vx(np)
                      vy(np)=-vy(np)
                      vz(np)=-vz(np)
                      x (np)=x(np)+vx(np)*dt
                      y (np)=y(np)+vy(np)*dt
                      z (np)=z(np)+vz(np)*dt
                      iphead(ixe,iye,ize,is)=link(np)
                      link(np)=iptemp(ixe,iye,ize,is)
                      iptemp(ixe,iye,ize,is)=np
                      np=iphead(ixe,iye,ize,is)
                    endif
                    goto 15
                  endif

                  if       (y(l) <   epsilon     ) then
                    
!                    y(l) = 2.*epsilon-y(l)
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

!                    y(l) = 2.*(ymax-epsilon) - y(l)
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

                  if       (z(l) <   epsilon     ) then
                    
!                    z(l) = 2.*epsilon - z(l)
                    z(l) = myranf()*epsilon

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vz_tmp=vpar(is)*u_array_injection(iv,jv)
                    vz(l)=+vz_tmp
                    vmag=sqrt(-log(1.-.999999*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*cos(th)
                    vy(l)=          vmag*vper(is)*sin(th)
                    x (l)=   myranf()*xmax
                    y (l)=yb+myranf()*(ye-yb)

                  endif

                  if       (z(l) >   zmax-epsilon   ) then

!                    z(l) = 2.*(zmax-epsilon) - z(l)
                    z(l) = zmax-myranf()*epsilon

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vz_tmp=vpar(is)*u_array_injection(iv,jv)
                    vz(l)=-vz_tmp
                    vmag=sqrt(-log(1.-.999999*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*cos(th)
                    vy(l)=          vmag*vper(is)*sin(th)
                    x (l)=     myranf()*xmax
                    y (l)=yb  +myranf()*(ye-yb)

                  endif

                  xpart = x(l)
                  ypart = y(l)
                  zpart = z(l)

 
                if (    xpart < 0..or.xpart > xmax&
                    .or.ypart < 0..or.ypart > ymax&       
                    .or.zpart < 0..or.zpart > zmax) then
 
                  if (xpart < 0.) then
                    nescape_yz(is)=nescape_yz(is)+1
                  else if (xpart > xmax) then
                    nescape_zy(is)=nescape_zy(is)+1
                  else if (ypart < 0.) then
                    nescape_xz(is)=nescape_xz(is)+1
                  else if (ypart > ymax) then
                    nescape_zx(is)=nescape_zx(is)+1
                  else if (zpart < 0.) then
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
 
                  if ((zpart <= ze.and.zpart >= zb).and.(ypart <= ye.and.ypart >= yb)) then

                    iphead(ixe,iye,ize,is)=link(np)
                    link(np)=iptemp(ixe,iye,ize,is)
                    iptemp(ixe,iye,ize,is)=np
                    np=iphead(ixe,iye,ize,is)

                  ELSE

                    if (ypart <= ye.and.ypart >= yb) then
                      iye_cc=jb
                    else
                      if (ypart > ye) then
                        iye_cc=je+1 
                      else
                        iye_cc=jb-1 
                      endif
                    endif

                    if (zpart <= ze.and.zpart >= zb) then
                      ize_cc=kb
                    else
                      if (zpart > ze) then
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
 15             CONTINUE
              ENDDO
              iphead(ixe,iye,ize,is)=iptemp(ixe,iye,ize,is)
              iptemp(ixe,iye,ize,is)=0
            ENDDO
          ENDDO
        ENDDO
!
!************************************************************************
!
!     exchange data among processes and compute to see how many
!     particles each process has to send to, and receive from, other
!     processes
!
!     nsendp(nbrleft) -> how many particles to be sent from myid to 
!                         left neighbor
!     nsendp(nbrrite) -> how many particles to be sent from myid to 
!                         right neighbor
!     nrecvp(nbrleft) -> how many particles to be sent to myid from 
!                         left neighbor
!     nrecvp(nbrrite) -> how many particles to be sent to myid from 
!                         right neighbor
!
!
!  exchange information about particle numbers in two steps. First, 
!  send to right and receive from left. Note that the processors on 
!  the physical boundaries only send (myid == nbrleft) or receive 
!  (myid == nbrrite)
!
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

        ppacket(1) = nsendtotp
        ppacket(2) = nrecvtotp
        ppacket(3) = n_fast_removed_local
        call MPI_ALLREDUCE(ppacket,ppacketg,3,MPI_INTEGER4,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
        nsendtot       = ppacketg(1)
        nrecvtot       = ppacketg(2)
        n_fast_removed = ppacketg(3)

        if (myid == 0) then
          write(6,*) " FINISHED COMPILING LISTS "
          write(6,*) " # OF PARTICLES TO BE SENT     = ",NSENDTOT
          write(6,*) " # OF PARTICLES TO BE RECEIVED = ",NRECVTOT
          write(6,*) " # OF PARTICLES REMOVED BECAUSE V > VLIMIT = ",n_fast_removed
        endif
        if (NSENDTOT.NE.NRECVTOT) THEN
          CALL MPI_FINALIZE(IERR)
          WRITE(*,*)"HERE TESTSSS"
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
          if (isendid(irepeat) == 1) then
            NP=IPSEND(IS)
            DO WHILE (NP.NE.0)
              nsendactualp=nsendactualp+1

!             Uniform mesh - Same as in version 5.0
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
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0
 
              ypart=y(np)
              zpart=z(np)
              if (ypart <= ye.and.ypart >= yb) then
                iye_cc=jb 
              else
                if (ypart > ye) then
                  iye_cc=je+1 
                else
                  iye_cc=jb-1 
                endif
              endif
              if (zpart <= ze.and.zpart >= zb) then
                ize_cc=kb 
              else
                if (zpart > ze) then
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
              i_source = idmap_yz(iye_cc,ize_cc)
              i_tag    = it
              call MPI_SEND(pdata,7,MPI_DOUBLE_PRECISION,&
                            i_source,i_tag,MPI_COMM_WORLD,IERR)
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

              i_tag=it
              call MPI_RECV(pdata,7,MPI_DOUBLE_PRECISION,&
                            MPI_ANY_SOURCE,I_TAG,        &
                            MPI_COMM_WORLD,STATUS2,IERR)


              if (ipstore == 0) then
                Storage_Error_p = 1
              else
              
              x (nprecv)=pdata(1)
              y (nprecv)=pdata(2)
              z (nprecv)=pdata(3)
              vx(nprecv)=pdata(4)
              vy(nprecv)=pdata(5)
              vz(nprecv)=pdata(6)
              qp(nprecv)=pdata(7)

!             Uniform mesh - Same as in version 5.0
!              ixe=hxi*x(nprecv)+1.5000000000000001d+00
!              iye=hyi*y(nprecv)+0.5000000000000001d+00
!              ize=hzi*z(nprecv)+0.5000000000000001d+00

!             Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000d+00
              rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000d+00
              rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000d+00
              ixe=rxe
              iye=rye
              ize=rze
              iye=iye-1             ! integer index in y direction starts at 0
              ize=ize-1             ! integer index in z direction starts at 0
 
              ipstore=link(nprecv)

              if ((ixe > nx+1  .or. ixe < 1 ) .or. (iye > je+1    .or. iye < jb-1) .or. (ize > ke+1    .or. ize < kb-1)) then
                Field_Diverge_p = 1
                ixe = min(max(iye,1   ),nx+1)
                iye = min(max(iye,jb-1),je+1)
                ize = min(max(ize,kb-1),ke+1)
              endif

              link(nprecv)=iphead(ixe,iye,ize,is)
              iphead(ixe,iye,ize,is)=nprecv

              endif

            enddo
          endif
        enddo

        dpacket(1) = Storage_Error_p
        dpacket(2) = Field_Diverge_p
        dpacket(3) = nsendactualp
        dpacket(4) = nrecvactualp
        call MPI_ALLREDUCE(dpacket,dpacketg,4,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,IERR)
        Storage_Error = dpacketg(1)
        Field_Diverge = dpacketg(2)
        nsendactual   = dpacketg(3)
        nrecvactual   = dpacketg(4)


        if (Storage_Error /= 0) then
           if (myid == 0) then
             write(6,*)" "
             write(6,*)" "
             write(6,*) "Particle storage allocation is exceeded."
             write(6,*) "3DHybrid is stopped"
             write(6,*)" "
             write(6,*)" "
           endif
           call MPI_FINALIZE(IERR)
           STOP
        endif

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
        if (myid == 0) then
          write(6,*) " FINISHED EXCHANGING PARTICLES "
          write(6,*) " # OF PARTICLES       SENT     = ",NSENDACTUAL
          write(6,*) " # OF PARTICLES       RECEIVED = ",NRECVACTUAL
        endif
!
!=======================================================================
!
        NPTOTP=0
        indx=1
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
          WRITE(6,*) " TOTAL # OF PARTICLES AFTER  PARMOV = ",NPTOT
        ENDIF
 999    CONTINUE
 
        call date_and_time(values=time_end_array(:,14))
        call accumulate_time_difference(time_begin_array(1,14),time_end_array(1,14),time_elapsed(14))
!
!*****************************************************************************
!
        call date_and_time(values=time_begin_array(:,15))
        if (testorbt) goto 10
        NPTOTP = 0
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

                q_p = qp(l)


!               Uniform mesh - Same as in version 5.0
!                rx=hxi*x(l)+1.5000000000000001d+00
!                ry=hyi*y(l)+0.5000000000000001d+00
!                rz=hzi*z(l)+0.5000000000000001d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                ix=max(1,min(ix,nx1))
!                iy=max(jb-1,min(iy,je))
!                iz=max(kb-1,min(iz,ke))
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

!               Nonuniform mesh - without using MESH_UNMAP
!                rx=hxi*x(l)+1.500000000000000d+00
!                ry=hyi*y(l)+1.500000000000000d+00
!                rz=hzi*z(l)+1.500000000000000d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                ix=ixc_2_c_map(ix)
!                iy=iyc_2_c_map(iy)
!                iz=izc_2_c_map(iz)
!                fx=(x(l)-meshX%xc(ix))/meshX%dxn(ix+1)
!                fy=(y(l)-meshY%xc(iy))/meshY%dxn(iy+1)
!                fz=(z(l)-meshZ%xc(iz))/meshZ%dxn(iz+1)

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


!                 Diagnosic test on particle cell index algorithm
!                  if (     fx < 0. .or. fx > 1.                &
!                      .or. fy < 0. .or. fy > 1.                &
!                      .or. fz < 0. .or. fz > 1.) then 
!                      write(6,*) " GATHER LOOP "
!                      write(6,*) fx,fy,fz
!                      write(6,*) " x;",x(l),meshX%xn(ix),meshX%xc(ix)
!                      write(6,*) " y;",y(l),meshY%xn(iy),meshY%xc(iy)
!                      write(6,*) " z;",z(l),meshZ%xn(iz),meshZ%xc(iz)
!                      write(6,*) " r; ",rx,ry,rz
!                      write(6,*) " i; ",ixe,iye,ize
!                      write(6,*) " ixc_map; ",ix,iy,iz
!                  endif
 
                w1=q_p*(1.-fx)*(1.-fy)*(1.-fz)
                w2=q_p*fx     *(1.-fy)*(1.-fz)
                w3=q_p*(1.-fx)*fy     *(1.-fz)
                w4=q_p*fx     *fy     *(1.-fz)
                w5=q_p*(1.-fx)*(1.-fy)*fz
                w6=q_p*fx     *(1.-fy)*fz
                w7=q_p*(1.-fx)*fy     *fz
                w8=q_p*fx     *fy     *fz
 
                dns(ix  ,iy  ,iz  ,is)=dns(ix  ,iy  ,iz  ,is)+w1
                dns(ixp1,iy  ,iz  ,is)=dns(ixp1,iy  ,iz  ,is)+w2
                dns(ix  ,iyp1,iz  ,is)=dns(ix  ,iyp1,iz  ,is)+w3
                dns(ixp1,iyp1,iz  ,is)=dns(ixp1,iyp1,iz  ,is)+w4
                dns(ix  ,iy  ,izp1,is)=dns(ix  ,iy  ,izp1,is)+w5
                dns(ixp1,iy  ,izp1,is)=dns(ixp1,iy  ,izp1,is)+w6
                dns(ix  ,iyp1,izp1,is)=dns(ix  ,iyp1,izp1,is)+w7
                dns(ixp1,iyp1,izp1,is)=dns(ixp1,iyp1,izp1,is)+w8
 
                vxs(ix  ,iy  ,iz  ,is)=vxs(ix  ,iy  ,iz  ,is)+w1*vx(l)
                vxs(ixp1,iy  ,iz  ,is)=vxs(ixp1,iy  ,iz  ,is)+w2*vx(l)
                vxs(ix  ,iyp1,iz  ,is)=vxs(ix  ,iyp1,iz  ,is)+w3*vx(l)
                vxs(ixp1,iyp1,iz  ,is)=vxs(ixp1,iyp1,iz  ,is)+w4*vx(l)
                vxs(ix  ,iy  ,izp1,is)=vxs(ix  ,iy  ,izp1,is)+w5*vx(l)
                vxs(ixp1,iy  ,izp1,is)=vxs(ixp1,iy  ,izp1,is)+w6*vx(l)
                vxs(ix  ,iyp1,izp1,is)=vxs(ix  ,iyp1,izp1,is)+w7*vx(l)
                vxs(ixp1,iyp1,izp1,is)=vxs(ixp1,iyp1,izp1,is)+w8*vx(l)

                vys(ix  ,iy  ,iz  ,is)=vys(ix  ,iy  ,iz  ,is)+w1*vy(l)
                vys(ixp1,iy  ,iz  ,is)=vys(ixp1,iy  ,iz  ,is)+w2*vy(l)
                vys(ix  ,iyp1,iz  ,is)=vys(ix  ,iyp1,iz  ,is)+w3*vy(l)
                vys(ixp1,iyp1,iz  ,is)=vys(ixp1,iyp1,iz  ,is)+w4*vy(l)
                vys(ix  ,iy  ,izp1,is)=vys(ix  ,iy  ,izp1,is)+w5*vy(l)
                vys(ixp1,iy  ,izp1,is)=vys(ixp1,iy  ,izp1,is)+w6*vy(l)
                vys(ix  ,iyp1,izp1,is)=vys(ix  ,iyp1,izp1,is)+w7*vy(l)
                vys(ixp1,iyp1,izp1,is)=vys(ixp1,iyp1,izp1,is)+w8*vy(l)

                vzs(ix  ,iy  ,iz  ,is)=vzs(ix  ,iy  ,iz  ,is)+w1*vz(l)
                vzs(ixp1,iy  ,iz  ,is)=vzs(ixp1,iy  ,iz  ,is)+w2*vz(l)
                vzs(ix  ,iyp1,iz  ,is)=vzs(ix  ,iyp1,iz  ,is)+w3*vz(l)
                vzs(ixp1,iyp1,iz  ,is)=vzs(ixp1,iyp1,iz  ,is)+w4*vz(l)
                vzs(ix  ,iy  ,izp1,is)=vzs(ix  ,iy  ,izp1,is)+w5*vz(l)
                vzs(ixp1,iy  ,izp1,is)=vzs(ixp1,iy  ,izp1,is)+w6*vz(l)
                vzs(ix  ,iyp1,izp1,is)=vzs(ix  ,iyp1,izp1,is)+w7*vz(l)
                vzs(ixp1,iyp1,izp1,is)=vzs(ixp1,iyp1,izp1,is)+w8*vz(l)
 
                np=link(np)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

 10     kspc=is

        nescapearr(1) = nescape(is)
        nescapearr(2) = nescape_yz(is)
        nescapearr(3) = nescape_zy(is)
        nescapearr(4) = nescape_xy(is)
        nescapearr(5) = nescape_yx(is)
        nescapearr(6) = nescape_zx(is)
        nescapearr(7) = nescape_xz(is)
        nescapearr(8) = npart(is)
        call MPI_ALLREDUCE(nescapearr,nescapearr_global,8,MPI_INTEGER4,&
                         MPI_SUM,COMM2D,IERR)
        nescape_global(is)    = nescapearr_global(1)
        nescape_yz_global(is) = nescapearr_global(2)
        nescape_zy_global(is) = nescapearr_global(3)
        nescape_xy_global(is) = nescapearr_global(4)
        nescape_yx_global(is) = nescapearr_global(5)
        nescape_zx_global(is) = nescapearr_global(6)
        nescape_xz_global(is) = nescapearr_global(7)
        npart_global(is)      = nescapearr_global(8)

        deltime2 = deltime2 + real(clock_time1-clock_time)
 
        call XREAL(DNS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(VXS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(VYS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(VZS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREALBCC(DNS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC(VXS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC(VYS(1,jb-1,kb-1,is),1,NX,NY,NZ)
        call XREALBCC(VZS(1,jb-1,kb-1,is),1,NX,NY,NZ)


        DO IIZ=KB-1,KE+1
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
 
 
      ENDDO                ! IS DO LOOP
!
!     end of main particle loop
!

      epacket(1) = nptotp
      epacket(2) = npleavingp
      call MPI_ALLREDUCE(epacket,epacketg,2,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
      nptot      = epacketg(1)
      npleaving  = epacketg(2)

      if (prntinfo) then
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
      endif

      if (myid == 0.and.prntinfo) then
        do is=1,nspec

          if (is == 1) then
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
          write(6,1005) ninj_global(1)+ninj_global(2),nescape_global(1)+nescape_global(2),npart_global(1)+npart_global(2)
        else
          write(6,1005) ninj_global(1),nescape_global(1),npart_global(1)
        endif

 1005   format(5x,'sum',4x,i8,2x,i8,2x,i10)
 1006   format(2i6,4x,'sum',4x,i8,2x,i8,2x,i10)
      endif


      ninj        = 0
      ninj_global = 0
 
      call date_and_time(values=time_end_array(:,19))
      call accumulate_time_difference(time_begin_array(1,19),time_end_array(1,19),time_elapsed(19))
 

      return
      end
!
!######################################################################
!
      subroutine sortit
!
!=======================================================================
!
      use parameter_mod
      use MESH2D
      double precision pstore(nplmax)
      integer*8 id, kb1, is, ix, iy, iz, ixe ,iye, ize, l, nttot, nplist
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
 
      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt
 
      id = 0
      kb1 = kb-1
      iptemp = 0
      porder = 0
      do is = 1,nspec
        do iz = kb-1,ke
          do iy = jb-1,je
            do ix = 1, nx1
              np = iphead(ix,iy,iz,is)
              DO WHILE (NP.NE.0)

!               Uniform mesh - Same as in version 5.0
!                ixe = hxi*x(np)+1.5000000000000001d+00
!                iye = hyi*y(np)+0.5000000000000001d+00
!                ize = hzi*z(np)+0.5000000000000001d+00

!                 Nonuniform mesh - using MESH_UNMAP
                  rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
                  rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
                  rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
                  ixe=rxe
                  iye=rye
                  ize=rze
                  iye=iye-1             ! integer index in y direction starts at 0
                  ize=ize-1             ! integer index in z direction starts at 0
 
                porder(np)=iptemp(ixe,iye,ize,is)
                iptemp(ixe,iye,ize,is)=np
                np = link(np)
              ENDDO
            enddo
          enddo
        enddo
      enddo

      l = 0
      nttot = 0
      do is = 1,nspec
        do iz = kb-1,ke
          do iy = jb-1,je
            do ix = 1, nx1
              np = iptemp(ix,iy,iz,is)
              nplist = 0
              DO WHILE (NP.NE.0)
                nplist = nplist+1
                l = l+1
                link(l) = np
                np = porder(np)
              ENDDO
              nttot = nttot + nplist
              iphead(ix,iy,iz,is) = nplist
            enddo
          enddo
        enddo
      enddo
 
 
      id = 0
      kb1 = kb-1
      do l = 1,nttot
        pstore(l) = vx(link(l))
      enddo
      do l = 1,nttot
        vx(l) = pstore(l)
      enddo


      do l = 1,nttot
        pstore(l) = vy(link(l))
      enddo
      do l = 1,nttot
        vy(l) = pstore(l)
      enddo      

	  
      do l = 1,nttot
        pstore(l) = vz(link(l))
      enddo
      do l = 1,nttot
        vz(l) = pstore(l)
      enddo      

	  
      do l = 1,nttot
        pstore(l) = x(link(l))
      enddo
      do l = 1,nttot
        x(l) = pstore(l)
      enddo      

	  
      do l = 1,nttot
        pstore(l) = y(link(l))
      enddo
      do l = 1,nttot
        y(l) = pstore(l)
      enddo      

	  
      do l = 1,nttot
        pstore(l) = z(link(l))
      enddo
      do l = 1,nttot
        z(l) = pstore(l)
      enddo

	  
      do l = 1,nttot
        pstore(l) = qp(link(l))
      enddo
      do l = 1,nttot
        qp(l) = pstore(l)
      enddo


      l=1
      do is = 1,nspec
        do iz = kb-1,ke
          do iy = jb-1,je
            do ix = 1, nx1
	      nplist = iphead(ix,iy,iz,is)
	      if (nplist.ne.0) then
                iphead(ix,iy,iz,is) = l
                do np = l, l+nplist-1
                  link(np) = np+1
                enddo
                link(l+nplist-1) = 0
                l = l + nplist
              else
                iphead(ix,iy,iz,is) = 0
              endif
            enddo
          enddo
        enddo
      enddo

      if (l-1.ne.nttot) then
        print *,'Problem in SORT: l-1 NE NTTOT'
        stop
      endif
      ipstore = nttot + 1
      do l = nttot+1, nplmax-1
        link(l) = l+1
      enddo
      link(nplmax) = 0

      return
      end

      subroutine clock_write(iunit,message,i2,i1,is,it)
      integer*8 iunit,i2,i1,is
      character*10 message
      write(iunit,"(i4,a10,e20.8)") it, real(i2-i1)/real(is)
      return
      end

!
!#######################################################################
!
      subroutine field
 
      use parameter_mod
  
!     Bypass field solver for timing on LANL machine
!      goto 999 


      call date_and_time(values=time_begin_array(:,21))
 
 
      call date_and_time(values=time_begin_array(:,9))
      call pressgrad(1)
      call date_and_time(values=time_end_array(:,9))
      call accumulate_time_difference(time_begin_array(1,9) &
     &                               ,time_end_array(1,9) &
     &                               ,time_elapsed(9))
 
 
      call date_and_time(values=time_begin_array(:,10))
      call bcalc
      call date_and_time(values=time_end_array(:,10))
      call accumulate_time_difference(time_begin_array(1,10) &
     &                               ,time_end_array(1,10) &
     &                               ,time_elapsed(10))
 
 
      call date_and_time(values=time_begin_array(:,9))
      call pressgrad(0)
      call date_and_time(values=time_end_array(:,9))
      call accumulate_time_difference(time_begin_array(1,9) &
     &                               ,time_end_array(1,9) &
     &                               ,time_elapsed(9))
 
 
      call date_and_time(values=time_begin_array(:,11))
      call ecalc( 0 )
      call date_and_time(values=time_end_array(:,11))
      call accumulate_time_difference(time_begin_array(1,11) &
     &                               ,time_end_array(1,11) &
     &                               ,time_elapsed(11))
 
 
      call date_and_time(values=time_begin_array(:,12))
      call focalc
      call date_and_time(values=time_end_array(:,12))
      call accumulate_time_difference(time_begin_array(1,12) &
     &                               ,time_end_array(1,12) &
     &                               ,time_elapsed(12))
 
 
      call date_and_time(values=time_end_array(:,21))
      call accumulate_time_difference(time_begin_array(1,21) &
     &                               ,time_end_array(1,21) &
     &                               ,time_elapsed(21))
 

 999  continue

      return
      end
!
!################################################################################
!
      subroutine pressgrad(iflag)

      use parameter_mod
      use MESH2D

      do k=kb,ke
        do j = jb,je
          do i=2,nx1
            dena=iflag*0.5*(den(i,j,k)+deno(i,j,k))&
                +(1.-iflag)*den(i,j,k)
            a=1/dena

!           Uniform mesh - Same as in version 5.0
!            dxa=a/(4.*hx)
!            dya=a/(4.*hy)
!            dza=a/(4.*hz)

!           Nonuniform mesh
!            dxa=a/(2.*(meshX%dxn(i  )+meshX%dxn(i+1)))
!            dya=a/(2.*(meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
!            dza=a/(2.*(meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0
            dxa=a/(4.*(meshX%dxn(i  )+meshX%dxn(i+1)))
            dya=a/(4.*(meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
            dza=a/(4.*(meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0


           dpedx(i,j,k)=((pe(i+1,j-1,k+1)+2.*pe(i+1,j,k+1)&
               +pe(i+1,j+1,k+1))/4.&
               +2.*(pe(i+1,j-1,k  )+2.*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4.&
               +   (pe(i+1,j-1,k-1)+2.*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4.&
               -   (pe(i-1,j-1,k+1)+2.*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4.&
               -2.*(pe(i-1,j-1,k  )+2.*pe(i-1,j,k  )+pe(i-1,j+1,k  ))/4.&
               -   (pe(i-1,j-1,k-1)+2.*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4.&
                     )*dxa
             dpedy(i,j,k)=((pe(i-1,j+1,k+1)+2.*pe(i,j+1,k+1)&
                +pe(i+1,j+1,k+1))/4.&
                +2.*(pe(i-1,j+1,k  )+2.*pe(i,j+1,k  )+pe(i+1,j+1,k  ))/4.&
                +   (pe(i-1,j+1,k-1)+2.*pe(i,j+1,k-1)+pe(i+1,j+1,k-1))/4.&
                -   (pe(i-1,j-1,k+1)+2.*pe(i,j-1,k+1)+pe(i+1,j-1,k+1))/4.&
                -2.*(pe(i-1,j-1,k  )+2.*pe(i,j-1,k  )+pe(i+1,j-1,k  ))/4.&
                -   (pe(i-1,j-1,k-1)+2.*pe(i,j-1,k-1)+pe(i+1,j-1,k-1))/4.&
                     )*dya
             dpedz(i,j,k)=((pe(i+1,j-1,k+1)+2.*pe(i+1,j,k+1)&
                +pe(i+1,j+1,k+1))/4.&
                +2.*(pe(i  ,j-1,k+1)+2.*pe(i  ,j,k+1)+pe(i  ,j+1,k+1))/4.&
                +   (pe(i-1,j-1,k+1)+2.*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4.&
                -   (pe(i+1,j-1,k-1)+2.*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4.&
                -2.*(pe(i  ,j-1,k-1)+2.*pe(i  ,j,k-1)+pe(i  ,j+1,k-1))/4.&
                -   (pe(i-1,j-1,k-1)+2.*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4.&
                     )*dza                    
          enddo
        enddo
      enddo
 

      return
      end
!
!################################################################################
!
      subroutine ecalc( iflag )

      use parameter_mod
      use MESH2D
      double precision:: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &
                        ,curr_tot
      double precision:: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8 
      double precision:: by1,by2,by3,by4,by5,by6,by7,by8 
      double precision:: bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8  
      double precision:: vixa, viya, viza, dena, a, dxa, dya, dza 
      double precision:: dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy 
      double precision:: curlbx_scalar,curlby_scalar,curlbz_scalar
      double precision:: bxav, byav, bzav  
      double precision:: dexdy, dexdz, deydx, deydz, dezdx,dezdy  
 
      if (eta_par == 0) then
      do k = kb,ke
        do j = jb,je
          do i=2,nx1
            bx1=bx(i+1,j+1,k)  +bdipole_x(i+1,j+1,k)
            bx2=bx(i  ,j+1,k)  +bdipole_x(i  ,j+1,k)
            bx3=bx(i  ,j  ,k)  +bdipole_x(i  ,j  ,k)
            bx4=bx(i+1,j  ,k)  +bdipole_x(i+1,j  ,k)
            bx5=bx(i+1,j+1,k+1)+bdipole_x(i+1,j+1,k+1)
            bx6=bx(i  ,j+1,k+1)+bdipole_x(i  ,j+1,k+1)
            bx7=bx(i  ,j  ,k+1)+bdipole_x(i  ,j  ,k+1)
            bx8=bx(i+1,j  ,k+1)+bdipole_x(i+1,j  ,k+1)
            by1=by(i+1,j+1,k)  +bdipole_y(i+1,j+1,k)
            by2=by(i  ,j+1,k)  +bdipole_y(i  ,j+1,k)
            by3=by(i  ,j  ,k)  +bdipole_y(i  ,j  ,k)
            by4=by(i+1,j  ,k)  +bdipole_y(i+1,j  ,k)
            by5=by(i+1,j+1,k+1)+bdipole_y(i+1,j+1,k+1)
            by6=by(i  ,j+1,k+1)+bdipole_y(i  ,j+1,k+1)
            by7=by(i  ,j  ,k+1)+bdipole_y(i  ,j  ,k+1)
            by8=by(i+1,j  ,k+1)+bdipole_y(i+1,j  ,k+1)
            bz1=bz(i+1,j+1,k)  +bdipole_z(i+1,j+1,k)
            bz2=bz(i  ,j+1,k)  +bdipole_z(i  ,j+1,k)
            bz3=bz(i  ,j  ,k)  +bdipole_z(i  ,j  ,k)
            bz4=bz(i+1,j  ,k)  +bdipole_z(i+1,j  ,k)
            bz5=bz(i+1,j+1,k+1)+bdipole_z(i+1,j+1,k+1)
            bz6=bz(i  ,j+1,k+1)+bdipole_z(i  ,j+1,k+1)
            bz7=bz(i  ,j  ,k+1)+bdipole_z(i  ,j  ,k+1)
            bz8=bz(i+1,j  ,k+1)+bdipole_z(i+1,j  ,k+1)

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
!            dxa=a/(4.*hx)
!            dya=a/(4.*hy)
!            dza=a/(4.*hz)

!          Nonuniform mesh
            dxa=a/(4.*meshX%dxc(i))
            dya=a/(4.*meshY%dxc(j+1))                  ! integer index in y direction starts at 0
            dza=a/(4.*meshZ%dxc(k+1))                  ! integer index in z direction starts at 0

            dbxdy= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                 +bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
                 -bx(i+1,j  ,k+1)-bx(i  ,j  ,k+1)&
                 -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
            dbxdz= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                 +bx(i  ,j  ,k+1)+bx(i+1,j  ,k+1)&
                 -bx(i+1,j+1,k  )-bx(i  ,j+1,k  )&
                 -bx(i  ,j  ,k  )-bx(i+1,j ,k  )
            dbydx= by(i+1,j+1,k+1)+by(i+1,j  ,k+1)&
                 +by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                 -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)&
                 -by(i  ,j  ,k  )-by(i  ,j+1,k  )
            dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)&
                 +by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
                 -by(i+1,j+1,k  )-by(i  ,j+1,k  )&
                 -by(i  ,j  ,k  )-by(i+1,j  ,k  )
            dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)&
                 +bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                 -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)&
                 -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
            dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)&
                 +bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                 -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)&
                 -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
            curlbx_scalar=dya*dbzdy-dza*dbydz
            curlby_scalar=dza*dbxdz-dxa*dbzdx
            curlbz_scalar=dxa*dbydx-dya*dbxdy
            bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar

! From the eta_par conditional
            tenx=eta(i,j,k)*xj
            teny=eta(i,j,k)*yj
            tenz=eta(i,j,k)*zj
! End content from the eta_par conditional

            ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                     -dpedx(i,j,k)+tenx/a
            ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                     -dpedy(i,j,k)+teny/a
            ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                     -dpedz(i,j,k)+tenz/a
          enddo
        enddo
      enddo
      else
       if (eta_par == 1) then
       do k = kb,ke
         do j = jb,je
           do i=2,nx1
            bx1=bx(i+1,j+1,k)  +bdipole_x(i+1,j+1,k)
            bx2=bx(i  ,j+1,k)  +bdipole_x(i  ,j+1,k)
            bx3=bx(i  ,j  ,k)  +bdipole_x(i  ,j  ,k)
            bx4=bx(i+1,j  ,k)  +bdipole_x(i+1,j  ,k)
            bx5=bx(i+1,j+1,k+1)+bdipole_x(i+1,j+1,k+1)
            bx6=bx(i  ,j+1,k+1)+bdipole_x(i  ,j+1,k+1)
            bx7=bx(i  ,j  ,k+1)+bdipole_x(i  ,j  ,k+1)
            bx8=bx(i+1,j  ,k+1)+bdipole_x(i+1,j  ,k+1)
            by1=by(i+1,j+1,k)  +bdipole_y(i+1,j+1,k)
            by2=by(i  ,j+1,k)  +bdipole_y(i  ,j+1,k)
            by3=by(i  ,j  ,k)  +bdipole_y(i  ,j  ,k)
            by4=by(i+1,j  ,k)  +bdipole_y(i+1,j  ,k)
            by5=by(i+1,j+1,k+1)+bdipole_y(i+1,j+1,k+1)
            by6=by(i  ,j+1,k+1)+bdipole_y(i  ,j+1,k+1)
            by7=by(i  ,j  ,k+1)+bdipole_y(i  ,j  ,k+1)
            by8=by(i+1,j  ,k+1)+bdipole_y(i+1,j  ,k+1)
            bz1=bz(i+1,j+1,k)  +bdipole_z(i+1,j+1,k)
            bz2=bz(i  ,j+1,k)  +bdipole_z(i  ,j+1,k)
            bz3=bz(i  ,j  ,k)  +bdipole_z(i  ,j  ,k)
            bz4=bz(i+1,j  ,k)  +bdipole_z(i+1,j  ,k)
            bz5=bz(i+1,j+1,k+1)+bdipole_z(i+1,j+1,k+1)
            bz6=bz(i  ,j+1,k+1)+bdipole_z(i  ,j+1,k+1)
            bz7=bz(i  ,j  ,k+1)+bdipole_z(i  ,j  ,k+1)
            bz8=bz(i+1,j  ,k+1)+bdipole_z(i+1,j  ,k+1)

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
!            dxa=a/(4.*hx)
!            dya=a/(4.*hy)
!            dza=a/(4.*hz)

!          Nonuniform mesh
            dxa=a/(4.*meshX%dxc(i))
            dya=a/(4.*meshY%dxc(j+1))                  ! integer index in y direction starts at 0
            dza=a/(4.*meshZ%dxc(k+1))                  ! integer index in z direction starts at 0

            dbxdy= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                 +bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
                 -bx(i+1,j  ,k+1)-bx(i  ,j  ,k+1)&
                 -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
            dbxdz= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                 +bx(i  ,j  ,k+1)+bx(i+1,j  ,k+1)&
                 -bx(i+1,j+1,k  )-bx(i  ,j+1,k  )&
                 -bx(i  ,j  ,k  )-bx(i+1,j ,k  )
            dbydx= by(i+1,j+1,k+1)+by(i+1,j  ,k+1)&
                 +by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                 -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)&
                 -by(i  ,j  ,k  )-by(i  ,j+1,k  )
            dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)&
                 +by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
                 -by(i+1,j+1,k  )-by(i  ,j+1,k  )&
                 -by(i  ,j  ,k  )-by(i+1,j  ,k  )
            dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)&
                 +bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                 -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)&
                 -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
            dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)&
                 +bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                 -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)&
                 -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
            curlbx_scalar=dya*dbzdy-dza*dbydz
            curlby_scalar=dza*dbxdz-dxa*dbzdx
            curlbz_scalar=dxa*dbydx-dya*dbxdy
            bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar

! From eta_par conditional
            bxx = bxav
            byy = byav
            bzz = bzav
            btot = sqrt(bxx**2 + byy**2 + bzz**2)
            tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
            tenx = tjdotb*bxx/btot
            teny = tjdotb*byy/btot
            tenz = tjdotb*bzz/btot
! End content from eta_par conditional

            ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                     -dpedx(i,j,k)+tenx/a
            ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                     -dpedy(i,j,k)+teny/a
            ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                     -dpedz(i,j,k)+tenz/a
           enddo
         enddo
       enddo
       else if (eta_par == 2) then
       do k = kb,ke
         do j = jb,je
           do i=2,nx1
            bx1=bx(i+1,j+1,k)  +bdipole_x(i+1,j+1,k)
            bx2=bx(i  ,j+1,k)  +bdipole_x(i  ,j+1,k)
            bx3=bx(i  ,j  ,k)  +bdipole_x(i  ,j  ,k)
            bx4=bx(i+1,j  ,k)  +bdipole_x(i+1,j  ,k)
            bx5=bx(i+1,j+1,k+1)+bdipole_x(i+1,j+1,k+1)
            bx6=bx(i  ,j+1,k+1)+bdipole_x(i  ,j+1,k+1)
            bx7=bx(i  ,j  ,k+1)+bdipole_x(i  ,j  ,k+1)
            bx8=bx(i+1,j  ,k+1)+bdipole_x(i+1,j  ,k+1)
            by1=by(i+1,j+1,k)  +bdipole_y(i+1,j+1,k)
            by2=by(i  ,j+1,k)  +bdipole_y(i  ,j+1,k)
            by3=by(i  ,j  ,k)  +bdipole_y(i  ,j  ,k)
            by4=by(i+1,j  ,k)  +bdipole_y(i+1,j  ,k)
            by5=by(i+1,j+1,k+1)+bdipole_y(i+1,j+1,k+1)
            by6=by(i  ,j+1,k+1)+bdipole_y(i  ,j+1,k+1)
            by7=by(i  ,j  ,k+1)+bdipole_y(i  ,j  ,k+1)
            by8=by(i+1,j  ,k+1)+bdipole_y(i+1,j  ,k+1)
            bz1=bz(i+1,j+1,k)  +bdipole_z(i+1,j+1,k)
            bz2=bz(i  ,j+1,k)  +bdipole_z(i  ,j+1,k)
            bz3=bz(i  ,j  ,k)  +bdipole_z(i  ,j  ,k)
            bz4=bz(i+1,j  ,k)  +bdipole_z(i+1,j  ,k)
            bz5=bz(i+1,j+1,k+1)+bdipole_z(i+1,j+1,k+1)
            bz6=bz(i  ,j+1,k+1)+bdipole_z(i  ,j+1,k+1)
            bz7=bz(i  ,j  ,k+1)+bdipole_z(i  ,j  ,k+1)
            bz8=bz(i+1,j  ,k+1)+bdipole_z(i+1,j  ,k+1)

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
!            dxa=a/(4.*hx)
!            dya=a/(4.*hy)
!            dza=a/(4.*hz)

!          Nonuniform mesh
            dxa=a/(4.*meshX%dxc(i))
            dya=a/(4.*meshY%dxc(j+1))                  ! integer index in y direction starts at 0
            dza=a/(4.*meshZ%dxc(k+1))                  ! integer index in z direction starts at 0

            dbxdy= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                 +bx(i  ,j+1,k  )+bx(i+1,j+1,k  )&
                 -bx(i+1,j  ,k+1)-bx(i  ,j  ,k+1)&
                 -bx(i  ,j  ,k  )-bx(i+1,j  ,k  )
            dbxdz= bx(i+1,j+1,k+1)+bx(i  ,j+1,k+1)&
                 +bx(i  ,j  ,k+1)+bx(i+1,j  ,k+1)&
                 -bx(i+1,j+1,k  )-bx(i  ,j+1,k  )&
                 -bx(i  ,j  ,k  )-bx(i+1,j ,k  )
            dbydx= by(i+1,j+1,k+1)+by(i+1,j  ,k+1)&
                 +by(i+1,j  ,k  )+by(i+1,j+1,k  )&
                 -by(i  ,j+1,k+1)-by(i  ,j  ,k+1)&
                 -by(i  ,j  ,k  )-by(i  ,j+1,k  )
            dbydz= by(i+1,j+1,k+1)+by(i  ,j+1,k+1)&
                 +by(i  ,j  ,k+1)+by(i+1,j  ,k+1)&
                 -by(i+1,j+1,k  )-by(i  ,j+1,k  )&
                 -by(i  ,j  ,k  )-by(i+1,j  ,k  )
            dbzdx= bz(i+1,j+1,k+1)+bz(i+1,j  ,k+1)&
                 +bz(i+1,j  ,k  )+bz(i+1,j+1,k  )&
                 -bz(i  ,j+1,k+1)-bz(i  ,j  ,k+1)&
                 -bz(i  ,j  ,k  )-bz(i  ,j+1,k  )
            dbzdy= bz(i+1,j+1,k+1)+bz(i  ,j+1,k+1)&
                 +bz(i  ,j+1,k  )+bz(i+1,j+1,k  )&
                 -bz(i+1,j  ,k+1)-bz(i  ,j  ,k+1)&
                 -bz(i  ,j  ,k  )-bz(i+1,j  ,k  )
            curlbx_scalar=dya*dbzdy-dza*dbydz
            curlby_scalar=dza*dbxdz-dxa*dbzdx
            curlbz_scalar=dxa*dbydx-dya*dbxdy
            bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar

! From eta_par conditional
            bxx = bxav
            byy = byav
            bzz = bzav
            btot = sqrt(bxx**2 + byy**2 + bzz**2)
            curr_tot = max(1.d-12,sqrt(xj**2 + yj**2 + zj**2))
            tenx = abs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
            tenx = min(resis,tenx)
            tenx = tenx*xj
            teny = abs(eta(i,j,k)*byy*yj/(btot*curr_tot))
            teny = min(resis,teny)
            teny = teny*yj
            tenz = abs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
            tenz = min(resis,tenz)
            tenz = tenz*zj
! End content from eta_par conditional

            ex(i,j,k)=(viza*byav-viya*bzav)+(curlby_scalar*bzav-curlbz_scalar*byav)&
                     -dpedx(i,j,k)+tenx/a
            ey(i,j,k)=(vixa*bzav-viza*bxav)+(curlbz_scalar*bxav-curlbx_scalar*bzav)&
                     -dpedy(i,j,k)+teny/a
            ez(i,j,k)=(viya*bxav-vixa*byav)+(curlbx_scalar*byav-curlby_scalar*bxav)&
                     -dpedz(i,j,k)+tenz/a
           enddo
         enddo
        enddo
       endif
      endif
      ex=ex*dipole_sphere_ex
      ey=ey*dipole_sphere_ey
      ez=ez*dipole_sphere_ez
!
!******************************************
!
!     boundary conditions
!
      call date_and_time(values=time_begin_array(:,18))
      call XREALBCC_PACK_E(EX,EY,EZ,1,NX,NY,NZ)
      call date_and_time(values=time_end_array(:,18))
      call accumulate_time_difference(time_begin_array(1,18) &
     &                               ,time_end_array(1,18) &
     &                               ,time_elapsed(18))

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
!***********************************************************
!
      do k=kb,ke+1
        do j = jb,je+1
          do i=2,nx2
            dexdy=  ex(i  ,j  ,k  ) + ex(i-1,j  ,k  )   &
                  + ex(i-1,j  ,k-1) + ex(i  ,j  ,k-1)   &
                  - ex(i  ,j-1,k  ) - ex(i-1,j-1,k  )   &
                  - ex(i-1,j-1,k-1) - ex(i  ,j-1,k-1)
            dexdz=  ex(i  ,j  ,k  ) + ex(i-1,j  ,k  )   &
                  + ex(i-1,j-1,k  ) + ex(i  ,j-1,k  )   &
                  - ex(i  ,j  ,k-1) - ex(i-1,j  ,k-1)   &
                  - ex(i-1,j-1,k-1) - ex(i  ,j-1,k-1)
            deydx=  ey(i  ,j  ,k  ) + ey(i  ,j-1,k  )   &
                  + ey(i  ,j-1,k-1) + ey(i  ,j  ,k-1)   &
                  - ey(i-1,j  ,k  ) - ey(i-1,j-1,k  )   &
                  - ey(i-1,j-1,k-1) - ey(i-1,j  ,k-1)
            deydz=  ey(i  ,j  ,k  ) + ey(i-1,j  ,k  )   &
                  + ey(i-1,j-1,k  ) + ey(i  ,j-1,k  )   &
                  - ey(i  ,j  ,k-1) - ey(i-1,j  ,k-1)   &
                  - ey(i-1,j-1,k-1) - ey(i  ,j-1,k-1)
            dezdx=  ez(i  ,j  ,k  ) + ez(i  ,j-1,k  )   &
                  + ez(i  ,j-1,k-1) + ez(i  ,j  ,k-1)   &
                  - ez(i-1,j  ,k  ) - ez(i-1,j-1,k  )   &
                  - ez(i-1,j-1,k-1) - ez(i-1,j  ,k-1)
            dezdy=  ez(i  ,j  ,k  ) + ez(i-1,j  ,k  )   &
                  + ez(i-1,j  ,k-1) + ez(i  ,j  ,k-1)   &
                  - ez(i  ,j-1,k  ) - ez(i-1,j-1,k  )   &
                  - ez(i-1,j-1,k-1) - ez(i  ,j-1,k-1)

!           Uniform mesh - Same as is in version 5.0
!            curlex(i,j,k)=dezdy/(4.*hy)-deydz/(4.*hz)
!            curley(i,j,k)=dexdz/(4.*hz)-dezdx/(4.*hx)
!            curlez(i,j,k)=deydx/(4.*hx)-dexdy/(4.*hy)

!          Nonuniform mesh
            curlex(i,j,k)=dezdy/(4.*meshY%dxn(j+1))-deydz/(4.*meshZ%dxn(k+1))        ! integer index in y and z directions start  at 0
            curley(i,j,k)=dexdz/(4.*meshZ%dxn(k+1))-dezdx/(4.*meshX%dxn(i  ))        ! integer index in z       direction  starts at 0
            curlez(i,j,k)=deydx/(4.*meshX%dxn(i  ))-dexdy/(4.*meshY%dxn(j+1))        ! integer index in y       direction  starts at 0

          enddo
        enddo
      enddo
!
!
      return
      end
!
!################################################################################
!
      subroutine bcalc

      use parameter_mod
      use MESH2D
      double precision:: tempx1(nxmax,jb-1:je+1,kb-1:ke+1)&
                        ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1)&
                        ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)
 
 
       call date_and_time(values=time_begin_array(:,22))
 
       dts=dt/real(iterb)
       dts2=dts/2.
       dts6=dts/6.
!
!***********************************
!   subcycle into iterb interations
!***********************************

      do 10 ii=1,iterb
 

!*******************************
!
!   save B at start of subcycle
!
!   Bs = B(n)
!*******************************
!
      bxs=bx
      bys=by
      bzs=bz
!
!*******************
!   R-K first part
!*******************
!
       call date_and_time(values=time_begin_array(:,16))
       call ecalc( 1 )
       call date_and_time(values=time_end_array(:,16))
       call accumulate_time_difference(time_begin_array(1,16) &
     &                                ,time_end_array(1,16) &
     &                                ,time_elapsed(16))
 
!***********************
!   B = B(n)+dt*K1/2
!***********************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               bx(i,j,k)=bxs(i,j,k)-dts2*curlex(i,j,k)
               by(i,j,k)=bys(i,j,k)-dts2*curley(i,j,k)
               bz(i,j,k)=bzs(i,j,k)-dts2*curlez(i,j,k)
             enddo
           enddo
         enddo
 

!******************
!   temp1 = K1
!******************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               tempx1(i,j,k)=curlex(i,j,k)
               tempy1(i,j,k)=curley(i,j,k)
               tempz1(i,j,k)=curlez(i,j,k)
             enddo
           enddo
         enddo

!***************
!   R-K part 2
!***************

       call date_and_time(values=time_begin_array(:,16))
       call ecalc( 1 )
       call date_and_time(values=time_end_array(:,16))
       call accumulate_time_difference(time_begin_array(1,16) &
     &                                ,time_end_array(1,16) &
     &                                ,time_elapsed(16))
 

!*********************
!  B = B(n)+dt*K2/2
!*********************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               bx(i,j,k)=bxs(i,j,k)-dts2*curlex(i,j,k)
               by(i,j,k)=bys(i,j,k)-dts2*curley(i,j,k)
               bz(i,j,k)=bzs(i,j,k)-dts2*curlez(i,j,k)
             enddo
           enddo
         enddo
 
 
!********************
!  temp2 = K2
!********************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               tempx1(i,j,k)=tempx1(i,j,k)+2.*curlex(i,j,k)
               tempy1(i,j,k)=tempy1(i,j,k)+2.*curley(i,j,k)
               tempz1(i,j,k)=tempz1(i,j,k)+2.*curlez(i,j,k)
             enddo
           enddo
         enddo

!*****************
!  R-K  part 3
!*****************

       call date_and_time(values=time_begin_array(:,16))
       call ecalc( 1 )
       call date_and_time(values=time_end_array(:,16))
       call accumulate_time_difference(time_begin_array(1,16) &
     &                                ,time_end_array(1,16) &
     &                                ,time_elapsed(16))

!*********************
!  B = B(n)+dt*K3
!*********************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               bx(i,j,k)=bxs(i,j,k)-dts*curlex(i,j,k)
               by(i,j,k)=bys(i,j,k)-dts*curley(i,j,k)
               bz(i,j,k)=bzs(i,j,k)-dts*curlez(i,j,k)
             enddo
           enddo
         enddo
 
 
!*********************
!  temp3 = K3
!*********************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               tempx1(i,j,k)=tempx1(i,j,k)+2.*curlex(i,j,k)
               tempy1(i,j,k)=tempy1(i,j,k)+2.*curley(i,j,k)
               tempz1(i,j,k)=tempz1(i,j,k)+2.*curlez(i,j,k)
             enddo
           enddo
         enddo


!***************
!   R-K  part 4
!***************
!
       call date_and_time(values=time_begin_array(:,16))
       call ecalc( 1 )
       call date_and_time(values=time_end_array(:,16))
       call accumulate_time_difference(time_begin_array(1,16) &
     &                                ,time_end_array(1,16) &
     &                                ,time_elapsed(16))
 

!*************************************
!  B = B(n) + dt*(K1+2K2+2K3+K4)/6
!*************************************

         do k=kb,ke+1
           do j = jb,je+1
             do i=2,nx2
               bx(i,j,k)=bxs(i,j,k)-dts6*(tempx1(i,j,k)+curlex(i,j,k))
               by(i,j,k)=bys(i,j,k)-dts6*(tempy1(i,j,k)+curley(i,j,k))
               bz(i,j,k)=bzs(i,j,k)-dts6*(tempz1(i,j,k)+curlez(i,j,k))
             enddo
           enddo
         enddo
 
 
!************************
!  end of iteration loop
!************************

  10     continue
         CALL XREALBCC_PACK_B(BX,BY,BZ,1,NX,NY,NZ)
 
 
!  set ghost cell values for B: it is only for diagnostics as these B
!  values are not used anywhere
 

      do k=kb-1,ke+1
        do j = jb-1,je+1
          bx(1  ,j,k)=bx(2  ,j,k)
          by(1  ,j,k)=by(2  ,j,k)
          bz(1  ,j,k)=bz(2  ,j,k)
        enddo
      enddo
 

       call date_and_time(values=time_end_array(:,22))
       call accumulate_time_difference(time_begin_array(1,22) &
     &                                ,time_end_array(1,22) &
     &                                ,time_elapsed(22))
 
 
       return
       end
!
!################################################################################
!
      subroutine focalc

      use parameter_mod
      use MESH2D
      double precision:: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &
                        ,curr_tot

      do k=kb,ke 
        do j = jb,je
          do i=2,nx1
            dbxdy= bx(i,j  ,k)+bx(i-1,j  ,k)&
                 +bx(i-1,j  ,k-1)+bx(i,j  ,k-1)&
                 -bx(i,j-1,k)-bx(i-1,j-1,k)&
                 -bx(i-1,j-1,k-1)-bx(i,j-1,k-1)
            dbxdz= bx(i,j,k  )+bx(i-1,j,k  )&
                 +bx(i-1,j-1,k  )+bx(i,j-1,k  )&
                 -bx(i,j,k-1)-bx(i-1,j,k-1)&
                 -bx(i-1,j-1,k-1)-bx(i,j-1,k-1)
            dbydx= by(i  ,j,k)+by(i  ,j-1,k)&
                 +by(i  ,j-1,k-1)+by(i  ,j,k-1)&
                 -by(i-1,j,k)-by(i-1,j-1,k)&
                 -by(i-1,j-1,k-1)-by(i-1,j,k-1)
            dbydz= by(i,j,k  )+by(i-1,j,k  )&
                 +by(i-1,j-1,k  )+by(i,j-1,k  )&
                 -by(i,j,k-1)-by(i-1,j,k-1)&
                 -by(i-1,j-1,k-1)-by(i,j-1,k-1)
            dbzdx= bz(i  ,j,k)+bz(i  ,j-1,k)&
                 +bz(i  ,j-1,k-1)+bz(i  ,j,k-1)&
                 -bz(i-1,j,k)-bz(i-1,j-1,k)&
                 -bz(i-1,j-1,k-1)-bz(i-1,j,k-1)
            dbzdy= bz(i,j  ,k)+bz(i-1,j  ,k)&
                 +bz(i-1,j  ,k-1)+bz(i,j  ,k-1)&
                 -bz(i,j-1,k)-bz(i-1,j-1,k)&
                 -bz(i-1,j-1,k-1)-bz(i,j-1,k-1)

!           Uniform mesh - Same as is in version 5.0
!            curlbx_scalar=dbzdy/(4.*hy)-dbydz/(4.*hz)
!            curlby_scalar=dbxdz/(4.*hz)-dbzdx/(4.*hx)
!            curlbz_scalar=dbydx/(4.*hx)-dbxdy/(4.*hy)


!          Nonuniform mesh
            curlbx_scalar=dbzdy/(4.*meshY%dxc(j+1))-dbydz/(4.*meshZ%dxc(k+1))
            curlby_scalar=dbxdz/(4.*meshZ%dxc(k+1))-dbzdx/(4.*meshX%dxc(i  ))
            curlbz_scalar=dbydx/(4.*meshX%dxc(i  ))-dbxdy/(4.*meshY%dxc(j+1))

!
! 6/25/2006 New eta_par option: tensor eta
!
            bx1=bx(i+1,j+1,k)  +bdipole_x(i+1,j+1,k)
            bx2=bx(i  ,j+1,k)  +bdipole_x(i  ,j+1,k)
            bx3=bx(i  ,j  ,k)  +bdipole_x(i  ,j  ,k)
            bx4=bx(i+1,j  ,k)  +bdipole_x(i+1,j  ,k)
            bx5=bx(i+1,j+1,k+1)+bdipole_x(i+1,j+1,k+1)
            bx6=bx(i  ,j+1,k+1)+bdipole_x(i  ,j+1,k+1)
            bx7=bx(i  ,j  ,k+1)+bdipole_x(i  ,j  ,k+1)
            bx8=bx(i+1,j  ,k+1)+bdipole_x(i+1,j  ,k+1)
            by1=by(i+1,j+1,k)  +bdipole_y(i+1,j+1,k)
            by2=by(i  ,j+1,k)  +bdipole_y(i  ,j+1,k)
            by3=by(i  ,j  ,k)  +bdipole_y(i  ,j  ,k)
            by4=by(i+1,j  ,k)  +bdipole_y(i+1,j  ,k)
            by5=by(i+1,j+1,k+1)+bdipole_y(i+1,j+1,k+1)
            by6=by(i  ,j+1,k+1)+bdipole_y(i  ,j+1,k+1)
            by7=by(i  ,j  ,k+1)+bdipole_y(i  ,j  ,k+1)
            by8=by(i+1,j  ,k+1)+bdipole_y(i+1,j  ,k+1)
            bz1=bz(i+1,j+1,k)  +bdipole_z(i+1,j+1,k)
            bz2=bz(i  ,j+1,k)  +bdipole_z(i  ,j+1,k)
            bz3=bz(i  ,j  ,k)  +bdipole_z(i  ,j  ,k)
            bz4=bz(i+1,j  ,k)  +bdipole_z(i+1,j  ,k)
            bz5=bz(i+1,j+1,k+1)+bdipole_z(i+1,j+1,k+1)
            bz6=bz(i  ,j+1,k+1)+bdipole_z(i  ,j+1,k+1)
            bz7=bz(i  ,j  ,k+1)+bdipole_z(i  ,j  ,k+1)
            bz8=bz(i+1,j  ,k+1)+bdipole_z(i+1,j  ,k+1)
            bxav=.125*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar
		    tenx=eta(i,j,k)*xj
		    teny=eta(i,j,k)*yj
		    tenz=eta(i,j,k)*zj
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
      call XREALBCC_PACK_E(fox,foy,foz,1,NX,NY,NZ)


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
          if (jb == 1) then
            fox(i,jb-1,k)=fox(i,jb,k)
            foy(i,jb-1,k)=foy(i,jb,k)
            foz(i,jb-1,k)=foz(i,jb,k)
          endif
          if (je == ny) then
            fox(i,je+1,k)=fox(i,je,k)
            foy(i,je+1,k)=foy(i,je,k)
            foz(i,je+1,k)=foz(i,je,k)
          endif
        enddo
      enddo

      if (kb == 1) then
        do j = jb-1,je+1
          do i=1,nx2
            fox(i,j,kb-1)=fox(i,j,kb)
            foy(i,j,kb-1)=foy(i,j,kb)
            foz(i,j,kb-1)=foz(i,j,kb)
          enddo
        enddo
      endif

      if (ke == nz) then
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
!!#######################################################################
!
      subroutine nsmth (a)
 
      use parameter_mod
      double precision, dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: temp, a

!  smoothing routine--assumes aperiodic in x
 
      call XREALBCC(a,0,NX,NY,NZ)
      do k=kb-1,ke+1
        do j = jb-1,je+1
          do i=1,nx2
            temp(i,j,k)=a(i,j,k)
          enddo
        enddo
      enddo

      do k=kb,ke
        do j = jb,je
          do i=2,nx1
            a(i,j,k)=temp(i,j,k)/8.&
      +( temp(i-1,j,k)+temp(i+1,j,k)+temp(i,j+1,k)+temp(i,j-1,k)&
      +temp(i,j,k+1)+temp(i,j,k-1))/16.&
      +( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k)&
      +temp(i-1,j-1,k)&
      +temp(i,j+1,k+1)+temp(i,j-1,k+1)+temp(i,j+1,k-1)+temp(i,j-1,k-1)&
      +temp(i+1,j,k+1)+temp(i-1,j,k+1)+temp(i+1,j,k-1)&
      +temp(i-1,j,k-1))/32.&
      +( temp(i+1,j+1,k+1)+temp(i-1,j+1,k+1)&
      +temp(i+1,j-1,k+1)+temp(i-1,j-1,k+1)&
      +temp(i+1,j+1,k-1)+temp(i-1,j+1,k-1)&
      +temp(i+1,j-1,k-1)+temp(i-1,j-1,k-1))/64.
          enddo
        enddo
      enddo
 
      if (jb == 1) then
        do k=kb,ke
          do i=2,nx1
            a(i,jb-1,k)=a(i,jb,k)
          enddo
        enddo
      endif
      if (je == ny) then
        do k=kb,ke
          do i=2,nx1
            a(i,je+1,k)=a(i,je,k)
          enddo
        enddo
      endif

      if (kb == 1) then
        do j = jb-1,je+1
          do i=2,nx1
            a(i,j,kb-1)=a(i,j,kb)
          enddo
        enddo
      endif

      if (ke == nz) then
        do j = jb-1,je+1
          do i=2,nx1
            a(i,j,ke+1)=a(i,j,ke)
          enddo
        enddo
      endif

      do k=kb,ke
        do j = jb-1,je+1
          a(1  ,j,k)=a(2  ,j,k)
          a(nx2,j,k)=a(nx1,j,k)
        enddo
      enddo

      call XREALBCC(a,0,NX,NY,NZ)

      return
      end
!
!#######################################################################
!
!
      double precision FUNCTION myranf()
 
      implicit none
      INTEGER*8:: idum
      INTEGER*8:: MBIG,MSEED,MZ
      double precision:: FAC
      INTEGER*8:: i,iff,ii,inext,inextp,k
      INTEGER*8:: mj,mk,ma(55)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0 &
                ,FAC=1.d-09)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      common /myrandom/ idum
 
      if (idum < 0.or.iff == 0)then
        iff=1
        mj=MSEED-iabs(int(idum,4))
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if (mk < MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if (ma(i) < MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if (inext == 56)inext=1
      inextp=inextp+1
      if (inextp == 56)inextp=1
      mj=ma(inext)-ma(inextp)
      if (mj < MZ)mj=mj+MBIG
      ma(inext)=mj
      myranf=dble(mj*FAC)
 
      end function myranf
!
!***********************************************************************
!
      subroutine wrtdatum(ndatum,datum,f_unit)
 
      use parameter_mod
      implicit none
      integer*8:: ndatum,f_unit
      integer:: ndatum_4
      double precision:: datum(ndatum),datum_tmp(ndatum)
 
      ndatum_4=ndatum
      do ipe=0,numprocs-1
        if (myid == 0) then
          if (ipe == 0) then
            write(f_unit) datum
          else
              call MPI_IRECV(datum_tmp,ndatum_4,MPI_DOUBLE_PRECISION&
                            ,ipe,0,MPI_COMM_WORLD,req(1),IERR)
              CALL MPI_WAIT(REQ(1),STATUS1,IERR)
              write(f_unit) datum_tmp
          endif
        else
          if (myid == ipe) then
            call MPI_ISEND(datum,ndatum_4,MPI_DOUBLE_PRECISION&
                          ,0,0,MPI_COMM_WORLD,req(2),IERR)
            CALL MPI_WAIT(REQ(2),STATUS2,IERR)
          endif
        endif
      enddo
 
      return
      end
!
!***********************************************************************
!
      subroutine readdatum(ndatum,datum,f_unit)
 
 
      use parameter_mod
      implicit none
      integer*8:: ndatum,f_unit
      integer:: ndatum_4
      double precision:: datum(ndatum),datum_tmp(ndatum)
 
 
      ndatum_4=ndatum
      do ipe=0,numprocs-1
        if (myid == 0) then
          if (ipe == 0) then
            read(f_unit) datum
          else
              read(f_unit) datum_tmp
              call MPI_ISEND(datum_tmp,ndatum_4,MPI_DOUBLE_PRECISION&
                            ,ipe,0,MPI_COMM_WORLD,req(1),IERR)
              CALL MPI_WAIT(REQ(1),STATUS1,IERR)
          endif
        else
          if (myid == ipe) then
            call MPI_IRECV(datum,ndatum_4,MPI_DOUBLE_PRECISION&
                          ,0,0,MPI_COMM_WORLD,req(2),IERR)
            CALL MPI_WAIT(REQ(2),STATUS2,IERR)
          endif
        endif
      enddo
 
      return
      end
!
!***********************************************************************
!
      subroutine accumulate_time_difference(time_begin,time_end &
     &                                     ,time_elapsed)
 
      implicit none
      integer,dimension(8):: time_begin,time_end
      double precision:: time_elapsed
 
      time_elapsed=time_elapsed &
                   +(time_end(3)-time_begin(3))*3600.*24. &
                   +(time_end(5)-time_begin(5))*3600. &
                   +(time_end(6)-time_begin(6))*60. &
                   +(time_end(7)-time_begin(7)) &
                   +(time_end(8)-time_begin(8))*0.001
 
      return
      end
!
! HXV (non-Yee mesh), global setup
!
!***********************************************************************
!
      subroutine init_global
!
!=======================================================================
!
!     Initialize particle position and velocity: assumed B is
!     in a 2-D plane
!
!=======================================================================
!
      use parameter_mod
      use MESH2D
      integer*8:: random_init,ibp1,ibp2,nptot_max,i,remake
      double precision:: my_random_init,myranf,by_avg,bz_avg,q_p
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi,profile_factor   &
                        ,b_dipole_sqr,x_p,y_p,z_p,x_p_logical,y_p_logical        &
                        ,z_p_logical,r_c
 

      remake = 0

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

 
      random_init=5000000*myid+1000
      IDUM=0
      do i=1,random_init
        my_random_init=myranf()
      enddo
 
!  initialize the random number generator with a different
!  seed for different processors
 
      ranval=myranf()
      pi=acos(-1.d+00)
      pifac=180./pi
      it=0
      itfin = 0
      nx1=nx+1
      nx2=nx+2
      ny1=ny+1
      ny2=ny+2
      nz1=nz+1
      nz2=nz+2
      hx=xmax/nx
      hy=ymax/ny
      hz=zmax/nz
      hxi=1./hx
      hyi=1./hy
      hzi=1./hz
 
 
 
!     Uniform mesh - Same as is in version 5.0
!      x_dipole=hx*(i_dipole+0.5)
!      y_dipole=hy*(j_dipole+0.5)
!      z_dipole=hz*(k_dipole+0.5)
 
!     Nonuniform mesh
      x_dipole=meshX%xc(i_dipole)
      y_dipole=meshY%xc(j_dipole)
      z_dipole=meshZ%xc(k_dipole)

 
!  zbglobal(ipe) gives the starting position, in z, for processor ipe
!  zeglobal(ipe) gives the ending position, in z, for ipe
 
!     Uniform mesh - Same as is in version 5.0
!      zb=(kb-1)*hz
!      ze= ke   *hz
!      do ipe=0,npes-1
!        zbglobal(ipe)=(kbglobal(ipe)-1)*hz
!        zeglobal(ipe)= keglobal(ipe)   *hz
!      enddo
!      yb=(jb-1)*hy
!      ye= je   *hy
!      do ipe=0,npes-1
!        ybglobal(ipe)=(jbglobal(ipe)-1)*hy
!        yeglobal(ipe)= jeglobal(ipe)   *hy
!      enddo

!     Nonuniform mesh
        zb=meshZ%xn(kb+1)
        ze=meshZ%xn(ke+2)
        do ipe=0,npes-1
          zbglobal(ipe)=meshZ%xn(kbglobal(ipe)+1)
          zeglobal(ipe)=meshZ%xn(keglobal(ipe)+2)
        enddo
        yb=meshY%xn(jb+1)
        ye=meshY%xn(je+2)
        do ipe=0,npes-1
          ybglobal(ipe)=meshY%xn(jbglobal(ipe)+1)
          yeglobal(ipe)=meshY%xn(jeglobal(ipe)+2)
        enddo
        volume_fraction = (ye-yb)*(ze-zb)/(ymax*zmax)

        xb        = 0.
        xe        = xmax
        xb_logical=MESH_UNMAP(meshX,xb)
        xe_logical=MESH_UNMAP(meshX,xe)
        yb_logical=MESH_UNMAP(meshY,yb)
        ye_logical=MESH_UNMAP(meshY,ye)
        zb_logical=MESH_UNMAP(meshZ,zb)
        ze_logical=MESH_UNMAP(meshZ,ze)


        do is=1,nspec
          npm=npx(is)*npy(is)*npz(is)*npes
          dfac(is)=real(ny*nz*nx)/real(npm)
          do ixe=1,nx2 
            do iye=jb-1,je+1
              do ize=kb-1,ke+1
!                 qp_cell(ixe,iye,ize,is) = (meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)/(hx*hy*hz))*dfac(is)*frac(is)
                 qp_cell(ixe,iye,ize,is) = meshX%dxc(ixe)*meshY%dxc(iye+1)*meshZ%dxc(ize+1)*dfac(is)*frac(is)
              enddo
            enddo
          enddo
        enddo


!

      if (uniform_loading_in_logical_grid) then
        nptotp=0
        do is=1,nspec
          nptotp = nptotp + npx(is)*npy(is)*npz(is)
        enddo
      else
        nptotp=0
        do is=1,nspec
          nptotp = nptotp + npx(is)*npy(is)*npz(is)*npes*volume_fraction
        enddo
      endif

      call MPI_ALLREDUCE(nptotp,nptot_max,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,IERR)
      if (nptot_max > nplmax) then

          remake = 1
          nplmax = 2*nptot_max

          deallocate (x)
          allocate   (x(nplmax))

          deallocate (y)
          allocate   (y(nplmax))

          deallocate (z)
          allocate   (z(nplmax))

          deallocate (vx)
          allocate   (vx(nplmax))

          deallocate (vy)
          allocate   (vy(nplmax))

          deallocate (vz)
          allocate   (vz(nplmax))

          deallocate (qp)
          allocate   (qp(nplmax))

          deallocate (link)
          allocate   (link(nplmax))

          deallocate (porder)
          allocate   (porder(nplmax))

      endif

      if (ndim /=1) then
        call dipole_field
      else
        call dipole_field_2d
      endif

      if (remake /= 0) call makelist

!        write(6,*) "myid = ",myid," FINISHED MAKELIST"
!        call MPI_FINALIZE(IERR)
!        STOP
 
      curion=0.5
      phibr=phib/pifac
      phicr = 0.5 *(180.0-phib)/ pifac
      do is=1,nspec
        ninj(is)=0
        ninj_global(is)=0
        npart(is)=0
        tx0(is)=btspec(is)/(2.*wpiwci**2)
        x0(is)=0.
        x1(is)=xmax
        if (is == 1) x1(is)=fxsho*xmax
        if (is == 2) x0(is)=fxsho*xmax
        call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      enddo
      te0=bete/(2.*wpiwci**2)
      vbal=1.

      if (b_in_xy) then
        bxc = cos(theta/pifac)/wpiwci
        byc = sin(theta/pifac)/wpiwci
        bzc = 0.
      else
        bxc = cos(theta/pifac)/wpiwci
        byc = 0.
        bzc = sin(theta/pifac)/wpiwci
      endif

! define cell positions
 
      nxa=nx*fxsho+1
      if (nxcel  ==  0)  then
        nxcelb = 1
      else
        nxcelb = nxcel
      endif
      nxb=nxa-nxcelb
      nxc=nxa+nxcelb
 
 
!      CALL MPI_FINALIZE(IERR)
!      STOP
 
      do is = 1, nspec
        if (is == 2) then
          iseed(1) = 2*myid + 8
          call random_seed(put = iseed(1:k))
          ranval=myranf()
        endif
        ipb1 = 1

!       Uniform mesh - Same as is in ver5.0
!        ipb2 = npx(is)*npy(is)*npz(is)

!       Nonuniform mesh
        if (uniform_loading_in_logical_grid) then
          ipb2 = npx(is)*npy(is)*npz(is)
        else
          ipb2 = npx(is)*npy(is)*npz(is)*npes*volume_fraction
        endif

        isp=is
        nk=0

        npm=npx(is)*npy(is)*npz(is)*npes
        dfac(is)=real(ny*nz)*(x1(is)-x0(is))/(hx*real(npm))
        vpar(is)=sqrt(btspec(is)/(wspec(is)*wpiwci*wpiwci))
        vper(is)=vpar(is)*sqrt(anisot(is))

        do ip=ipb1,ipb2
            
          if (ipb2 == ipb1) write(6,*) "myid = , # particles = ",myid,ipb1,ipb2

          if (uniform_loading_in_logical_grid) then
            X_P_LOGICAL  = XB_LOGICAL+(XE_LOGICAL-XB_LOGICAL)*myranf()
            Y_P_LOGICAL  = YB_LOGICAL+(YE_LOGICAL-YB_LOGICAL)*myranf()
            Z_P_LOGICAL  = ZB_LOGICAL+(ZE_LOGICAL-ZB_LOGICAL)*myranf()
            X_P          = MESH_MAP(meshX,X_P_LOGICAL)
            Y_P          = MESH_MAP(meshY,Y_P_LOGICAL)
            Z_P          = MESH_MAP(meshZ,Z_P_LOGICAL)
            IXE          = dtxi*X_P_LOGICAL+1.50000000000d+00
            IYE          = dtyi*Y_P_LOGICAL+1.50000000000d+00
            IZE          = dtzi*Z_P_LOGICAL+1.50000000000d+00
!            q_p          = (meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)/(hx*hy*hz))*dfac(is)*frac(is)
            q_p          = meshX%dxc(ixe)*meshY%dxc(iye)*meshZ%dxc(ize)*dfac(is)*frac(is)
          else
            X_P  = X0(IS)+(X1(IS)-X0(IS))*MYRANF()
            Y_P  = YB+(YE-YB)*myranf()   
            Z_P  = ZB+(ZE-ZB)*myranf()
            q_p  = hx*hy*hz*dfac(is)*frac(is)
          endif
         
          if (q_p == 0) write(6,*) "q_p = ",q_p,ixe,iye,ize
 
          if (ndim /=1) then
            r_particle=sqrt((x_p-x_dipole)**2+(y_p-y_dipole)**2+(z_p-z_dipole)**2)
          else
            r_particle=sqrt((x_p-x_dipole)**2+(y_p-y_dipole)**2)
          endif

          if (r_particle < dipole_sphere_radius) goto 10
 
          np     = ipstore
          x(np)  = x_p
          y(np)  = y_p
          z(np)  = z_p
          qp(np) = q_p


!         Nonuniform mesh - using MESH_UNMAP
          rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000d+00
          rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000d+00
          rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000d+00
          ixe=rxe
          iye=rye
          ize=rze
          iye=iye-1             ! integer index in y direction starts at 0
          ize=ize-1             ! integer index in z direction starts at 0
!

 
! Modify ion temperature according to B profile for pressure balance
! (this is nearly but not quite complete for a jump in quantities
!  across the layer)


!
          if (profile_power /= 0) then
            if (meshX%xc(ixe) >= x_dipole) then
              profile_factor = 1.
            else

              if (ndim /= 1) then
                r_c = sqrt((meshX%xc(ixe)-x_dipole)**2+(meshY%xc(iye+1)-y_dipole)**2+(meshZ%xc(ize+1)-z_dipole)**2)
              else
                r_c = sqrt((meshX%xc(ixe)-x_dipole)**2+(meshY%xc(iye+1)-y_dipole)**2)
              endif

              if (r_c > dipole_sphere_radius+moat_zone) then
                b_dipole_sqr=(bdipole_x(ixe,iye+1,ize+1)**2+bdipole_y(ixe+1,iye+1,ize+1)**2+bdipole_z(ixe+1,iye+1,ize+1)**2)    &
                             *wpiwci**2
                profile_factor = max(denmin,min(1.,(2./sqrt(b_dipole_sqr))**profile_power))
              else
                profile_factor = 1.
              endif

            endif

          else

            profile_factor = 1.

          endif
          qp(np) = qp(np)*profile_factor
!
 
          vthfac = 1./sqrt(profile_factor)
          vcurr = 0.0
 
 
          xran=myranf()
          vmag=sqrt(-alog(1.-.999999*xran))
          xran=myranf()
          th=2.*pi*xran
          vxa=vpar(is)*vmag*cos(th)*vthfac
 

          xran=myranf()
          vmag=sqrt(-alog(1.-.999999*xran))
          xran=myranf()
          th=2.*pi*xran
          vya=vper(is)*vmag*cos(th)*vthfac
          vza=vper(is)*vmag*sin(th)*vthfac
 
 
          vx(np)=vxa+vxbar(is)
          vy(np)=vya
          vz(np)=vza

          ipstore=link(np)
          link(np)=iphead(ixe,iye,ize,is)
          iphead(ixe,iye,ize,is)=np
 
 10       CONTINUE
 
 
        enddo
      enddo
 
!      CALL MPI_FINALIZE(IERR)
!      STOP
 
      bx=bxc*dipole_sphere_ex;by=byc*dipole_sphere_ey;bz=bzc*dipole_sphere_ez
      ex=0.;ey=0.;ez=0.
 
 
      exc= vzbar(1)*byc-vybar(1)*(bzc+bz_IMF/wpiwci)
      eyc=+vxbar(1)*(bzc+bz_IMF/wpiwci)-vzbar(1)*bxc
      ezc= vybar(1)*bxc-vxbar(1)*byc

      if (myid == 0) then
        write(6,*) " bxc = ",bxc*wpiwci
        write(6,*) " byc = ",byc*wpiwci
        write(6,*) " bzc = ",bzc*wpiwci
      endif

      if (ndim /=1) then
        do k=kb-1,ke  
          do j=jb-1,je  
            do i=1,nx1
               if (i/=1) then
                bz_avg= 0.125*( bz       (i+1,j+1,k  ) &
                               +bz       (i  ,j+1,k  ) &
                               +bz       (i+1,j  ,k  ) &
                               +bz       (i  ,j  ,k  ) &
                               +bz       (i+1,j+1,k+1) &
                               +bz       (i  ,j+1,k+1) &
                               +bz       (i+1,j  ,k+1) &
                               +bz       (i  ,j  ,k+1) &
                              )                        &
                       +0.125*( bdipole_z(i+1,j+1,k  ) &
                               +bdipole_z(i  ,j+1,k  ) &
                               +bdipole_z(i+1,j  ,k  ) &
                               +bdipole_z(i  ,j  ,k  ) &
                               +bdipole_z(i+1,j+1,k+1) &
                               +bdipole_z(i  ,j+1,k+1) &
                               +bdipole_z(i+1,j  ,k+1) &
                               +bdipole_z(i  ,j  ,k+1) &
                              )
                by_avg= 0.125*( by       (i+1,j+1,k  ) &
                               +by       (i  ,j+1,k  ) &
                               +by       (i+1,j  ,k  ) &
                               +by       (i  ,j  ,k  ) &
                               +by       (i+1,j+1,k+1) &
                               +by       (i  ,j+1,k+1) &
                               +by       (i+1,j  ,k+1) &
                               +by       (i  ,j  ,k+1) &
                              )                        &
                       +0.125*( bdipole_y(i+1,j+1,k  ) &
                               +bdipole_y(i  ,j+1,k  ) &
                               +bdipole_y(i+1,j  ,k  ) &
                               +bdipole_y(i  ,j  ,k  ) &
                               +bdipole_y(i+1,j+1,k+1) &
                               +bdipole_y(i  ,j+1,k+1) &
                               +bdipole_y(i+1,j  ,k+1) &
                               +bdipole_y(i  ,j  ,k+1) &
                              )
               else
                  bx_avg= bxc
                  by_avg= byc
                  bz_avg= bzc+bz_IMF/wpiwci
               endif
                ex(i,j,k)= vzbar(1)*by_avg-vybar(1)*bz_avg
                ey(i,j,k)=+vxbar(1)*bz_avg-vzbar(1)*bx_avg
                ez(i,j,k)= vybar(1)*bx_avg-vxbar(1)*by_avg
            enddo
          enddo
        enddo
      else
        do k=kb-1,ke  
          do j=jb-1,je  
            do i=1,nx1
               if (i/=1) then
                 bz_avg= 0.250*( bz       (i+1,j+1,k+1) &
                                +bz       (i  ,j+1,k+1) &
                                +bz       (i+1,j  ,k+1) &
                                +bz       (i  ,j  ,k+1) &
                               )                        &
                        +0.250*( bdipole_z(i+1,j+1,k+1) &
                                +bdipole_z(i  ,j+1,k+1) &
                                +bdipole_z(i+1,j  ,k+1) &
                                +bdipole_z(i  ,j  ,k+1) &
                               )                        
                 by_avg= 0.250*( by       (i+1,j+1,k+1) &
                                +by       (i  ,j+1,k+1) &
                                +by       (i+1,j  ,k+1) &
                                +by       (i  ,j  ,k+1) &
                               )                        &
                        +0.250*( bdipole_y(i+1,j+1,k+1) &
                                +bdipole_y(i  ,j+1,k+1) &
                                +bdipole_y(i+1,j  ,k+1) &
                                +bdipole_y(i  ,j  ,k+1) &
                               )
                else
                  bx_avg= bxc
                  by_avg= byc
                  bz_avg= bzc+bz_IMF/wpiwci
                endif
                ex(i,j,k)= vzbar(1)*by_avg-vybar(1)*bz_avg
                ey(i,j,k)=+vxbar(1)*bz_avg-vzbar(1)*bx_avg

                ez(i,j,k)= vybar(1)*bx_avg-vxbar(1)*by_avg
            enddo
          enddo
        enddo
      endif
 
 
      ey(nx2,:,:)=ey(nx1,:,:)
      ez(nx2,:,:)=ez(nx1,:,:)

      if (ndim /= 1) then
        call xrealbcc(ey,1,nx,ny,nz)
        call xrealbcc(ez,1,nx,ny,nz)
      else
        call xrealbcc_pack_e_2d(ex,ey,ez,1,nx,ny,nz)
      endif
 
 
      ex=ex*dipole_sphere_ex
      ey=ey*dipole_sphere_ey
      ez=ez*dipole_sphere_ez
 
!***********************
!   set the friction force and resitivity
!   friction force can be zero at t=0
!**********************
 
      do k=kb-1,ke+1
        do j=jb-1,je+1
          do i=1,nx2
            fox(i,j,k)=0.
            foy(i,j,k)=0.
            foz(i,j,k)=0.
            eta(i,j,k)=resis
          enddo
        enddo
      enddo
      dtsav=dt
      dt=0.
      call trans
 

      dt=dtsav
      if (myid == 0) write(6,*) " BEFORE FIELD (CALLED BY INIT)"
      if (.not.testorbt) then
        do field_subcycle=1,n_subcycles
          if (ndim /= 1) then
            call field
          else
            call field_2d
          endif
        enddo
      endif
 18      CONTINUE

      if (myid == 0) write(6,*) " AFTER  FIELD (CALLED BY INIT)"
         if (ndim /= 1) then
           call etacalc      ! Dietmar's resistivity
         else
           call etacalc_2d   ! Dietmar's resistivity
         endif
 
 999  CONTINUE
 
      return
      end
!
!***********************************************************************
!
      subroutine caltemp2_global
 
      use parameter_mod
      use MESH2D
      double precision:: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,xx,xy,xz,yy,yz,zz
      integer*8 ix,iy,iz,ixp1,iyp1,izp1
 
      call date_and_time(values=time_begin_array(:,23))

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt


      tpar=0.
      tperp=0.

      p_xx=0.;p_xy=0.;p_xz=0.;p_yy=0.;p_yz=0.;p_zz=0.

      if (nspec >= 2) then
         rfrac = frac(2)/frac(1)
      else
         rfrac = 0.
      endif
 
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

!  begin advance of particle position and velocity
!  If dt=0, skip
!
              DO WHILE (NP.NE.0)
                L=NP

!               Uniform mesh - Same as is in version 5.0
!                rx=hxi*x(l)+1.5000000000000001
!                ry=hyi*y(l)+0.5000000000000001d+00
!                rz=hzi*z(l)+0.5000000000000001d+00
!                ix=rx
!                iy=ry
!                iz=rz
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

!               Nonuniform mesh - using MESH_UNMAP
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

                w1=(1.-fx)*(1.-fy)*(1.-fz)
                w2=    fx *(1.-fy)*(1.-fz)
                w3=(1.-fx)*    fy *(1.-fz)
                w4=    fx*     fy *(1.-fz)
                w5=(1.-fx)*(1.-fy)*    fz
                w6=    fx *(1.-fy)*    fz
                w7=(1.-fx)*    fy*     fz
                w8=    fx*     fy*     fz

                dns1= dns(ix  ,iy  ,iz  ,1)*w1+dns(ixp1,iy  ,iz  ,1)*w2  &
                +     dns(ix  ,iyp1,iz  ,1)*w3+dns(ixp1,iyp1,iz  ,1)*w4  &
                +     dns(ix  ,iy  ,izp1,1)*w5+dns(ixp1,iy  ,izp1,1)*w6  &
                +     dns(ix  ,iyp1,izp1,1)*w7+dns(ixp1,iyp1,izp1,1)*w8

                dns2= 0.

                dnst = dns1 + dns2
      
                vxavg1=vxs(ix  ,iy  ,iz  ,1)*w1+vxs(ixp1,iy  ,iz  ,1)*w2  &
                +      vxs(ix  ,iyp1,iz  ,1)*w3+vxs(ixp1,iyp1,iz  ,1)*w4  &
                +      vxs(ix  ,iy  ,izp1,1)*w5+vxs(ixp1,iy  ,izp1,1)*w6  &
                +      vxs(ix  ,iyp1,izp1,1)*w7+vxs(ixp1,iyp1,izp1,1)*w8

                vxavg2= 0.

                vxavg = (dns1*vxavg1 + dns2*vxavg2)/dnst

                vyavg1=vys(ix  ,iy  ,iz  ,1)*w1+vys(ixp1,iy  ,iz  ,1)*w2  &
                +      vys(ix  ,iyp1,iz  ,1)*w3+vys(ixp1,iyp1,iz  ,1)*w4  &
                +      vys(ix  ,iy  ,izp1,1)*w5+vys(ixp1,iy  ,izp1,1)*w6  &
                +      vys(ix  ,iyp1,izp1,1)*w7+vys(ixp1,iyp1,izp1,1)*w8

                vyavg2=0.
  
                vyavg = (dns1*vyavg1 + dns2*vyavg2)/dnst

                vzavg1=vzs(ix  ,iy  ,iz  ,1)*w1+vzs(ixp1,iy  ,iz  ,1)*w2  &
                +      vzs(ix  ,iyp1,iz  ,1)*w3+vzs(ixp1,iyp1,iz  ,1)*w4  &
                +      vzs(ix  ,iy  ,izp1,1)*w5+vzs(ixp1,iy  ,izp1,1)*w6  &
                +      vzs(ix  ,iyp1,izp1,1)*w7+vzs(ixp1,iyp1,izp1,1)*w8

                vzavg2=0.

                vzavg = (dns1*vzavg1 + dns2*vzavg2)/dnst

                vxa=vx(l)-vxavg
                vya=vy(l)-vyavg
                vza=vz(l)-vzavg

                bxa  =bx(ix  ,iy  ,iz  )*w1+bx(ixp1,iy  ,iz  )*w2  &
                +     bx(ix  ,iyp1,iz  )*w3+bx(ixp1,iyp1,iz  )*w4  &
                +     bx(ix  ,iy  ,izp1)*w5+bx(ixp1,iy  ,izp1)*w6  &
                +     bx(ix  ,iyp1,izp1)*w7+bx(ixp1,iyp1,izp1)*w8

                bya  =by(ix  ,iy  ,iz  )*w1+by(ixp1,iy  ,iz  )*w2  &
                +     by(ix  ,iyp1,iz  )*w3+by(ixp1,iyp1,iz  )*w4  &
                +     by(ix  ,iy  ,izp1)*w5+by(ixp1,iy  ,izp1)*w6  &
                +     by(ix  ,iyp1,izp1)*w7+by(ixp1,iyp1,izp1)*w8

                bza  =bz(ix  ,iy  ,iz  )*w1+bz(ixp1,iy  ,iz  )*w2  &
                +     bz(ix  ,iyp1,iz  )*w3+bz(ixp1,iyp1,iz  )*w4  &
                +     bz(ix  ,iy  ,izp1)*w5+bz(ixp1,iy  ,izp1)*w6  &
                +     bz(ix  ,iyp1,izp1)*w7+bz(ixp1,iyp1,izp1)*w8

                btota=sqrt(bxa**2+bya**2+bza**2)
                if (btota < 1.e-20) btota=1.e-20
                wpar=(vxa*bxa+vya*bya+vza*bza)/btota
                wperp2=vxa**2+vya**2+vza**2-wpar**2
                xx=vxa*vxa
                xy=vxa*vya
                xz=vxa*vza
                yy=vya*vya
                yz=vya*vza
                zz=vza*vza

                tpar (ix  ,iy  ,iz  ,is)=tpar (ix  ,iy  ,iz  ,is)+qp(np)*w1*wpar*wpar
                tpar (ixp1,iy  ,iz  ,is)=tpar (ixp1,iy  ,iz  ,is)+qp(np)*w2*wpar*wpar 
                tpar (ix  ,iyp1,iz  ,is)=tpar (ix  ,iyp1,iz  ,is)+qp(np)*w3*wpar*wpar 
                tpar (ixp1,iyp1,iz  ,is)=tpar (ixp1,iyp1,iz  ,is)+qp(np)*w4*wpar*wpar 
                tpar (ix  ,iy  ,izp1,is)=tpar (ix  ,iy  ,izp1,is)+qp(np)*w5*wpar*wpar 
                tpar (ixp1,iy  ,izp1,is)=tpar (ixp1,iy  ,izp1,is)+qp(np)*w6*wpar*wpar 
                tpar (ix  ,iyp1,izp1,is)=tpar (ix  ,iyp1,izp1,is)+qp(np)*w7*wpar*wpar 
                tpar (ixp1,iyp1,izp1,is)=tpar (ixp1,iyp1,izp1,is)+qp(np)*w8*wpar*wpar 
                tperp(ix  ,iy  ,iz  ,is)=tperp(ix  ,iy  ,iz  ,is)+qp(np)*w1*wperp2 
                tperp(ixp1,iy  ,iz  ,is)=tperp(ixp1,iy  ,iz  ,is)+qp(np)*w2*wperp2 
                tperp(ix  ,iyp1,iz  ,is)=tperp(ix  ,iyp1,iz  ,is)+qp(np)*w3*wperp2 
                tperp(ixp1,iyp1,iz  ,is)=tperp(ixp1,iyp1,iz  ,is)+qp(np)*w4*wperp2
                tperp(ix  ,iy  ,izp1,is)=tperp(ix  ,iy  ,izp1,is)+qp(np)*w5*wperp2 
                tperp(ixp1,iy  ,izp1,is)=tperp(ixp1,iy  ,izp1,is)+qp(np)*w6*wperp2
                tperp(ix  ,iyp1,izp1,is)=tperp(ix  ,iyp1,izp1,is)+qp(np)*w7*wperp2 
                tperp(ixp1,iyp1,izp1,is)=tperp(ixp1,iyp1,izp1,is)+qp(np)*w8*wperp2
                dpedx(ix  ,iy  ,iz  )=dpedx(ix  ,iy  ,iz  )+qp(np)*w1
                dpedx(ixp1,iy  ,iz  )=dpedx(ixp1,iy  ,iz  )+qp(np)*w2 
                dpedx(ix  ,iyp1,iz  )=dpedx(ix  ,iyp1,iz  )+qp(np)*w3 
                dpedx(ixp1,iyp1,iz  )=dpedx(ixp1,iyp1,iz  )+qp(np)*w4 
                dpedx(ix  ,iy  ,izp1)=dpedx(ix  ,iy  ,izp1)+qp(np)*w5 
                dpedx(ixp1,iy  ,izp1)=dpedx(ixp1,iy  ,izp1)+qp(np)*w6 
                dpedx(ix  ,iyp1,izp1)=dpedx(ix  ,iyp1,izp1)+qp(np)*w7 
                dpedx(ixp1,iyp1,izp1)=dpedx(ixp1,iyp1,izp1)+qp(np)*w8 
 
                p_xx (ix  ,iy  ,iz  ,is)=p_xx (ix  ,iy  ,iz  ,is)+qp(np)*w1*xx
                p_xx (ixp1,iy  ,iz  ,is)=p_xx (ixp1,iy  ,iz  ,is)+qp(np)*w2*xx 
                p_xx (ix  ,iyp1,iz  ,is)=p_xx (ix  ,iyp1,iz  ,is)+qp(np)*w3*xx 
                p_xx (ixp1,iyp1,iz  ,is)=p_xx (ixp1,iyp1,iz  ,is)+qp(np)*w4*xx 
                p_xx (ix  ,iy  ,izp1,is)=p_xx (ix  ,iy  ,izp1,is)+qp(np)*w5*xx 
                p_xx (ixp1,iy  ,izp1,is)=p_xx (ixp1,iy  ,izp1,is)+qp(np)*w6*xx 
                p_xx (ix  ,iyp1,izp1,is)=p_xx (ix  ,iyp1,izp1,is)+qp(np)*w7*xx 
                p_xx (ixp1,iyp1,izp1,is)=p_xx (ixp1,iyp1,izp1,is)+qp(np)*w8*xx 
 
                p_xy (ix  ,iy  ,iz  ,is)=p_xy (ix  ,iy  ,iz  ,is)+qp(np)*w1*xy
                p_xy (ixp1,iy  ,iz  ,is)=p_xy (ixp1,iy  ,iz  ,is)+qp(np)*w2*xy 
                p_xy (ix  ,iyp1,iz  ,is)=p_xy (ix  ,iyp1,iz  ,is)+qp(np)*w3*xy 
                p_xy (ixp1,iyp1,iz  ,is)=p_xy (ixp1,iyp1,iz  ,is)+qp(np)*w4*xy 
                p_xy (ix  ,iy  ,izp1,is)=p_xy (ix  ,iy  ,izp1,is)+qp(np)*w5*xy 
                p_xy (ixp1,iy  ,izp1,is)=p_xy (ixp1,iy  ,izp1,is)+qp(np)*w6*xy 
                p_xy (ix  ,iyp1,izp1,is)=p_xy (ix  ,iyp1,izp1,is)+qp(np)*w7*xy 
                p_xy (ixp1,iyp1,izp1,is)=p_xy (ixp1,iyp1,izp1,is)+qp(np)*w8*xy 
 
                p_xz (ix  ,iy  ,iz  ,is)=p_xz (ix  ,iy  ,iz  ,is)+qp(np)*w1*xz
                p_xz (ixp1,iy  ,iz  ,is)=p_xz (ixp1,iy  ,iz  ,is)+qp(np)*w2*xz 
                p_xz (ix  ,iyp1,iz  ,is)=p_xz (ix  ,iyp1,iz  ,is)+qp(np)*w3*xz 
                p_xz (ixp1,iyp1,iz  ,is)=p_xz (ixp1,iyp1,iz  ,is)+qp(np)*w4*xz 
                p_xz (ix  ,iy  ,izp1,is)=p_xz (ix  ,iy  ,izp1,is)+qp(np)*w5*xz 
                p_xz (ixp1,iy  ,izp1,is)=p_xz (ixp1,iy  ,izp1,is)+qp(np)*w6*xz 
                p_xz (ix  ,iyp1,izp1,is)=p_xz (ix  ,iyp1,izp1,is)+qp(np)*w7*xz 
                p_xz (ixp1,iyp1,izp1,is)=p_xz (ixp1,iyp1,izp1,is)+qp(np)*w8*xz 
 
                p_yy (ix  ,iy  ,iz  ,is)=p_yy (ix  ,iy  ,iz  ,is)+qp(np)*w1*yy
                p_yy (ixp1,iy  ,iz  ,is)=p_yy (ixp1,iy  ,iz  ,is)+qp(np)*w2*yy 
                p_yy (ix  ,iyp1,iz  ,is)=p_yy (ix  ,iyp1,iz  ,is)+qp(np)*w3*yy 
                p_yy (ixp1,iyp1,iz  ,is)=p_yy (ixp1,iyp1,iz  ,is)+qp(np)*w4*yy 
                p_yy (ix  ,iy  ,izp1,is)=p_yy (ix  ,iy  ,izp1,is)+qp(np)*w5*yy 
                p_yy (ixp1,iy  ,izp1,is)=p_yy (ixp1,iy  ,izp1,is)+qp(np)*w6*yy 
                p_yy (ix  ,iyp1,izp1,is)=p_yy (ix  ,iyp1,izp1,is)+qp(np)*w7*yy 
                p_yy (ixp1,iyp1,izp1,is)=p_yy (ixp1,iyp1,izp1,is)+qp(np)*w8*yy 
 
                p_yz (ix  ,iy  ,iz  ,is)=p_yz (ix  ,iy  ,iz  ,is)+qp(np)*w1*yz
                p_yz (ixp1,iy  ,iz  ,is)=p_yz (ixp1,iy  ,iz  ,is)+qp(np)*w2*yz 
                p_yz (ix  ,iyp1,iz  ,is)=p_yz (ix  ,iyp1,iz  ,is)+qp(np)*w3*yz 
                p_yz (ixp1,iyp1,iz  ,is)=p_yz (ixp1,iyp1,iz  ,is)+qp(np)*w4*yz 
                p_yz (ix  ,iy  ,izp1,is)=p_yz (ix  ,iy  ,izp1,is)+qp(np)*w5*yz 
                p_yz (ixp1,iy  ,izp1,is)=p_yz (ixp1,iy  ,izp1,is)+qp(np)*w6*yz 
                p_yz (ix  ,iyp1,izp1,is)=p_yz (ix  ,iyp1,izp1,is)+qp(np)*w7*yz 
                p_yz (ixp1,iyp1,izp1,is)=p_yz (ixp1,iyp1,izp1,is)+qp(np)*w8*yz 
 
                p_zz (ix  ,iy  ,iz  ,is)=p_zz (ix  ,iy  ,iz  ,is)+qp(np)*w1*zz
                p_zz (ixp1,iy  ,iz  ,is)=p_zz (ixp1,iy  ,iz  ,is)+qp(np)*w2*zz 
                p_zz (ix  ,iyp1,iz  ,is)=p_zz (ix  ,iyp1,iz  ,is)+qp(np)*w3*zz 
                p_zz (ixp1,iyp1,iz  ,is)=p_zz (ixp1,iyp1,iz  ,is)+qp(np)*w4*zz 
                p_zz (ix  ,iy  ,izp1,is)=p_zz (ix  ,iy  ,izp1,is)+qp(np)*w5*zz 
                p_zz (ixp1,iy  ,izp1,is)=p_zz (ixp1,iy  ,izp1,is)+qp(np)*w6*zz 
                p_zz (ix  ,iyp1,izp1,is)=p_zz (ix  ,iyp1,izp1,is)+qp(np)*w7*zz 
                p_zz (ixp1,iyp1,izp1,is)=p_zz (ixp1,iyp1,izp1,is)+qp(np)*w8*zz 
 
                np=link(np)
              ENDDO
            ENDDO
          ENDDO
        ENDDO


        call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(dpedx(1,jb-1,kb-1   ),NX,NY,NZ)

        call XREAL(p_xx (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_xy (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_xz (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_yy (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_yz (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(p_zz (1,jb-1,kb-1,is),NX,NY,NZ)

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

        DO IIZ=KB-1,KE+1
          DO IIY=JB-1,JE+1
            DO IIX=1,NX2
              p_xx(iix,iiy,iiz,is) = p_xx(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_xy(iix,iiy,iiz,is) = p_xy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_xz(iix,iiy,iiz,is) = p_xz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_yy(iix,iiy,iiz,is) = p_yy(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_yz(iix,iiy,iiz,is) = p_yz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
              p_zz(iix,iiy,iiz,is) = p_zz(iix,iiy,iiz,is) / (meshX%dxc(iix)*meshY%dxc(iiy+1)*meshZ%dxc(iiz+1))
            ENDDO
          ENDDO
        ENDDO
        p_xx(:,:,:,is)=p_xx(:,:,:,is)/(tx0(is)*frac(is))
        p_xy(:,:,:,is)=p_xy(:,:,:,is)/(tx0(is)*frac(is))
        p_xz(:,:,:,is)=p_xz(:,:,:,is)/(tx0(is)*frac(is))
        p_yy(:,:,:,is)=p_yy(:,:,:,is)/(tx0(is)*frac(is))
        p_yz(:,:,:,is)=p_yz(:,:,:,is)/(tx0(is)*frac(is))
        p_zz(:,:,:,is)=p_zz(:,:,:,is)/(tx0(is)*frac(is))

      ENDDO
      call date_and_time(values=time_end_array(:,26))
      call accumulate_time_difference(time_begin_array(1,26) &
     &                               ,time_end_array(1,26) &
     &                                ,time_elapsed(26))
 

!      do is=1,nspec
!        call date_and_time(values=time_begin_array(:,24))
!        call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
!        call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
!        call date_and_time(values=time_end_array(:,24))
!        call accumulate_time_difference(time_begin_array(1,24) &
!     &                                 ,time_end_array(1,24) &
!     &                                 ,time_elapsed(24))
! 
!        call date_and_time(values=time_begin_array(:,25))
!        call XREALBCC(tpar (1,jb-1,kb-1,is),1,NX,NY,NZ)
!        call XREALBCC(tperp(1,jb-1,kb-1,is),1,NX,NY,NZ)
!        call date_and_time(values=time_end_array(:,25))
!        call accumulate_time_difference(time_begin_array(1,25) &
!     &                                 ,time_end_array(1,25) &
!     &                                 ,time_elapsed(25))
! 
!      enddo


!      call date_and_time(values=time_begin_array(:,26))
!      do is=1,nspec
!        do k=kb-1,ke+1
!          do j = jb-1,je+1
!            do i=1,nx2
!              if (is == 1) then
!                dns1=dns(i,j,k,1)/(dfac(1)*frac(1))
!                dns2=0.
!                denum=dns1+rfrac*dns2
!              else
!                denum=dns(i,j,k,is)/(dfac(is)*frac(is))
!              endif
!              if (denum < denmin)  then
!               tpar(i,j,k,is)=1.e-5
!               tperp(i,j,k,is)=1.e-5
!              else
!               denum=denum*tx0(is)
!               tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
!               tperp(i,j,k,is)=0.5*tperp(i,j,k,is)*wspec(is)/denum
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
!      call date_and_time(values=time_end_array(:,26))
!      call accumulate_time_difference(time_begin_array(1,26) &
!     &                               ,time_end_array(1,26) &
!     &                               ,time_elapsed(26))

 
       call date_and_time(values=time_end_array(:,23))
       call accumulate_time_difference(time_begin_array(1,23) &
     &                                ,time_end_array(1,23) &
     &                                ,time_elapsed(23))
 
      return
      end
!
!#######################################################################
!
      subroutine xrealbcc_pack_b(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)

      use parameter_mod
      implicit none
      integer*8 ibnd,i,j,nx1m,ny1m,nz1m,k
      double precision a_x(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                      ,packed_data_xy_recv(nxmax,jb-1:je+1,3)
 
      do j=jb-1,je+1
        do i=1,nxmax
          packed_data_xy_send(i,j,1)=a_x(i,j,ke)
          packed_data_xy_send(i,j,2)=a_y(i,j,ke)
          packed_data_xy_send(i,j,3)=a_z(i,j,ke)
        enddo
      enddo
   
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrtop ,0,&
                        packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrbot ,0,&
                        mpi_comm_world,status,ierr)
      if (kb == 1) then
        if (ibnd == 1) then
          a_x(:,:,kb-1)= a_x(:,:,kb)
          a_y(:,:,kb-1)= a_y(:,:,kb)
          a_z(:,:,kb-1)= a_z(:,:,kb)
        else
          do j=jb-1,je+1
            do i=1,nxmax
              a_x(i,j,kb-1)=packed_data_xy_recv(i,j,1)
              a_y(i,j,kb-1)=packed_data_xy_recv(i,j,2)
              a_z(i,j,kb-1)=packed_data_xy_recv(i,j,3)
            enddo
          enddo
        endif
      else
        do j=jb-1,je+1
          do i=1,nxmax
            a_x(i,j,kb-1)=packed_data_xy_recv(i,j,1)
            a_y(i,j,kb-1)=packed_data_xy_recv(i,j,2)
            a_z(i,j,kb-1)=packed_data_xy_recv(i,j,3)
          enddo
        enddo
      endif
 
 
      do j=jb-1,je+1
        do i=1,nxmax
          packed_data_xy_send(i,j,1)=a_x(i,j,kb)
          packed_data_xy_send(i,j,2)=a_y(i,j,kb)
          packed_data_xy_send(i,j,3)=a_z(i,j,kb)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrbot ,1,&
                        packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrtop ,1,&
                        mpi_comm_world,status,ierr)
      if (ke == nz1m) then
        if (ibnd == 1) then
          a_x(:,:,ke+1)= a_x(:,:,ke)
          a_y(:,:,ke+1)= a_y(:,:,ke)
          a_z(:,:,ke+1)= a_z(:,:,ke)
        else
          do j=jb-1,je+1
            do i=1,nxmax
              a_x(i,j,ke+1)=packed_data_xy_recv(i,j,1)
              a_y(i,j,ke+1)=packed_data_xy_recv(i,j,2)
              a_z(i,j,ke+1)=packed_data_xy_recv(i,j,3)
            enddo
          enddo
        endif
      else
        do j=jb-1,je+1
          do i=1,nxmax
            a_x(i,j,ke+1)=packed_data_xy_recv(i,j,1)
            a_y(i,j,ke+1)=packed_data_xy_recv(i,j,2)
            a_z(i,j,ke+1)=packed_data_xy_recv(i,j,3)
          enddo
        enddo
      endif
 
 
      do k=kb-1,ke+1
        do i=1,nxmax
          packed_data_xz_send(i,k,1)=a_x(i,je,k)
          packed_data_xz_send(i,k,2)=a_y(i,je,k)
          packed_data_xz_send(i,k,3)=a_z(i,je,k)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrrite,0,&
                        packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrleft,0,&
                        mpi_comm_world,status,ierr)
      if (jb == 1) then
        if (ibnd == 1)  then
          a_x(:,jb-1,:)= a_x(:,jb,:)
          a_y(:,jb-1,:)= a_y(:,jb,:)
          a_z(:,jb-1,:)= a_z(:,jb,:)
        else
          do k=kb-1,ke+1
            do i=1,nxmax
              a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
              a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
              a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
            enddo
          enddo
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
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrleft,0,&
                        packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrrite,0,&
                        mpi_comm_world,status,ierr)
      if (je == ny1m) then
        if (ibnd == 1)  then
          a_x(:,je+1,:)= a_x(:,je,:)
          a_y(:,je+1,:)= a_y(:,je,:)
          a_z(:,je+1,:)= a_z(:,je,:)
        else
          do k=kb-1,ke+1
            do i=1,nxmax
              a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
              a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
              a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
            enddo
          enddo
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
      subroutine get_cleanup_status(maxchar)
 
      use parameter_mod
 
      if (myid==0) then
!        open(unit=1,file='.cleanup_status',status='old')
        open(unit=1,file=trim(adjustl(data_directory))//'.cleanup_status',status='old')
        read(1,*) cleanup_status
        close(unit=1)
      endif
      call MPI_BCAST(cleanup_status,maxchar,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
!
!=======================================================================
!
      end subroutine get_cleanup_status
!
!#######################################################################
!
      subroutine xrealbcc_pack_e(a_x,a_y,a_z, ibnd, nx1m, ny1m,nz1m)

      use parameter_mod
      implicit none
      integer*8 ibnd,i,j,nx1m,ny1m,nz1m,k
      double precision a_x(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,a_y(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,a_z(nxmax,jb-1:je+1,kb-1:ke+1)&
                      ,packed_data_xz_send(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xz_recv(nxmax,kb-1:ke+1,3) &
                      ,packed_data_xy_send(nxmax,jb-1:je+1,3) &
                      ,packed_data_xy_recv(nxmax,jb-1:je+1,3)
 
      do j=jb-1,je+1
        do i=1,nxmax
          packed_data_xy_send(i,j,1)=a_x(i,j,ke)
          packed_data_xy_send(i,j,2)=a_y(i,j,ke)
          packed_data_xy_send(i,j,3)=a_z(i,j,ke)
        enddo
      enddo
   
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrtop ,0,&
                        packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrbot ,0,&
                        mpi_comm_world,status,ierr)
      if (kb == 1) then
        if (ibnd == 1) then
          a_x(:,:,kb-1)= a_x(:,:,kb)
          a_y(:,:,kb-1)= a_y(:,:,kb)
          a_z(:,:,kb-1)=+a_z(:,:,kb)
          a_x(:,:,kb-1)= exc
          a_y(:,:,kb-1)= eyc
          a_z(:,:,kb-1)=+ezc
        else
          do j=jb-1,je+1
            do i=1,nxmax
              a_x(i,j,kb-1)=packed_data_xy_recv(i,j,1)
              a_y(i,j,kb-1)=packed_data_xy_recv(i,j,2)
              a_z(i,j,kb-1)=packed_data_xy_recv(i,j,3)
            enddo
          enddo
        endif
      else
        do j=jb-1,je+1
          do i=1,nxmax
            a_x(i,j,kb-1)=packed_data_xy_recv(i,j,1)
            a_y(i,j,kb-1)=packed_data_xy_recv(i,j,2)
            a_z(i,j,kb-1)=packed_data_xy_recv(i,j,3)
          enddo
        enddo
      endif
 
 
      do j=jb-1,je+1
        do i=1,nxmax
          packed_data_xy_send(i,j,1)=a_x(i,j,kb)
          packed_data_xy_send(i,j,2)=a_y(i,j,kb)
          packed_data_xy_send(i,j,3)=a_z(i,j,kb)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),MPI_DOUBLE_PRECISION,nbrbot ,1,&
                        packed_data_xy_recv,size(packed_data_xy_recv),MPI_DOUBLE_PRECISION,nbrtop ,1,&
                        mpi_comm_world,status,ierr)
      if (ke == nz1m) then
        if (ibnd == 1) then
          a_x(:,:,ke+1)= a_x(:,:,ke)
          a_y(:,:,ke+1)= a_y(:,:,ke)
          a_z(:,:,ke+1)=+a_z(:,:,ke)
          a_x(:,:,ke+1)= exc
          a_y(:,:,ke+1)= eyc
          a_z(:,:,ke+1)=+ezc
        else
          do j=jb-1,je+1
            do i=1,nxmax
              a_x(i,j,ke+1)=packed_data_xy_recv(i,j,1)
              a_y(i,j,ke+1)=packed_data_xy_recv(i,j,2)
              a_z(i,j,ke+1)=packed_data_xy_recv(i,j,3)
            enddo
          enddo
        endif
      else
        do j=jb-1,je+1
          do i=1,nxmax
            a_x(i,j,ke+1)=packed_data_xy_recv(i,j,1)
            a_y(i,j,ke+1)=packed_data_xy_recv(i,j,2)
            a_z(i,j,ke+1)=packed_data_xy_recv(i,j,3)
          enddo
        enddo
      endif
 
 
      do k=kb-1,ke+1
        do i=1,nxmax
          packed_data_xz_send(i,k,1)=a_x(i,je,k)
          packed_data_xz_send(i,k,2)=a_y(i,je,k)
          packed_data_xz_send(i,k,3)=a_z(i,je,k)
        enddo
      enddo
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrrite,0,&
                        packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrleft,0,&
                        mpi_comm_world,status,ierr)
      if (jb == 1) then
        if (ibnd == 1)  then
          a_x(:,jb-1,:)= a_x(:,jb,:)
!          a_y(:,jb-1,:)=-a_y(:,jb,:)
          a_y(:,jb-1,:)= a_y(:,jb,:)

          a_z(:,jb-1,:)= a_z(:,jb,:)

          a_x(:,jb-1,:)= exc
          a_y(:,jb-1,:)= eyc
          a_z(:,jb-1,:)= ezc
        else
          do k=kb-1,ke+1
            do i=1,nxmax
              a_x(i,jb-1,k)=packed_data_xz_recv(i,k,1)
              a_y(i,jb-1,k)=packed_data_xz_recv(i,k,2)
              a_z(i,jb-1,k)=packed_data_xz_recv(i,k,3)
            enddo
          enddo
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
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),MPI_DOUBLE_PRECISION,nbrleft,0,&
                        packed_data_xz_recv,size(packed_data_xz_recv),MPI_DOUBLE_PRECISION,nbrrite,0,&
                        mpi_comm_world,status,ierr)
      if (je == ny1m) then
        if (ibnd == 1)  then
          a_x(:,je+1,:)= a_x(:,je,:)
!          a_y(:,je+1,:)=-a_y(:,je,:)
          a_y(:,je+1,:)= a_y(:,je,:)
          a_z(:,je+1,:)= a_z(:,je,:)

          a_x(:,je+1,:)= exc
          a_y(:,je+1,:)= eyc
          a_z(:,je+1,:)= ezc
        else
          do k=kb-1,ke+1
            do i=1,nxmax
              a_x(i,je+1,k)=packed_data_xz_recv(i,k,1)
              a_y(i,je+1,k)=packed_data_xz_recv(i,k,2)
              a_z(i,je+1,k)=packed_data_xz_recv(i,k,3)
            enddo
          enddo
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
      subroutine MESH_INTERPOLATED_3D(nonuniform_mesh,uniform_mesh,nonuniform_mesh_global)
 
      use parameter_mod
      use MESH2D

      double precision:: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,w1,w2,w3,w4,w5,w6,w7,w8
      integer*8:: ix,iy,iz,ixp1,iyp1,izp1,i,j,k,jmin,jmax,kmin,kmax
      double precision,dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(in) :: nonuniform_mesh
      double precision,dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(out):: uniform_mesh
      double precision,dimension(nxmax,0:ny+1,0:nz+1), intent(out):: nonuniform_mesh_global
      double precision,dimension(nxmax,0:ny+1,0:nz+1):: nonuniform_mesh_local
      double precision:: xc_uniform_pos,yc_uniform_pos,zc_uniform_pos

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt

      uniform_mesh          = 0.
      nonuniform_mesh_local = 0.

      if (jb == 1) then
        jmin = 0
      else 
        jmin = jb
      endif
      if (je == ny) then
        jmax = ny+1
      else 
        jmax = je
      endif
      if (kb == 1) then
        kmin = 0
      else 
        kmin = kb
      endif
      if (ke == nz) then
        kmax = nz+1
      else 
        kmax = ke
      endif

      do i=1,nxmax
        do j=jmin,jmax
          do k=kmin,kmax
            nonuniform_mesh_local(i,j,k)=nonuniform_mesh(i,j,k)
          enddo
        enddo
      enddo
      call MPI_ALLREDUCE( nonuniform_mesh_local,nonuniform_mesh_global,size(nonuniform_mesh_local)        &
                         ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

      do i=2,nx+1
        xc_uniform_pos = (i-1.5)*hx
        rx   = dtxi*MESH_UNMAP(meshX,xc_uniform_pos)+1.50000000000d+00
        ix   = rx
        fx   = rx-ix
        ixp1 = ix+1
        do j=jb,je
          yc_uniform_pos = (j-0.5)*hy
          ry   = dtyi*MESH_UNMAP(meshY,yc_uniform_pos)+1.50000000000d+00
          iy   = ry
          fy   = ry-iy
          iy   = iy-1             ! integer index in y direction starts at 0
          iyp1 = iy+1
          do k=kb,ke
            zc_uniform_pos = (k-0.5)*hz
            rz   = dtzi*MESH_UNMAP(meshZ,zc_uniform_pos)+1.50000000000d+00
            iz   = rz
            fz   = rz-iz
            iz   = iz-1             ! integer index in z direction starts at 0
            izp1 = iz+1
       
            w1=(1.-fx)*(1.-fy)*(1.-fz)
            w2=fx     *(1.-fy)*(1.-fz)
            w3=(1.-fx)*fy     *(1.-fz)
            w4=fx     *fy     *(1.-fz)
            w5=(1.-fx)*(1.-fy)*fz
            w6=fx     *(1.-fy)*fz
            w7=(1.-fx)*fy     *fz
            w8=fx     *fy     *fz
 
            uniform_mesh(i,j,k) =  w1 * nonuniform_mesh_global(ix  ,iy  ,iz  )     &
                                 + w2 * nonuniform_mesh_global(ixp1,iy  ,iz  )     &
                                 + w3 * nonuniform_mesh_global(ix  ,iyp1,iz  )     &
                                 + w4 * nonuniform_mesh_global(ixp1,iyp1,iz  )     &
                                 + w5 * nonuniform_mesh_global(ix  ,iy  ,izp1)     &
                                 + w6 * nonuniform_mesh_global(ixp1,iy  ,izp1)     &
                                 + w7 * nonuniform_mesh_global(ix  ,iyp1,izp1)     &
                                 + w8 * nonuniform_mesh_global(ixp1,iyp1,izp1)
          enddo
        enddo
      enddo

      return
      end
!
!#######################################################################
!
      subroutine particle_in_volume_write
 
      use parameter_mod
      use MESH2D
      integer*8 count_kbq,time_begin(8),time_end(8)
      integer*8 nptotp_kbq,npart_kbq(2),np_ijk,Storage_Error_p,Storage_Error
      data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
      data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
      data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
      integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                 ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
      double precision:: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8,xpart,ypart,zpart,r_particle
      double precision:: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi
      double precision, dimension(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bx_av,by_av,bz_av
      double precision:: v_limit,eps2,rx0,ry0,rz0,rrat,sqrr,outer_radius,myranf,twopi,fluxran,vxa,vyz,vza
      INTEGER*8:: L, EXIT_CODE_P, EXIT_CODE
      integer*8:: n_fast_removed,n_fast_removed_local,nptot_max,Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p
      double precision:: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min,x_disp,y_disp,z_disp          &
                        ,y_disp_max_p,x_disp_max_p,z_disp_max_p,y_disp_max,x_disp_max,z_disp_max
      double precision:: disp_max_p(3),disp_max(3),tx,ty,tz,v_x,v_y,v_z  
      INTEGER*4 :: nescapearr(8),nescapearr_global(8)
      INTEGER*4 :: ppacket(3),ppacketg(3),dpacket(4),dpacketg(4)
      INTEGER*8 :: epacket(2),epacketg(2),indx,loop
      INTEGER*8,dimension(:),allocatable :: nparr
      double precision, dimension(3,nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bxyz_av
      double precision:: TEX1,TEX2,TEX3,TEX4,TEX5,TEX6,TEX7,TEX8  
      double precision:: TEY1,TEY2,TEY3,TEY4,TEY5,TEY6,TEY7,TEY8  
      double precision:: TEZ1,TEZ2,TEZ3,TEZ4,TEZ5,TEZ6,TEZ7,TEZ8  
      double precision:: mX_xa,mX_ta,mX_ca1,mX_ca2,mX_xb,mX_dtdx,mX_tb,mX_cb1,mX_cb2
      double precision:: mY_xa,mY_ta,mY_ca1,mY_ca2,mY_xb,mY_dtdx,mY_tb,mY_cb1,mY_cb2
      double precision:: mZ_xa,mZ_ta,mZ_ca1,mZ_ca2,mZ_xb,mZ_dtdx,mZ_tb,mZ_cb1,mZ_cb2
      character(len=160):: filename
      integer:: N_IN_VOLUME
      integer iErr1,iErr2,file,eStrLen,subArray,stat(MPI_STATUS_SIZE),mode
      character eStr*(1024)
      integer:: N_TOTAL,N_MAX,RECNUM,LENREC,FILENUM,IP
      real*4, dimension(:), allocatable:: xw,yw,zw,vxw,vyw,vzw,qw,vwpar,vwperp1,vwperp2

      dtxi = 1./meshX%dt
      dtyi = 1./meshY%dt
      dtzi = 1./meshZ%dt
 
      N_IN_VOLUME = 0
      DO IS=1,NSPEC
        NPTOTP=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                NP=IPHEAD(ixe,iye,ize,is)
                DO WHILE (NP.NE.0)
                  NPTOTP=NPTOTP+1
                  IF ( X(NP) <= XBOX_R .AND. X(NP) >= XBOX_L .AND.         &
                       Y(NP) <= YBOX_R .AND. Y(NP) >= YBOX_L .AND.         &
                       Z(NP) <= ZBOX_R .AND. Z(NP) >= ZBOX_L                 ) N_IN_VOLUME = N_IN_VOLUME + 1
                  NP=LINK(NP)
                ENDDO
            enddo
          enddo
        enddo
      ENDDO

      call MPI_ALLREDUCE(N_IN_VOLUME,N_TOTAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(N_IN_VOLUME,N_MAX  ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
      allocate ( XW(N_MAX), YW(N_MAX), ZW(N_MAX), VXW(N_MAX), VYW(N_MAX), VZW(N_MAX), QW(N_MAX)       &
                ,VWPAR(N_MAX),VWPERP1(N_MAX),VWPERP2(N_MAX))

      N_IN_VOLUME = 0
      DO IS=1,NSPEC
        NPTOTP=0
        do ize=kb-1,ke
          do iye=jb-1,je
            do ixe=1,nx1   
                NP=IPHEAD(ixe,iye,ize,is)
                L=NP
                DO WHILE (NP.NE.0)
                  NPTOTP=NPTOTP+1
                  IF ( X(NP) <= XBOX_R .AND. X(NP) >= XBOX_L .AND.         &
                       Y(NP) <= YBOX_R .AND. Y(NP) >= YBOX_L .AND.         &
                       Z(NP) <= ZBOX_R .AND. Z(NP) >= ZBOX_L                 ) THEN
                      N_IN_VOLUME = N_IN_VOLUME + 1
                      XW(N_IN_VOLUME)  = X(NP)
                      YW(N_IN_VOLUME)  = Y(NP)
                      ZW(N_IN_VOLUME)  = Z(NP)
                      VXW(N_IN_VOLUME) = VX(NP)
                      VYW(N_IN_VOLUME) = VY(NP)
                      VZW(N_IN_VOLUME) = VZ(NP)
                      QW(N_IN_VOLUME)  = QP(NP)


!                   Uniform mesh - Same as is in version 5.0
!                    rx=hxi*x(l)+1.5000000000000001
!                    ry=hyi*y(l)+0.5000000000000001d+00
!                    rz=hzi*z(l)+0.5000000000000001d+00
!                    ix=rx
!                    iy=ry
!                    iz=rz
!                    fx=rx-ix
!                    fy=ry-iy
!                    fz=rz-iz

!                   Nonuniform mesh - using MESH_UNMAP
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

                    w1=(1.-fx)*(1.-fy)*(1.-fz)
                    w2=    fx *(1.-fy)*(1.-fz)
                    w3=(1.-fx)*    fy *(1.-fz)
                    w4=    fx*     fy *(1.-fz)
                    w5=(1.-fx)*(1.-fy)*    fz
                    w6=    fx *(1.-fy)*    fz
                    w7=(1.-fx)*    fy*     fz
                    w8=    fx*     fy*     fz

                    vxa=vx(l)
                    vya=vy(l)
                    vza=vz(l)

                    bxa  =bx(ix  ,iy  ,iz  )*w1+bx(ixp1,iy  ,iz  )*w2  &
                    +     bx(ix  ,iyp1,iz  )*w3+bx(ixp1,iyp1,iz  )*w4  &
                    +     bx(ix  ,iy  ,izp1)*w5+bx(ixp1,iy  ,izp1)*w6  &
                    +     bx(ix  ,iyp1,izp1)*w7+bx(ixp1,iyp1,izp1)*w8

                    bya  =by(ix  ,iy  ,iz  )*w1+by(ixp1,iy  ,iz  )*w2  &
                    +     by(ix  ,iyp1,iz  )*w3+by(ixp1,iyp1,iz  )*w4  &
                    +     by(ix  ,iy  ,izp1)*w5+by(ixp1,iy  ,izp1)*w6  &
                    +     by(ix  ,iyp1,izp1)*w7+by(ixp1,iyp1,izp1)*w8

                    bza  =bz(ix  ,iy  ,iz  )*w1+bz(ixp1,iy  ,iz  )*w2  &
                    +     bz(ix  ,iyp1,iz  )*w3+bz(ixp1,iyp1,iz  )*w4  &
                    +     bz(ix  ,iy  ,izp1)*w5+bz(ixp1,iy  ,izp1)*w6  &
                    +     bz(ix  ,iyp1,izp1)*w7+bz(ixp1,iyp1,izp1)*w8

                    btota=sqrt(bxa**2+bya**2+bza**2)
                    if (btota < 1.e-20) btota=1.e-20

                    par_x   = bxa / btota
                    par_y   = bya / btota
                    par_z   = bza / btota
                    arb_x   = 0.
                    arb_y   = 0.
                    arb_z   = 1.
                    perp1_x = arb_y*par_z   - arb_z*par_y
                    perp1_y = arb_z*par_x   - arb_x*par_z
                    perp1_z = arb_x*par_y   - arb_y*par_x
                    perp2_x = perp1_y*par_z - perp1_z*par_y
                    perp2_y = perp1_z*par_x - perp1_x*par_z
                    perp2_z = perp1_x*par_y - perp1_y*par_x

                    vwpar  (N_IN_VOLUME)  = vxa*par_x   + vya*par_y   + vza*par_z
                    vwperp1(N_IN_VOLUME)  = vxa*perp1_x + vya*perp1_y + vza*perp1_z
                    vwperp2(N_IN_VOLUME)  = vxa*perp2_x + vya*perp2_y + vza*perp2_z
                     

                  ENDIF
                  NP=LINK(NP)
                ENDDO
            enddo
          enddo
        enddo
      ENDDO
      
      lenrec=10*recl_for_real
      if (myid == 0) then
          filenum = 101
          open (filenum,                                                       &
                file= trim(adjustl(data_directory))//'particle_'//             &
                trim(adjustl(cycle_ascii))//'.bin',                            &
                form='unformatted',                                            &
                action='write',access='direct', status='unknown',recl=lenrec)
         write(6,*) " FILE_NAME = ",trim(adjustl(data_directory))//'particle_'//             &
                trim(adjustl(cycle_ascii))//'.bin'
          recnum = 1
          write(filenum,rec=recnum) N_TOTAL
          recnum = recnum + 1
          IF (N_IN_VOLUME /= 0) THEN
            write(filenum,rec=recnum) N_IN_VOLUME
            recnum = recnum + 1
            DO IP=1,N_IN_VOLUME
              write(filenum,rec=recnum) XW(IP),YW(IP),ZW(IP),VXW(IP),VYW(IP),VZW(IP),QW(IP)                &
                                       ,VWPAR(IP),VWPERP1(IP),VWPERP2(IP)
              recnum = recnum + 1
            ENDDO
          ENDIF
      endif

      if (myid == 0) then
        do ipe=1,numprocs-1
           CALL MPI_RECV(N_IN_VOLUME,1,MPI_INTEGER,ipe,ipe,MPI_COMM_WORLD,status,ierror)
           IF (N_IN_VOLUME /= 0) THEN
             CALL MPI_RECV(XW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(YW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(ZW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(VXW    ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(VYW    ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(VZW    ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(QW     ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(VWPAR  ,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(VWPERP1,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             CALL MPI_RECV(VWPERP2,N_IN_VOLUME,MPI_REAL4,ipe,ipe+numprocs,MPI_COMM_WORLD,status,ierror)
             write(filenum,rec=recnum) N_IN_VOLUME
             recnum = recnum + 1
             DO IP=1,N_IN_VOLUME
               write(filenum,rec=recnum) XW(IP),YW(IP),ZW(IP),VXW(IP),VYW(IP),VZW(IP),QW(IP)                &
                                        ,VWPAR(IP),VWPERP1(IP),VWPERP2(IP)
               recnum = recnum + 1
             ENDDO
           ENDIF
        enddo
        CLOSE(UNIT=FILENUM)
      else
         CALL MPI_SEND(N_IN_VOLUME,1,MPI_INTEGER,0,MYID,MPI_COMM_WORLD,ierror)
         IF (N_IN_VOLUME /= 0) THEN
           CALL MPI_SEND(XW     ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(YW     ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(ZW     ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(VXW    ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(VYW    ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(VZW    ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(QW     ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(VWPAR  ,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(VWPERP1,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
           CALL MPI_SEND(VWPERP2,N_IN_VOLUME,MPI_REAL4,0,myid+numprocs,MPI_COMM_WORLD,ierror)
         ENDIF
      endif

      deallocate (XW,YW,ZW,VXW,VYW,VZW,QW,VWPAR,VWPERP1,VWPERP2)

      return
      end
