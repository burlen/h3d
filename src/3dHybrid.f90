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
      integer*8 :: time_begin(8),time_end(8),time_per_cycle,input_error,is
      integer*8 :: clock_time_re1, itstart, itfinish, nwrtrestart, irecnum
      integer*8 :: nzl, nyl, nplmax6, lenrec, ixe, iye, ize
      integer*8 :: i, j, k, jbt, jet, kbt, ket
      integer*4 :: numvars
      integer*8, pointer :: itest
      real(kind=CUSTOM_REAL) :: rnorm
      real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable:: uniform_mesh,nonuniform_mesh_global	!LAURA
      real*8 :: time1,time2
      character (len=240)::fileNameX
      integer:: iwrite
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
      
      external get_environment_variable
!
!***********************************************************************
!
       namelist /datum/ denmin,resis,iterb,testorbt,norbskip                &
       ,restart,nspec,nx,xmax,ny,ymax,nz,zmax,npx,npy,npz                   &
       ,dt,nprint,nwrtdata,qspec,wspec,nskipx                               &
       ,nskipy,nskipz,bxc,byc,bzc,frac,vxbar,vybar,vzbar,anisot             &
       ,btspec,bete,fxsho,nxcel,phib,netax,netay,netaz,ave1,ave2            &
       ,wpiwci,theta,restrt_write,Yee,b_in_xy,global,dtwci,dipole_moment    &
       ,dipole_sphere_radius,i_dipole,j_dipole,k_dipole,absorbing_dipole    &
       ,harris,rcorr,ishape,teti,nwrtrestart,ieta,etamin,etamax             &
       ,R_MP,R_obstacle_to_MP,eta_par,bz_IMF,gama,QUOTA                     &
       ,maximum_simulation_time,n_subcycles,buffer_zone                     &
       ,xaa,xbb,nax,nbx,yaa,ybb,nay,nby,zaa,zbb,naz,nbz,setup_mesh          &
       ,profile_power,uniform_loading_in_logical_grid,moat_zone,            &
       MPI_IO_format,insitu

        time_elapsed=0._CUSTOM_REAL;time_begin_array=0;time_end_array=0
        buffer_zone=0._CUSTOM_REAL
        notime=1
        bxc=0._CUSTOM_REAL;byc=0._CUSTOM_REAL;bzc=0._CUSTOM_REAL
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
      iwt=0;nskipx=1;nskipy=1;nskipz=1;testorbt=.false.;pi=dacos(-1.E+00_CUSTOM_REAL);frac=1.E+00_CUSTOM_REAL;t_stopped=0._CUSTOM_REAL !LAURA
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

      nwrtrestart=5000

      call MPI_BCAST(denmin                 ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(resis                  ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(iterb                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(testorbt               ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(MPI_IO_format          ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(norbskip               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(restart                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nspec                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nx                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ny                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nz                     ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(xmax                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(ymax                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(zmax                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(npx                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(npy                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(npz                    ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dt                     ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nprint                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nwrtdata               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(qspec                  ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(wspec                  ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nskipx                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nskipy                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nskipz                 ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(bxc                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(byc                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(bzc                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(frac                   ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(vxbar                  ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(vybar                  ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(vzbar                  ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(anisot                 ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(gama                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(btspec                 ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(bete                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(fxsho                  ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nxcel                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(phib                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(netax                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(netay                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(netaz                  ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ave1                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(ave2                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(wpiwci                 ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(theta                  ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(restrt_write           ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(b_in_xy                ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(global                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dtwci                  ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(dipole_moment          ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(dipole_sphere_radius   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(i_dipole               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(j_dipole               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(k_dipole               ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(absorbing_dipole       ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(rcorr                  ,5     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(ishape                 ,5     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(teti                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nwrtrestart            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(ieta                   ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(etamin                 ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(etamax                 ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(R_MP                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(R_obstacle_to_MP       ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(eta_par                ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(bz_IMF                 ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(Yee                    ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(b_in_xy                ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(global                 ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(absorbing_dipole       ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(harris                 ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dtwci  ,1              ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(dipole_moment          ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(dipole_sphere_radius   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(R_MP                   ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(R_obstacle_to_MP       ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nwrtrestart            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(QUOTA                  ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(maximum_simulation_time,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(n_subcycles            ,1     ,MPI_INTEGER8         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(buffer_zone            ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(xaa                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(xbb                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nax                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nbx                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(yaa                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(ybb                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(nay                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nby                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(zaa                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(zbb                    ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
      call MPI_BCAST(naz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nbz                    ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(setup_mesh             ,1     ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(profile_power          ,1     ,MPI_INTEGER8        ,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(uniform_loading_in_logical_grid,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(moat_zone              ,1     ,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,IERR)	!LAURA
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
!    coprocessorinitialize(python insitu script filename,
!                          length of python insitu script filename)
!
!----------------------------------------------------------------------
      if (insitu) then
         call coprocessorinitialize("insitu.py", 9)
      endif
!----------------------------------------------------------------------
! end Initialize
!----------------------------------------------------------------------
!
!
       vxbar=vxbar/wpiwci
       vybar=vybar/wpiwci
       vzbar=vzbar/wpiwci
       dt=dtwci*wpiwci
       if (nz == 1) then
         dipole_moment=dsqrt(2._CUSTOM_REAL)*vxbar(1)*(R_MP**2_CUSTOM_REAL)
       else
         dipole_moment=dsqrt(2._CUSTOM_REAL)*vxbar(1)*(R_MP**3_CUSTOM_REAL)
       endif
       dipole_sphere_radius=R_MP*R_obstacle_to_MP

       if (myid==0) then
         write(6,*) " dipole_moment        = ",dipole_moment
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
      call MPI_TYPE_VECTOR(nzl+2,nx+2,(nx+2)*(nyl+2),CUSTOM_MPI_TYPE,stridery,IERR)	!LAURA
      call MPI_TYPE_COMMIT(stridery,IERR)
      call MPI_TYPE_VECTOR(nyl+2,nx+2,nx+2          ,CUSTOM_MPI_TYPE,STRIDERZ,IERR)	!LAURA
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


!      if (myid == 0) then
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

!        open(unit=11,file='mesh_vertices.dat',status='unknown',form='formatted')
!
!        write(11,*) meshX%nl+1,meshY%nl+1,meshZ%nl+1
!        do i=2,meshX%nl+2 
!          write(11,*) meshX%xn(i)
!        enddo
!        do i=2,meshX%nl+2 
!          write(11,*) meshX%dxc(i)
!        enddo
!
!        do i=2,meshY%nl+2 
!          write(11,*) meshY%xn(i)
!        enddo
!        do i=2,meshY%nl+2 
!          write(11,*) meshY%dxc(i)
!        enddo
!
!        do i=2,meshZ%nl+2 
!          write(11,*) meshZ%xn(i)
!        enddo
!        do i=2,meshZ%nl+2 
!          write(11,*) meshZ%dxc(i)
!        enddo
!         
!        close(unit=11)
!      endif

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
      clock_time_re1=(time_begin(5)*3600._CUSTOM_REAL+time_begin(6)*60._CUSTOM_REAL+time_begin(7)+time_begin(8)*0.001_CUSTOM_REAL)
 
!      if (myid == 0) inquire(file='restart_index.dat',exist=restart)
      if (myid == 0) inquire(file=trim(adjustl(restart_directory))//'restart_index.dat',exist=restart)
      call MPI_BCAST(restart         ,1     ,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
   
      if (restart) then
 
        if (myid == 0) then
!          open(unit=222,file='restart_index.dat' ,status='old')
          open(unit=222,file=trim(adjustl(restart_directory))//'restart_index.dat' ,status='old')
          read(222,*) restart_index,itfin
          close(222)
        endif

        call MPI_BCAST(restart_index,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
        call MPI_BCAST(itfin        ,1,MPI_INTEGER8,0,MPI_COMM_WORLD,IERR)
        nplmax6 = 6*nplmax

 
!       comment out for timing on LANL machine

        do iwrite = 0,npes_over_60 
         if (mod(myid,npes_over_60 + 1).eq.iwrite) then
           call restrtrw(-1.0_4,itstart)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,IERR)
        enddo

        if (restart_index == 1) then
          restart_index=2
        else
          restart_index=1
        endif

!       Uniform mesh - Same as is in version 5.0
        yb=dble(jb-1)*hy
        ye=dble(je)   *hy
        zb=dble(kb-1)*hz
        ze=dble(ke)   *hz

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
          dfac(is)=dble(ny*nz*nx)/dble(npm)
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
          xc_uniform(i) = hx*(dble(i)-1.5_CUSTOM_REAL)
          xv_uniform(i) = hx*(dble(i)-2.0_CUSTOM_REAL)
        enddo
        do j=1,nymax
          yc_uniform(j) = hy*(dble(j)-0.5_CUSTOM_REAL)
          yv_uniform(j) = hy*(dble(j)-1.0_CUSTOM_REAL)
        enddo
        do k=1,nzmax
          zc_uniform(k) = hz*(dble(k)-0.5_CUSTOM_REAL)
          zv_uniform(k) = hz*(dble(k)-1.0_CUSTOM_REAL)
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
          open (101,                                                                                  &
                file= trim(adjustl(data_directory))//'bx_dipole.gda', &
                form='unformatted',                                                                   &
                action='write',access='direct', status='unknown',recl=lenrec)
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,101,irecnum,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(bdipole_y,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bdipole_y
        if (MPI_IO_format) then
	  fileNameX = trim(adjustl(data_directory))//'by_dipole.gda'
          call wrtfile         (uniform_mesh,rnorm,fileNameX,irecnum,ny,nz)
        else
          open (101,                                                                                  &
                file= trim(adjustl(data_directory))//'bx_dipole.gda', &
                form='unformatted',                                                                   &
                action='write',access='direct', status='unknown',recl=lenrec)
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,101,irecnum,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(bdipole_z,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bdipole_z
        if (MPI_IO_format) then
	  fileNameX = trim(adjustl(data_directory))//'bz_dipole.gda'
          call wrtfile         (uniform_mesh,rnorm,fileNameX,irecnum,ny,nz)
        else
          open (101,                                                                                  &
                file= trim(adjustl(data_directory))//'bx_dipole.gda', &
                form='unformatted',                                                                   &
                action='write',access='direct', status='unknown',recl=lenrec)
          call wrtfile_NON_MPIO(uniform_mesh,rnorm,101,irecnum,ny,nz)
        endif

      call date_and_time(values=time_end)
      clock_time_init=( time_end(5)*3600._CUSTOM_REAL+time_end(6)*60._CUSTOM_REAL+time_end(7)+time_end(8)*0.001_CUSTOM_REAL)
      if (myid == 0) then
        print *,'load time = ',dble(clock_time_init-clock_time_re1)  
      endif
      clock_time_old = clock_time_init

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
      time_elapsed=0._CUSTOM_REAL;time_begin_array=0;time_end_array=0

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
           call restrtrw(1.0_4,itstart)
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
          write(file_unit_time,"(i4,' begin    ',e20.8)") it,dble(clock_time-clock_time_init)
        endif
        if (myid == 0.and.mod(it,10) == 0) then
          print *,"it = ",it
          print *,'system time (delta) = ',dble(clock_time - clock_time_old)
          print *,'system time (total) = ',dble(clock_time - clock_time_init)
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
          write(file_unit_time,"(i4,' trans    ',e20.8)") it,dble(clock_time-clock_time_init)
        endif
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
        if (notime == 0) then
          write(file_unit_time,"(i4,' injctpar ',e20.8)") it,dble(clock_time-clock_time_init)
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
          write(file_unit_time,"(i4,' sortit   ',e20.8)") it,dble(clock_time-clock_time_init)
        endif

 
        call date_and_time(values=time_begin_array(:,4))
        if (mod(it,10) == 0) call sortit    !  sort the particles
        call date_and_time(values=time_end_array(:,4))
 
 
        call date_and_time(values=time_end)
        clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 

        if (notime == 0) then
          write(file_unit_time,"(i4,' field    ',e20.8)") it,dble(clock_time-clock_time_init)
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

         if (myid == 0 .and. .not. MPI_IO_format) call openfiles

          call date_and_time(values=time_begin_array(:,6))
          if (ndim /= 1) then
            call caltemp2_global
          else
            call caltemp2_global_2d
          endif
          call date_and_time(values=time_end_array(:,6))
          numvars = 12
          irecnum=1
          call dataout(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,tperp,   &
                      nxmax,nymax,nzmax,file_unit,myid,          &
	              numvars,irecnum,kb,ke,numprocs,wpiwci,jb,je,ny,nz,nylmax,nzlmax,nspecm, &
                      bdipole_x,bdipole_y,bdipole_z,eta,eta_times_b_dot_j,eta_par,            &
                      uniform_mesh,nonuniform_mesh_global,trim(adjustl(data_directory)), trim(adjustl(cycle_ascii)),      &
                      MPI_IO_format)
         if (myid == 0 .and. .not. MPI_IO_format) then
           do j=1,13
             close(file_unit(j))
           enddo
         endif

          if (cleanup_status == 'CLEANUP_STATUS=TRUE') goto 999
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
         call addscalars(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,      &
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
!       maximum_simulation_time
        nwrtrestart = 5000
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
           call restrtrw(1.0_4,itstart)
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
!=======================================================================
!
      call MPI_FINALIZE(IERR)
      stop
      end
!
!***********************************************************************
!
      subroutine addscalars(bx,by,bz,den,ex,ey,ez,vix,viy,viz,tpar,   &
                            tperp,eta,nxmax,nylmax,nzlmax,            &
                            kb,ke,jb,je,nspecm)
!======================================================================
!
! insitu Patrick O'Leary 2/21/2013
!
!  addscalars - adding raw scalars to insitu pipeline.
!
!======================================================================
      use PRECISION
      implicit none
      integer*8 nxmax,nylmax,nzlmax,kb,ke,jb,je,nspecm
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: &
          bx,by,bz,den,ex,ey,ez,vix,viy,viz,eta
      real(kind=CUSTOM_REAL),                                         &
          dimension(nxmax,jb-1:je+1,kb-1:ke+1,nspecm) ::              &
          tpar,tperp
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  AddScalars - At this point we add raw scalars, otherwise we would 
!               have create new memory for the REAL*4, and scale by
!               other scalars and a real*4 normalizer.
!
!----------------------------------------------------------------------
        call adduniformgridscalar('bx',2,bx,                          &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('by',2,by,                          &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('bz',2,bz,                          &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('den',3,den,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('ex',2,ex,                          &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('ey',2,ey,                          &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('ez',2,ez,                          &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('vix',3,vix,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('viy',3,viy,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('viz',3,viz,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('tpar',4,tpar(:,:,:,1),             &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('tperp',5,tperp(:,:,:,1),           &
                                  nxmax*(nylmax+2)*(nzlmax+2))
        call adduniformgridscalar('eta',3,eta,                        &
                                  nxmax*(nylmax+2)*(nzlmax+2))
!----------------------------------------------------------------------
! end AddScalars
!----------------------------------------------------------------------
        return
      end subroutine addscalars
!
!***********************************************************************
!
      subroutine dataout( bx, by, bz, den, ex, ey, ez, vix, viy, viz,                             &
                          tpar, tperp, nxmax,nymax,nzmax,file_unit,myid,                          &
                          numvars,irecnum,kb,ke,numprocs,wpiwci,jb,je,ny,nz,nylmax,nzlmax,nspecm, &
                          bdipole_x,bdipole_y,bdipole_z,eta,eta_times_b_dot_j,eta_par,            &
                          uniform_mesh,nonuniform_mesh_global,data_directory, cycle_ascii,MPI_IO_format)
use PRECISION    !LAURA
      implicit none !LAURA
      integer*4 myid, numvars, numprocs
      integer*8 nxmax, nymax, nzmax, irecnum, kb, ke, jb,je,nylmax,nzlmax,nspecm
      integer*8, dimension(numvars+1) :: file_unit
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: bx, by, bz, den,                  &	!LAURA
      ex, ey, ez, vix, viy, viz,bdipole_x,bdipole_y,bdipole_z,eta,eta_times_b_dot_j
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:je+1,kb-1:ke+1,nspecm) :: tpar, tperp	!LAURA
      integer*8 ny,nz,eta_par
      real(kind=CUSTOM_REAL) rnorm, wpiwci	!LAURA
      integer*8 irecdel, ir1, ir2, irec_start , idebug,IERR
!      real(kind=CUSTOM_REAL):: uniform_mesh(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)	!LAURA
      real(kind=CUSTOM_REAL):: uniform_mesh(nxmax,jb-1:je+1,kb-1:ke+1)	!LAURA
      real(kind=CUSTOM_REAL):: nonuniform_mesh_global(nxmax,0:ny+1,0:nz+1)	!LAURA
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
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(1),irec_start,ny,nz)
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
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(2),irec_start,ny,nz)
        endif
        by=by-bdipole_y
        rnorm = wpiwci
        bz=bz+bdipole_z
!        call MESH_INTERPOLATED_3D(bz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=bz
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'bz_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(3),irec_start,ny,nz)
        endif
        bz=bz-bdipole_z
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(den,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=den
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'den_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(4),irec_start,ny,nz)
        endif
        rnorm = wpiwci**2
!        call MESH_INTERPOLATED_3D(ex,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=ex 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'ex_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(5),irec_start,ny,nz)
        endif
        rnorm = wpiwci**2
!        call MESH_INTERPOLATED_3D(ey,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=ey 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'ey_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(6),irec_start,ny,nz)
        endif
        rnorm = wpiwci**2
!        call MESH_INTERPOLATED_3D(ez,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=ez 
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'ez_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(7),irec_start,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(vix,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=vix
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'vix_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(8),irec_start,ny,nz)
        endif
        rnorm = wpiwci
!        call MESH_INTERPOLATED_3D(viy,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=viy
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'viy_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(9),irec_start,ny,nz)
        endif
        rnorm = wpiwci 
!        call MESH_INTERPOLATED_3D(viz,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=viz
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'viz_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(10),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(tpar,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=tpar(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'tpar_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(11),irec_start,ny,nz)
        endif
        rnorm = 1.
!        call MESH_INTERPOLATED_3D(tperp,uniform_mesh,nonuniform_mesh_global)
        uniform_mesh=tperp(:,:,:,1)
        if (MPI_IO_format) then
          fileName= trim(trim(adjustl(data_directory))//'tperp_'//trim(adjustl(cycle_ascii)))//'.gda'
!          call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
        else
!          call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(12),irec_start,ny,nz)
        endif
        if (eta_par == 0) then
          rnorm = 1.
!          call MESH_INTERPOLATED_3D(eta,uniform_mesh,nonuniform_mesh_global)
          uniform_mesh=eta
          if (MPI_IO_format) then
            fileName= trim(trim(adjustl(data_directory))//'eta_'//trim(adjustl(cycle_ascii)))//'.gda'
!            call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
          else
!            call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(13),irec_start,ny,nz)
          endif
        else
          rnorm = 1.
!          call MESH_INTERPOLATED_3D(eta_times_b_dot_j,uniform_mesh,nonuniform_mesh_global)
          uniform_mesh=eta_times_b_dot_j
          if (MPI_IO_format) then
            fileName= trim(trim(adjustl(data_directory))//'eta_par_'//trim(adjustl(cycle_ascii)))//'.gda'
!            call wrtfile         (uniform_mesh,rnorm,trim(adjustl(fileName)),irec_start,ny,nz)
          else
!            call wrtfile_NON_MPIO(uniform_mesh,rnorm,file_unit(13),irec_start,ny,nz)
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
      implicit none !LAURA
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
      implicit none !LAURA
      integer*8  file_unit_ref,lenrec !LAURA
      integer*8 j !LAURA
 
      if (.not.testorbt) then
!        lenrec=nxmax*recl_for_real
        lenrec=(nxmax-2)*recl_for_real
        file_unit_ref = 400
        do j=1,12
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
      endif
      return
      end subroutine openfiles
!
!***********************************************************************
!
      subroutine opendiagfiles
 
      use parameter_mod
      implicit none !LAURA
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
      implicit none !LAURA
      integer*8 f_unit,itstart,noresete,one !LAURA
      real*4 rw !LAURA
      one=1 !LAURA
 
      if (rw == +1.0_4) then
 
        t_stopped = t_stopped + (it-itstart+1)*dtwci
        f_unit=215+myid
!        open(unit=f_unit,file='restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        open(unit=f_unit,file=trim(adjustl(restart_directory))//'restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        ,form='unformatted',status='unknown')
      
        write(f_unit) x,y,z,vx,vy,vz,qp,link,porder

        write(f_unit) ninj,ninj_global,nescape,nescape_global,npart, &
        npart_global,qleft,qrite,maximum_simulation_time

        write(f_unit) x0,x1,tx0,vpar,vper,vbal,bbal,rcorr,teti,ishape

        write(f_unit) btspec, qspec, wspec, frac, vxbar, vybar,      &
        vzbar, anisot, denmin, resis, wpiwci, bete, fxsho,ave1,      &
        ave2,phib,theta, xmax,ymax,zmax,bxc,byc,bzc,gama,            &
        npx, npy, npz,                                               &
        iterb,norbskip,restrt_write,nxcel,netax,netay,netaz,nspec,   &
        nx,ny,nz,nskipx,nskipy,                                      &
        nskipz, testorbt, restart,etamin,etamax,ieta,eta_par,bz_IMF

        write(f_unit) hx,hy,hz,hxi,hyi,hzi                           &
       ,pi,efld,bfld,time,te0                                        &
       ,prntinfo,wrtdat,itfin,iwt                                    &
       ,nx1,nx2,ny1,ny2,nz1,nz2,it                                   &
       ,ipstore,nptot,npleaving,npentering,myid_stop                 &
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

        write(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp

        close(unit=f_unit)


      else if (rw == -1.0_4) then
 
        f_unit=215+myid
!        open(unit=f_unit,file='restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        open(unit=f_unit,file=trim(adjustl(restart_directory))//'restfld_'//trim(adjustl(myid_char))//'.bin'//restart_index_suffix(restart_index)&
        ,form='unformatted',status='unknown')
      
        read(f_unit) x,y,z,vx,vy,vz,qp,link,porder

        read(f_unit) ninj,ninj_global,nescape,nescape_global,npart,  &
        npart_global,qleft,qrite,t_stopped

        read(f_unit) x0,x1,tx0,vpar,vper,vbal,bbal,rcorr,teti,ishape

        read(f_unit) btspec, qspec, wspec, frac, vxbar, vybar,       &
        vzbar, anisot, denmin, resis, wpiwci, bete, fxsho,ave1,      &
        ave2,phib,theta, xmax,ymax,zmax,bxc,byc,bzc,gama,            &
        npx, npy, npz,                                               &
        iterb,norbskip,restrt_write,nxcel,netax,netay,netaz,nspec,   &
        nx,ny,nz,nskipx,nskipy,                                      &
        nskipz, testorbt, restart,etamin,etamax,ieta,eta_par,bz_IMF

        read(f_unit) hx,hy,hz,hxi,hyi,hzi                            &
       ,pi,efld,bfld,time,te0                                        &
       ,prntinfo,wrtdat,itfin,iwt                     &
       ,nx1,nx2,ny1,ny2,nz1,nz2,it                                   &
       ,ipstore,nptot,npleaving,npentering,myid_stop                 &
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

        read(f_unit) tpar,tperp,dns,vxs,vys,vzs,iphead,iptemp

        close(unit=f_unit)

! Reset electic field and fox,y,z 

         noresete = 1
         if (noresete == 0) then
           deno=den;vixo=vix;viyo=viy;vizo=viz
           call ecalc( one )
           call focalc
         endif

      endif
!
!***********************************************************************
!
      return
      end
!
!***********************************************************************
!
      subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
      implicit none !LAURA
 
      integer   numprocs,myid
      integer*8 n
      integer*8 s, e
      integer*8 nlocal
      integer*8 deficit
 
      nlocal  = n / numprocs
      s       = myid * nlocal + 1
      deficit = mod(int(n),int(numprocs))
      s       = s + min(int(myid),int(deficit))
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

      implicit none  !LAURA
      integer*8:: n,nskip,ncount,i !LAURA
      
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
      implicit none  !LAURA
      integer*8:: is,i,j,k,jbmin,jbmax,kbmin,kbmax
      integer*8:: time_begin(8),time_end(8)
      integer*8:: one,zero
      real(kind=CUSTOM_REAL):: dttmp,dns_tmp	!LAURA
      one=1 !LAURA
      zero=0 !LAURA
 
      call date_and_time(values=time_begin_array(:,20))
 
      do is=1,nspec
        DO K=KB-1,KE+1
          do j=jb-1,je+1
            do i=1,nx2
              dns(i,j,k,is)=1.e-10_CUSTOM_REAL
              vxs(i,j,k,is)=0._CUSTOM_REAL
              vys(i,j,k,is)=0._CUSTOM_REAL
              vzs(i,j,k,is)=0._CUSTOM_REAL
              if (is == 1) then
                deno(i,j,k)=den(i,j,k)
                vixo(i,j,k)=vix(i,j,k)
                viyo(i,j,k)=viy(i,j,k)
                vizo(i,j,k)=viz(i,j,k)
                den(i,j,k)=0._CUSTOM_REAL
                vix(i,j,k)=0._CUSTOM_REAL
                viy(i,j,k)=0._CUSTOM_REAL
                viz(i,j,k)=0._CUSTOM_REAL
              endif
            enddo
          enddo
        enddo
      enddo
 
      call date_and_time(values=time_end)
      clock_time=( time_end(5)*3600.+time_end(6)*60.+time_end(7)+time_end(8)*0.001)
 
      if (notime == 0) then
        write(file_unit_time,"(i4,' parmovin ',e20.8)") it,dble(clock_time-clock_time_init)
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
        write(file_unit_time,"(i4,' prmovout ',e20.8)") it,dble(clock_time-clock_time_init)
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
!              if (dns_tmp <= denmin) dns_tmp=1.E+10_CUSTOM_REAL	!LAURA

!             Nonuniform mesh
!              if (dns_tmp*cell_volume_ratio <= denmin) dns_tmp=1.E+10_CUSTOM_REAL	!LAURA
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
      call XREALBCC(PE ,one,NX,NY,NZ) !LAURA
      call XREALBCC(DEN,one,NX,NY,NZ) !LAURA
      call XREALBCC(VIX,one,NX,NY,NZ) !LAURA
      call XREALBCC(VIY,one,NX,NY,NZ) !LAURA
      call XREALBCC(VIZ,one,NX,NY,NZ) !LAURA
      else 
      call XREALBCC_2D(PE ,one,NX,NY,NZ) !LAURA
      call XREALBCC_2D(DEN,one,NX,NY,NZ) !LAURA
      call XREALBCC_2D(VIX,one,NX,NY,NZ) !LAURA
      call XREALBCC_2D(VIY,one,NX,NY,NZ) !LAURA
      call XREALBCC_2D(VIZ,one,NX,NY,NZ) !LAURA
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
      implicit none !LAURA
      integer*8 i,j,nx1m,ny1m,nz1m,k
      real(kind=CUSTOM_REAL) a(nxmax,jb-1:je+1,kb-1:ke+1),tmp(nxmax,jb-1:je+1,kb-1:ke+1)	!LAURA
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
      subroutine xreal_4d(a,is1m,nx1m,ny1m,nz1m,nspec1m)
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
      implicit none !LAURA
      integer*8 i,j,nx1m,ny1m,nz1m,k,is1m,nspec1m
      real(kind=CUSTOM_REAL) a(nxmax,jb-1:je+1,kb-1:ke+1,nspec1m),tmp(nxmax,jb-1:je+1,kb-1:ke+1,nspec1m)        !LAURA
      a(2   ,:,:,is1m)=a(2   ,:,:,is1m)+a(1   ,:,:,is1m)
      a(nx1m+1,:,:,is1m)=a(nx1m+1,:,:,is1m)+a(nx1m+2,:,:,is1m)
      call MPI_SENDRECV(a    (1    ,je+1,kb-1,is1m),1,stridery,nbrrite,0, &
                        tmp  (1    ,jb-1,kb-1,is1m),1,stridery,nbrleft,0, &
                        mpi_comm_world,status,ierr)
      if (jb == 1) tmp(:,jb-1,:,is1m)=a(:,jb-1,:,is1m)
      a(:,jb,:,is1m)=a(:,jb,:,is1m)+tmp  (:,jb-1,:,is1m)
      call MPI_SENDRECV(a    (1    ,jb-1,kb-1,is1m),1,stridery,nbrleft,1, &
                        tmp  (1    ,je+1,kb-1,is1m),1,stridery,nbrrite,1, &
                        mpi_comm_world,status,ierr)
      if (je == ny1m) tmp(:,je+1,:,is1m)=a(:,je+1,:,is1m)
      a(:,je,:,is1m)=a(:,je,:,is1m)+tmp  (:,je+1,:,is1m)
      call MPI_SENDRECV(a    (1    ,jb-1,ke+1,is1m),1,striderz,nbrtop ,0, &
                        tmp  (1    ,jb-1,kb-1,is1m),1,striderz,nbrbot ,0, &
                        mpi_comm_world,status,ierr)
      if (kb == 1) tmp(:,:,kb-1,is1m)=a(:,:,kb-1,is1m)
      a(:,:,kb,is1m)=a(:,:,kb,is1m)+tmp  (:,:,kb-1,is1m)
      call MPI_SENDRECV(a    (1    ,jb-1,kb-1,is1m),1,striderz,nbrbot ,1, &
                        tmp  (1    ,jb-1,ke+1,is1m),1,striderz,nbrtop ,1, &
                        mpi_comm_world,status,ierr)
      if (ke == nz1m) tmp(:,:,ke+1,is1m)=a(:,:,ke+1,is1m)
      a(:,:,ke,is1m)=a(:,:,ke,is1m)+tmp  (:,:,ke+1,is1m)


      a(1   ,:,:,is1m)=a(2   ,:,:,is1m)
      a(nx1m+2,:,:,is1m)=a(nx1m+1,:,:,is1m)


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
      implicit none !LAURA
      integer*8 ibnd,i,j,nx1m,ny1m,nz1m
      real(kind=CUSTOM_REAL) a(nxmax,jb-1:je+1,kb-1:ke+1),tmp(nxmax,jb-1:je+1,kb-1:ke+1)	!LAURA
 
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
      implicit none !LAURA
      real(KIND=CUSTOM_REAL) :: efldp,bfldp !LAURA
      integer*8 :: i,j,k !LAURA
 
      efldp=0._CUSTOM_REAL
      bfldp=0._CUSTOM_REAL
      do k=kb,ke
        do j=jb,je
          do i=2,nx1
            efldp=efldp+ex(i,j,k)**2_CUSTOM_REAL+ey(i,j,k)**2_CUSTOM_REAL+ez(i,j,k)**2_CUSTOM_REAL
            bfldp=bfldp+bx(i,j,k)**2_CUSTOM_REAL+by(i,j,k)**2_CUSTOM_REAL+bz(i,j,k)**2_CUSTOM_REAL
          enddo
        enddo
      enddo
      efldp=efldp*hx*hy*hz*0.5_CUSTOM_REAL
      bfldp=bfldp*hx*hy*hz*0.5_CUSTOM_REAL
      call MPI_ALLREDUCE(efldp,efld,1,CUSTOM_MPI_TYPE,&	!LAURA
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(bfldp,bfld,1,CUSTOM_MPI_TYPE,&	!LAURA
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      return
      end
!
!#######################################################################
!
      subroutine parmov
 
      use parameter_mod
      use MESH2D
      implicit none !LAURA
      real(KIND=CUSTOM_REAL) :: FOX1, FOX2, FOX3, FOX4, FOX5, FOX6  !LAURA
      real(KIND=CUSTOM_REAL) :: FOX7, FOX8, FOY1, FOY2, FOY3, FOY4  !LAURA
      real(KIND=CUSTOM_REAL) :: FOY5, FOY6, FOY7, FOY8, FOZ1, FOZ2  !LAURA
      real(KIND=CUSTOM_REAL) :: FOZ3, FOZ4, FOZ5, FOZ6, FOZ7, FOZ8  !LAURA
      real(KIND=CUSTOM_REAL) :: EPSILON, D_RANF, DELTIME1, DELTIME2, H, HH   !LAURA
      INTEGER*8 :: I, IIZE, IIYE, IIXE,ICOUNT, II, IIZ, IIY, IIX   !LAURA
      real(KIND=CUSTOM_REAL) :: EX1, EX2, EX3, EX4, EX5, EX6, EX7  !LAURA
      real(KIND=CUSTOM_REAL) :: EX8, EY1, EY2, EY3, EY4, EY5, EY6  !LAURA
      real(KIND=CUSTOM_REAL) :: EY7, EY8, EZ1, EZ2, EZ3, EZ4, EZ5  !LAURA
      real(KIND=CUSTOM_REAL) :: EZ6, EZ7, EZ8, BX1, BX2, BX3, BX4  !LAURA
      real(KIND=CUSTOM_REAL) :: BX5, BX6, BX7, BX8, BY1, BY2, BY3  !LAURA
      real(KIND=CUSTOM_REAL) :: BY4, BY5, BY6, BY7, BY8, BZ1, BZ2  !LAURA
      real(KIND=CUSTOM_REAL) :: BZ3, BZ4, BZ5, BZ6, BZ7, BZ8, EXA  !LAURA
      real(KIND=CUSTOM_REAL) :: EYA, EZA, BXA, BYA, BZA, FF  !LAURA
      real(KIND=CUSTOM_REAL) :: X_DIPOLE, Y_DIPOLE, Z_DIPOLE   !LAURA
      real(KIND=CUSTOM_REAL) :: WMULT, W1E, W2E, W3E, W4E, W5E, W6E   !LAURA
      real(KIND=CUSTOM_REAL) :: W7E, W8E, VEX, VEY, VEZ   !LAURA
      real(KIND=CUSTOM_REAL) :: P2XS, P2YS, P2ZS, XPART, YPART   !LAURA
      real(KIND=CUSTOM_REAL) :: ZPART, R_PARTICLE, VY_TMP, VMAG   !LAURA
      real(KIND=CUSTOM_REAL) :: TH, VZ_TMP, Q_P   !LAURA
      INTEGER*8 :: IYE_CC, IZE_CC, IREPEAT, NPRECVTMP, ITMP, NPRECV   !LAURA
      INTEGER*8 :: K, J, NPLEAVINGP, IS,IREPEATP, IV, JV, KSPC   !LAURA
      integer*8:: one !LAURA
      integer*8 count_kbq,time_begin(8),time_end(8)
      integer*8 nptotp_kbq,npart_kbq(2),np_ijk,Storage_Error_p,Storage_Error
      data fox1,fox2,fox3,fox4,fox5,fox6,fox7,fox8/0,0,0,0,0,0,0,0/
      data foy1,foy2,foy3,foy4,foy5,foy6,foy7,foy8/0,0,0,0,0,0,0,0/
      data foz1,foz2,foz3,foz4,foz5,foz6,foz7,foz8/0,0,0,0,0,0,0,0/
      integer*8:: nsendactual,nsendactualp,nrecvactualp,nrecvactual,jj,kk,ix,iy,iz,ixe,iye,ize           &
                 ,ixep1,iyep1,izep1,ixp1,iyp1,izp1
      real(kind=CUSTOM_REAL):: pdata(7),rx,ry,rz,fx,fy,fz,w1,w2,w3,w4,w5,w6,w7,w8	!LAURA
      real(kind=CUSTOM_REAL):: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi	!LAURA
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) :: bx_av,by_av,bz_av	!LAURA
      real(kind=CUSTOM_REAL):: v_limit,eps2,rx0,ry0,rz0,rrat,sqrr,outer_radius,myranf,twopi,fluxran,vxa,vyz,vza	!LAURA
      INTEGER*8:: L, EXIT_CODE_P, EXIT_CODE
      integer*8:: n_fast_removed,n_fast_removed_local,nptot_max,Courant_Violation,Courant_Violation_p,Field_Diverge,Field_Diverge_p
      real(kind=CUSTOM_REAL):: hxmin,hxmax,hymin,hymax,hzmin,hzmax,cell_size_min,x_disp,y_disp,z_disp          &	!LAURA
                        ,y_disp_max_p,x_disp_max_p,z_disp_max_p,y_disp_max,x_disp_max,z_disp_max
 
      one=1 !LAURA
      call date_and_time(values=time_begin_array(:,19))

      Storage_Error_p = 0
      Field_Diverge_p = 0
      Courant_Violation_p = 0
      x_disp_max_p        = 0
      y_disp_max_p        = 0
      z_disp_max_p        = 0

      dtxi = 1._CUSTOM_REAL/meshX%dt
      dtyi = 1._CUSTOM_REAL/meshY%dt
      dtzi = 1._CUSTOM_REAL/meshZ%dt

 
      epsilon= buffer_zone
      d_ranf=1._CUSTOM_REAL/1001._CUSTOM_REAL
      twopi=2._CUSTOM_REAL*dacos(-1._CUSTOM_REAL)

      outer_radius=dipole_sphere_radius+min(hx,hy,hz)
      eps2=1.E-25_CUSTOM_REAL	!LAURA
 
 
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
!        v_limit=1.E+10_CUSTOM_REAL	!LAURA
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
 
      bx_av=0._CUSTOM_REAL;by_av=0._CUSTOM_REAL;bz_av=0._CUSTOM_REAL
      DO K = KB-1,KE
        DO J = JB-1,JE
          DO I = 1, NX1
            bx_av(i,j,k)=0.125_CUSTOM_REAL*( bx(i  ,j  ,k  )+bdipole_x(i  ,j  ,k  )             &
                                +bx(i+1,j  ,k  )+bdipole_x(i+1,j  ,k  )             &
                                +bx(i  ,j+1,k  )+bdipole_x(i  ,j+1,k  )             &
                                +bx(i+1,j+1,k  )+bdipole_x(i+1,j+1,k  )             &
                                +bx(i  ,j  ,k+1)+bdipole_x(i  ,j  ,k+1)             &
                                +bx(i+1,j  ,k+1)+bdipole_x(i+1,j  ,k+1)             &
                                +bx(i  ,j+1,k+1)+bdipole_x(i  ,j+1,k+1)             &
                                +bx(i+1,j+1,k+1)+bdipole_x(i+1,j+1,k+1)             &
                               )
            by_av(i,j,k)=0.125_CUSTOM_REAL*( by(i  ,j  ,k  )+bdipole_y(i  ,j  ,k  )             &
                                +by(i+1,j  ,k  )+bdipole_y(i+1,j  ,k  )             &
                                +by(i  ,j+1,k  )+bdipole_y(i  ,j+1,k  )             &
                                +by(i+1,j+1,k  )+bdipole_y(i+1,j+1,k  )             &
                                +by(i  ,j  ,k+1)+bdipole_y(i  ,j  ,k+1)             &
                                +by(i+1,j  ,k+1)+bdipole_y(i+1,j  ,k+1)             &
                                +by(i  ,j+1,k+1)+bdipole_y(i  ,j+1,k+1)             &
                                +by(i+1,j+1,k+1)+bdipole_y(i+1,j+1,k+1)             &
                               )
            bz_av(i,j,k)=0.125_CUSTOM_REAL*( bz(i  ,j  ,k  )+bdipole_z(i  ,j  ,k  )             &
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
      CALL XREALBCC_PACK_B(BX_AV,BY_AV,BZ_AV,one,NX,NY,NZ)  !LAURA

      if (myid == 0) WRITE(6,*) " CALLING PARMOVE, NSPEC = ",NSPEC
      if (myid == 0) WRITE(6,*) " ABSORBING_DIPOLE       = ",absorbing_dipole
      if (myid == 0) WRITE(6,*) " DIPOLE_SPHERE_RADIUS   = ",dipole_sphere_radius
!
!=======================================================================
!
!  initalize diagnostic variables that keep track of
!  particle number, injection, and escape
!
      deltime1 = 0.0_CUSTOM_REAL
      deltime2 = 0.0_CUSTOM_REAL
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
        hh=.5_CUSTOM_REAL*h
 
 
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
!                  rxe=hxi*x(l)+1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                  rye=hyi*y(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                  rze=hzi*z(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
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
!                  rxe=hxi*x(l)+1.500000000000000E+00_CUSTOM_REAL	!LAURA
!                  rye=hyi*y(l)+1.500000000000000E+00_CUSTOM_REAL	!LAURA
!                  rze=hzi*z(l)+1.500000000000000E+00_CUSTOM_REAL	!LAURA
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
                  rxe=dtxi*MESH_UNMAP(meshX,x(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                  rye=dtyi*MESH_UNMAP(meshY,y(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                  rze=dtzi*MESH_UNMAP(meshZ,z(l))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
 

                  w1e=(1._CUSTOM_REAL-fxe)*(1._CUSTOM_REAL-fye)*(1._CUSTOM_REAL-fze)
                  w2e=fxe*(1._CUSTOM_REAL-fye)*(1._CUSTOM_REAL-fze)
                  w3e=(1._CUSTOM_REAL-fxe)*fye*(1._CUSTOM_REAL-fze)
                  w4e=fxe*fye*(1._CUSTOM_REAL-fze)
                  w5e=(1._CUSTOM_REAL-fxe)*(1._CUSTOM_REAL-fye)*fze
                  w6e=fxe*(1._CUSTOM_REAL-fye)*fze
                  w7e=(1._CUSTOM_REAL-fxe)*fye*fze
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
 
                  ff=2._CUSTOM_REAL/(1._CUSTOM_REAL+hh*hh*(bxa**2_CUSTOM_REAL+bya**2_CUSTOM_REAL+bza**2_CUSTOM_REAL))
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

        call MPI_ALLREDUCE(x_disp_max_p,x_disp_max,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,IERR)	!LAURA
        call MPI_ALLREDUCE(y_disp_max_p,y_disp_max,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,IERR)	!LAURA
        call MPI_ALLREDUCE(z_disp_max_p,z_disp_max,1,CUSTOM_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD,IERR)	!LAURA
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
                  r_particle=dsqrt((xpart-x_dipole)**2_CUSTOM_REAL+(ypart-y_dipole)**2_CUSTOM_REAL+(zpart-z_dipole)**2_CUSTOM_REAL)
                else
                  r_particle=dsqrt((xpart-x_dipole)**2_CUSTOM_REAL+(ypart-y_dipole)**2_CUSTOM_REAL)
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
                      rrat = (2.0_CUSTOM_REAL *dipole_sphere_radius -sqrr) /sqrr
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
                    
                    y(l) = 2._CUSTOM_REAL*epsilon-y(l)

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vy_tmp=vpar(is)*u_array_injection(iv,jv)
                    vy(l)=+vy_tmp
                    vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*dcos(th)
                    vz(l)=          vmag*vper(is)*dsin(th)
                    x (l)=   myranf()*xmax
                    z (l)=zb+myranf()*(ze-zb)

                  endif

                  if       (y(l) >   ymax-epsilon   ) then

                    y(l) = 2._CUSTOM_REAL*(ymax-epsilon) - y(l)

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vy_tmp=vpar(is)*u_array_injection(iv,jv)
                    vy(l)=-vy_tmp
                    vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*dcos(th)
                    vz(l)=          vmag*vper(is)*dsin(th)
                    x (l)=     myranf()*xmax
                    z (l)=zb  +myranf()*(ze-zb)

                  endif

                  if       (z(l) <   epsilon     ) then
                    
                    z(l) = 2._CUSTOM_REAL*epsilon - z(l)

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vz_tmp=vpar(is)*u_array_injection(iv,jv)
                    vz(l)=+vz_tmp
                    vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*dcos(th)
                    vy(l)=          vmag*vper(is)*dsin(th)
                    x (l)=   myranf()*xmax
                    y (l)=yb+myranf()*(ye-yb)

                  endif

                  if       (z(l) >   zmax-epsilon   ) then

                    z(l) = 2._CUSTOM_REAL*(zmax-epsilon) - z(l)

                    fluxran=myranf()
                    iv=int((fluxran/d_ranf)+0.5)
                    iv=max(1,min(1001,iv))
                    jv=0
                    vz_tmp=vpar(is)*u_array_injection(iv,jv)
                    vz(l)=-vz_tmp
                    vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*myranf()))
                    th=twopi*myranf()
                    vx(l)=vxbar(is)+vmag*vper(is)*dcos(th)
                    vy(l)=          vmag*vper(is)*dsin(th)
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
        call MPI_ALLREDUCE(nsendtotp,nsendtot,1,MPI_INTEGER8,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(nrecvtotp,nrecvtot,1,MPI_INTEGER8,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
        call MPI_ALLREDUCE(n_fast_removed_local,n_fast_removed,1,MPI_INTEGER8,MPI_SUM,&
                           MPI_COMM_WORLD,IERR)
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
              call MPI_SEND(pdata,7,CUSTOM_MPI_TYPE,&	!LAURA
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
              call MPI_RECV(pdata,7,CUSTOM_MPI_TYPE,&	!LAURA
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
!              ixe=hxi*x(nprecv)+1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!              iye=hyi*y(nprecv)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!              ize=hzi*z(nprecv)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA

!             Nonuniform mesh - using MESH_UNMAP
              rxe=dtxi*MESH_UNMAP(meshX,x(nprecv))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              rye=dtyi*MESH_UNMAP(meshY,y(nprecv))+1.50000000000E+00_CUSTOM_REAL	!LAURA
              rze=dtzi*MESH_UNMAP(meshZ,z(nprecv))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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

        call MPI_ALLREDUCE(Storage_Error_p,Storage_Error,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
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
        if (myid == 0) then
          write(6,*) " FINISHED EXCHANGING PARTICLES "
          write(6,*) " # OF PARTICLES       SENT     = ",NSENDACTUAL
          write(6,*) " # OF PARTICLES       RECEIVED = ",NRECVACTUAL
        endif
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
!                rx=hxi*x(l)+1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                ry=hyi*y(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                rz=hzi*z(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
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
!                rx=hxi*x(l)+1.500000000000000E+00_CUSTOM_REAL	!LAURA
!                ry=hyi*y(l)+1.500000000000000E+00_CUSTOM_REAL	!LAURA
!                rz=hzi*z(l)+1.500000000000000E+00_CUSTOM_REAL	!LAURA
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
 
                w1=q_p*(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)*(1._CUSTOM_REAL-fz)
                w2=q_p*fx     *(1._CUSTOM_REAL-fy)*(1._CUSTOM_REAL-fz)
                w3=q_p*(1._CUSTOM_REAL-fx)*fy     *(1._CUSTOM_REAL-fz)
                w4=q_p*fx     *fy     *(1._CUSTOM_REAL-fz)
                w5=q_p*(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)*fz
                w6=q_p*fx     *(1._CUSTOM_REAL-fy)*fz
                w7=q_p*(1._CUSTOM_REAL-fx)*fy     *fz
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

        call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape(is),nescape_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape_yz(is),nescape_yz_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape_zy(is),nescape_zy_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape_xy(is),nescape_xy_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape_yx(is),nescape_yx_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape_zx(is),nescape_zx_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)
        call MPI_ALLREDUCE(nescape_xz(is),nescape_xz_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,COMM2D,IERR)

        deltime2 = deltime2 + dble(clock_time1-clock_time)
 
        call XREAL(DNS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(VXS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(VYS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL(VZS(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREALBCC(DNS(1,jb-1,kb-1,is),one,NX,NY,NZ) !LAURA
        call XREALBCC(VXS(1,jb-1,kb-1,is),one,NX,NY,NZ) !LAURA
        call XREALBCC(VYS(1,jb-1,kb-1,is),one,NX,NY,NZ) !LAURA
        call XREALBCC(VZS(1,jb-1,kb-1,is),one,NX,NY,NZ) !LAURA


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
      call MPI_ALLREDUCE(nptotp,nptot,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
      call MPI_ALLREDUCE(npleavingp,npleaving,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)

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
      implicit none !LAURA
      real(kind=CUSTOM_REAL) pstore(nplmax)	!LAURA
      integer*8 id, kb1, is, ix, iy, iz, ixe ,iye, ize, l, nttot, nplist
      real(kind=CUSTOM_REAL):: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi	!LAURA
 
      dtxi = 1._CUSTOM_REAL/meshX%dt
      dtyi = 1._CUSTOM_REAL/meshY%dt
      dtzi = 1._CUSTOM_REAL/meshZ%dt
 
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
!                ixe = hxi*x(np)+1.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                iye = hyi*y(np)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                ize = hzi*z(np)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA

!                 Nonuniform mesh - using MESH_UNMAP
                  rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                  rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
                  rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
      implicit none !LAURA
      integer*8 iunit,i2,i1,is,it !LAURA
      character*10 message
      write(iunit,"(i4,a10,e20.8)") it, dble(i2-i1)/dble(is)
      return
      end

!
!#######################################################################
!
      subroutine field
 
      use parameter_mod
      implicit none !LAURA
      integer*8 :: one,zero !LAURA
      one=1 !LAURA
      zero=0 !LAURA
  
!     Bypass field solver for timing on LANL machine
!      goto 999 


      call date_and_time(values=time_begin_array(:,21))
 
 
      call date_and_time(values=time_begin_array(:,9))
      call pressgrad(one)
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
      call pressgrad(zero)
      call date_and_time(values=time_end_array(:,9))
      call accumulate_time_difference(time_begin_array(1,9) &
     &                               ,time_end_array(1,9) &
     &                               ,time_elapsed(9))
 
 
      call date_and_time(values=time_begin_array(:,11))
      call ecalc( zero )
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
      implicit none !LAURA
      integer*8 :: iflag,i,k,j !LAURA
      real(kind=CUSTOM_REAL):: dena,a,dxa,dza,dya !LAURA

      do k=kb,ke
        do j = jb,je
          do i=2,nx1
            dena=dble(iflag)*0.5_CUSTOM_REAL*(den(i,j,k)+deno(i,j,k))&
                +(1._CUSTOM_REAL-dble(iflag))*den(i,j,k)
            a=1._CUSTOM_REAL/dena

!           Uniform mesh - Same as in version 5.0
!            dxa=a/(4._CUSTOM_REAL*hx)
!            dya=a/(4._CUSTOM_REAL*hy)
!            dza=a/(4._CUSTOM_REAL*hz)

!           Nonuniform mesh
            dxa=a/(2._CUSTOM_REAL*(meshX%dxn(i  )+meshX%dxn(i+1)))
            dya=a/(2._CUSTOM_REAL*(meshY%dxn(j+1)+meshY%dxn(j+2)))  ! integer index in y direction starts at 0
            dza=a/(2._CUSTOM_REAL*(meshZ%dxn(k+1)+meshZ%dxn(k+2)))  ! integer index in z direction starts at 0


           dpedx(i,j,k)=((pe(i+1,j-1,k+1)+2._CUSTOM_REAL*pe(i+1,j,k+1)&
               +pe(i+1,j+1,k+1))/4._CUSTOM_REAL&
               +2._CUSTOM_REAL*(pe(i+1,j-1,k  )+2._CUSTOM_REAL*pe(i+1,j,k  )+pe(i+1,j+1,k  ))/4._CUSTOM_REAL&
               +   (pe(i+1,j-1,k-1)+2._CUSTOM_REAL*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4._CUSTOM_REAL&
               -   (pe(i-1,j-1,k+1)+2._CUSTOM_REAL*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4._CUSTOM_REAL&
               -2._CUSTOM_REAL*(pe(i-1,j-1,k  )+2._CUSTOM_REAL*pe(i-1,j,k  )+pe(i-1,j+1,k  ))/4._CUSTOM_REAL&
               -   (pe(i-1,j-1,k-1)+2._CUSTOM_REAL*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4._CUSTOM_REAL&
                     )*dxa
             dpedy(i,j,k)=((pe(i-1,j+1,k+1)+2._CUSTOM_REAL*pe(i,j+1,k+1)&
                +pe(i+1,j+1,k+1))/4._CUSTOM_REAL&
                +2._CUSTOM_REAL*(pe(i-1,j+1,k  )+2._CUSTOM_REAL*pe(i,j+1,k  )+pe(i+1,j+1,k  ))/4._CUSTOM_REAL&
                +   (pe(i-1,j+1,k-1)+2._CUSTOM_REAL*pe(i,j+1,k-1)+pe(i+1,j+1,k-1))/4._CUSTOM_REAL&
                -   (pe(i-1,j-1,k+1)+2._CUSTOM_REAL*pe(i,j-1,k+1)+pe(i+1,j-1,k+1))/4._CUSTOM_REAL&
                -2._CUSTOM_REAL*(pe(i-1,j-1,k  )+2._CUSTOM_REAL*pe(i,j-1,k  )+pe(i+1,j-1,k  ))/4._CUSTOM_REAL&
                -   (pe(i-1,j-1,k-1)+2._CUSTOM_REAL*pe(i,j-1,k-1)+pe(i+1,j-1,k-1))/4._CUSTOM_REAL&
                     )*dya
             dpedz(i,j,k)=((pe(i+1,j-1,k+1)+2._CUSTOM_REAL*pe(i+1,j,k+1)&
                +pe(i+1,j+1,k+1))/4._CUSTOM_REAL&
                +2._CUSTOM_REAL*(pe(i  ,j-1,k+1)+2._CUSTOM_REAL*pe(i  ,j,k+1)+pe(i  ,j+1,k+1))/4._CUSTOM_REAL&
                +   (pe(i-1,j-1,k+1)+2._CUSTOM_REAL*pe(i-1,j,k+1)+pe(i-1,j+1,k+1))/4._CUSTOM_REAL&
                -   (pe(i+1,j-1,k-1)+2._CUSTOM_REAL*pe(i+1,j,k-1)+pe(i+1,j+1,k-1))/4._CUSTOM_REAL&
                -2._CUSTOM_REAL*(pe(i  ,j-1,k-1)+2._CUSTOM_REAL*pe(i  ,j,k-1)+pe(i  ,j+1,k-1))/4._CUSTOM_REAL&
                -   (pe(i-1,j-1,k-1)+2._CUSTOM_REAL*pe(i-1,j,k-1)+pe(i-1,j+1,k-1))/4._CUSTOM_REAL&
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
      implicit none !LAURA
      integer*8 :: one,iflag,i,k,j !LAURA
      real(kind=CUSTOM_REAL):: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &	!LAURA
                        ,curr_tot
      real(kind=CUSTOM_REAL):: bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8  !LAURA
      real(kind=CUSTOM_REAL):: by1,by2,by3,by4,by5,by6,by7,by8  !LAURA
      real(kind=CUSTOM_REAL):: bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8  !LAURA
      real(kind=CUSTOM_REAL):: vixa, viya, viza, dena, a, dxa, dya, dza  !LAURA
      real(kind=CUSTOM_REAL):: dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy  !LAURA
      real(kind=CUSTOM_REAL):: curlbx_scalar, curlby_scalar, curlbz_scalar, bxav, byav, bzav  !LAURA
      real(kind=CUSTOM_REAL):: dexdy, dexdz, deydx, deydz, dezdx,dezdy  !LAURA

      one=1 !LAURA
 
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
 
            vixa=(1._CUSTOM_REAL-dble(iflag))*(1.5_CUSTOM_REAL*vix(i,j,k)-0.5_CUSTOM_REAL*vixo(i,j,k))&
                +dble(iflag)*vix(i,j,k)
            viya=(1._CUSTOM_REAL-dble(iflag))*(1.5_CUSTOM_REAL*viy(i,j,k)-0.5_CUSTOM_REAL*viyo(i,j,k))&
                +dble(iflag)*viy(i,j,k)
            viza=(1._CUSTOM_REAL-dble(iflag))*(1.5_CUSTOM_REAL*viz(i,j,k)-0.5_CUSTOM_REAL*vizo(i,j,k))&
                +dble(iflag)*viz(i,j,k)

            dena=dble(iflag)*0.5_CUSTOM_REAL*(den(i,j,k)+deno(i,j,k))&
                +(1._CUSTOM_REAL-dble(iflag))*den(i,j,k)
            a=1._CUSTOM_REAL/dena

!           Uniform mesh - Same as is in version 5.0
!            dxa=a/(4.*hx)
!            dya=a/(4.*hy)
!            dza=a/(4.*hz)

!          Nonuniform mesh
            dxa=a/(4._CUSTOM_REAL*meshX%dxc(i))
            dya=a/(4._CUSTOM_REAL*meshY%dxc(j+1))                  ! integer index in y direction starts at 0
            dza=a/(4._CUSTOM_REAL*meshZ%dxc(k+1))                  ! integer index in z direction starts at 0

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
            bxav=.125_CUSTOM_REAL*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125_CUSTOM_REAL*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125_CUSTOM_REAL*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
!
! 6/25/2006 New eta_par option: tensor eta
!
            xj = curlbx_scalar
            yj = curlby_scalar
            zj = curlbz_scalar
            if (eta_par == 0) then
		    tenx=eta(i,j,k)*xj
		    teny=eta(i,j,k)*yj
		    tenz=eta(i,j,k)*zj
		 else
		    bxx = bxav
		    byy = byav
		    bzz = bzav
		    btot = dsqrt(bxx**2_CUSTOM_REAL + byy**2_CUSTOM_REAL + bzz**2_CUSTOM_REAL)
		 
		    if (eta_par == 1) then

		       tjdotb = eta(i,j,k)*(bxx*xj + byy*yj + bzz*zj)/btot
  
               tenx = tjdotb*bxx/btot
		       teny = tjdotb*byy/btot
		       tenz = tjdotb*bzz/btot

		    else if (eta_par == 2) then

		       curr_tot = max(1.E-12_CUSTOM_REAL,dsqrt(xj**2_CUSTOM_REAL + yj**2_CUSTOM_REAL + zj**2_CUSTOM_REAL))	!LAURA
  
               tenx = dabs(eta(i,j,k)*bxx*xj/(btot*curr_tot))
		       tenx = min(resis,tenx)
		       tenx = tenx*xj
               teny = dabs(eta(i,j,k)*byy*yj/(btot*curr_tot))
		       teny = min(resis,teny)
		       teny = teny*yj
               tenz = dabs(eta(i,j,k)*bzz*zj/(btot*curr_tot))
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
      call XREALBCC_PACK_E(EX,EY,EZ,one,NX,NY,NZ)
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
            curlex(i,j,k)=dezdy/(4._CUSTOM_REAL*meshY%dxn(j+1))-deydz/(4._CUSTOM_REAL*meshZ%dxn(k+1))        ! integer index in y and z directions start  at 0
            curley(i,j,k)=dexdz/(4._CUSTOM_REAL*meshZ%dxn(k+1))-dezdx/(4._CUSTOM_REAL*meshX%dxn(i  ))        ! integer index in z       direction  starts at 0
            curlez(i,j,k)=deydx/(4._CUSTOM_REAL*meshX%dxn(i  ))-dexdy/(4._CUSTOM_REAL*meshY%dxn(j+1))        ! integer index in y       direction  starts at 0

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
      implicit none !LAURA
      real(kind=CUSTOM_REAL):: tempx1(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
                        ,tempy1(nxmax,jb-1:je+1,kb-1:ke+1)&
                        ,tempz1(nxmax,jb-1:je+1,kb-1:ke+1)
      real(kind=CUSTOM_REAL):: dts,dts2,dts6  !LAURA
      integer*8:: ii,k,j,i
      integer*8:: one !LAURA
      one=1  !LAURA
 
 
       call date_and_time(values=time_begin_array(:,22))
 
       dts=dt/dble(iterb)
       dts2=dts/2._CUSTOM_REAL
       dts6=dts/6._CUSTOM_REAL
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
       call ecalc( one )  !LAURA
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
       call ecalc( one )  !LAURA
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
               tempx1(i,j,k)=tempx1(i,j,k)+2._CUSTOM_REAL*curlex(i,j,k)
               tempy1(i,j,k)=tempy1(i,j,k)+2._CUSTOM_REAL*curley(i,j,k)
               tempz1(i,j,k)=tempz1(i,j,k)+2._CUSTOM_REAL*curlez(i,j,k)
             enddo
           enddo
         enddo

!*****************
!  R-K  part 3
!*****************

       call date_and_time(values=time_begin_array(:,16))
       call ecalc( one )  !LAURA
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
               tempx1(i,j,k)=tempx1(i,j,k)+2._CUSTOM_REAL*curlex(i,j,k)
               tempy1(i,j,k)=tempy1(i,j,k)+2._CUSTOM_REAL*curley(i,j,k)
               tempz1(i,j,k)=tempz1(i,j,k)+2._CUSTOM_REAL*curlez(i,j,k)
             enddo
           enddo
         enddo


!***************
!   R-K  part 4
!***************
!
       call date_and_time(values=time_begin_array(:,16))
       call ecalc( one ) !LAURA
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
         CALL XREALBCC_PACK_B(BX,BY,BZ,one,NX,NY,NZ)  !LAURA
 
 
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
      implicit none !LAURA
      real(kind=CUSTOM_REAL):: tenx,teny,tenz,xj,yj,zj,bxx,byy,bzz,btot,tjdotb &	!LAURA
                        ,curr_tot
      integer*8 :: i,j,k,one  !LAURA
      real(kind=CUSTOM_REAL):: dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy !LAURA
      real(kind=CUSTOM_REAL):: curlbx_scalar, curlby_scalar, curlbz_scalar  !LAURA
      real(kind=CUSTOM_REAL):: bx1, bx2, bx3, bx4, bx5, bx6, bx7, bx8  !LAURA
      real(kind=CUSTOM_REAL):: by1, by2, by3, by4, by5, by6, by7, by8  !LAURA
      real(kind=CUSTOM_REAL):: bz1, bz2, bz3, bz4, bz5, bz6, bz7, bz8  !LAURA
      real(kind=CUSTOM_REAL):: bxav,byav,bzav  !LAURA
      one = 1
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
            curlbx_scalar=dbzdy/(4._CUSTOM_REAL*meshY%dxc(j+1))-dbydz/(4._CUSTOM_REAL*meshZ%dxc(k+1))
            curlby_scalar=dbxdz/(4._CUSTOM_REAL*meshZ%dxc(k+1))-dbzdx/(4._CUSTOM_REAL*meshX%dxc(i  ))
            curlbz_scalar=dbydx/(4._CUSTOM_REAL*meshX%dxc(i  ))-dbxdy/(4._CUSTOM_REAL*meshY%dxc(j+1))

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
            bxav=.125_CUSTOM_REAL*(bx1+bx2+bx3+bx4+bx5+bx6+bx7+bx8)
            byav=.125_CUSTOM_REAL*(by1+by2+by3+by4+by5+by6+by7+by8)
            bzav=.125_CUSTOM_REAL*(bz1+bz2+bz3+bz4+bz5+bz6+bz7+bz8)
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
      call XREALBCC_PACK_E(fox,foy,foz,one,NX,NY,NZ)


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
      implicit none !LAURA
      integer*8 :: i,j,k  !LAURA
      real(kind=CUSTOM_REAL), dimension(nxmax,jb-1:je+1,kb-1:ke+1) :: temp, a	!LAURA
      integer*8 :: zero  !LAURA
      zero=0  !LAURA

!  smoothing routine--assumes aperiodic in x
 
      call XREALBCC(a,zero,NX,NY,NZ)  !LAURA
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
            a(i,j,k)=temp(i,j,k)/8._CUSTOM_REAL&
      +( temp(i-1,j,k)+temp(i+1,j,k)+temp(i,j+1,k)+temp(i,j-1,k)&
      +temp(i,j,k+1)+temp(i,j,k-1))/16._CUSTOM_REAL&
      +( temp(i+1,j+1,k)+temp(i+1,j-1,k)+temp(i-1,j+1,k)&
      +temp(i-1,j-1,k)&
      +temp(i,j+1,k+1)+temp(i,j-1,k+1)+temp(i,j+1,k-1)+temp(i,j-1,k-1)&
      +temp(i+1,j,k+1)+temp(i-1,j,k+1)+temp(i+1,j,k-1)&
      +temp(i-1,j,k-1))/32._CUSTOM_REAL&
      +( temp(i+1,j+1,k+1)+temp(i-1,j+1,k+1)&
      +temp(i+1,j-1,k+1)+temp(i-1,j-1,k+1)&
      +temp(i+1,j+1,k-1)+temp(i-1,j+1,k-1)&
      +temp(i+1,j-1,k-1)+temp(i-1,j-1,k-1))/64._CUSTOM_REAL
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

      call XREALBCC(a,zero,NX,NY,NZ)  !LAURA

      return
      end
!
!#######################################################################
!
!
      FUNCTION myranf()
 
use PRECISION    !LAURA
      implicit none
      INTEGER*8:: idum
      INTEGER*8:: MBIG,MSEED,MZ
      real(kind=CUSTOM_REAL):: myranf,FAC	!LAURA
      INTEGER*8:: i,iff,ii,inext,inextp,k
      INTEGER*8:: mj,mk,ma(55)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0 &
                ,FAC=1.E-09_CUSTOM_REAL)	!LAURA
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      common /myrandom/ idum
 
      if (idum < 0.or.iff == 0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(int(21*i),55)
          ma(ii)=mk
          mk=mj-mk
          if (mk < MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(int(i+30),55))
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
use PRECISION    !LAURA
      implicit none
      integer*8:: ndatum,f_unit
      integer:: ndatum_4
      real(kind=CUSTOM_REAL):: datum(ndatum),datum_tmp(ndatum)	!LAURA
 
      ndatum_4=ndatum
      do ipe=0,numprocs-1
        if (myid == 0) then
          if (ipe == 0) then
            write(f_unit) datum
          else
              call MPI_IRECV(datum_tmp,ndatum_4,CUSTOM_MPI_TYPE&	!LAURA
                            ,ipe,0,MPI_COMM_WORLD,req(1),IERR)
              CALL MPI_WAIT(REQ(1),STATUS1,IERR)
              write(f_unit) datum_tmp
          endif
        else
          if (myid == ipe) then
            call MPI_ISEND(datum,ndatum_4,CUSTOM_MPI_TYPE&	!LAURA
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
use PRECISION    !LAURA
      implicit none
      integer*8:: ndatum,f_unit
      integer:: ndatum_4
      real(kind=CUSTOM_REAL):: datum(ndatum),datum_tmp(ndatum)	!LAURA
 
 
      ndatum_4=ndatum
      do ipe=0,numprocs-1
        if (myid == 0) then
          if (ipe == 0) then
            read(f_unit) datum
          else
              read(f_unit) datum_tmp
              call MPI_ISEND(datum_tmp,ndatum_4,CUSTOM_MPI_TYPE&	!LAURA
                            ,ipe,0,MPI_COMM_WORLD,req(1),IERR)
              CALL MPI_WAIT(REQ(1),STATUS1,IERR)
          endif
        else
          if (myid == ipe) then
            call MPI_IRECV(datum,ndatum_4,CUSTOM_MPI_TYPE&	!LAURA
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
 
use PRECISION    !LAURA
      implicit none
      integer,dimension(8):: time_begin,time_end
      real(kind=CUSTOM_REAL):: time_elapsed	!LAURA
 
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
      implicit none !LAURA
      integer*8:: random_init,ibp1,ibp2,nptot_max,i,remake
      integer*8 :: one !LAURA
      real(kind=CUSTOM_REAL):: my_random_init,myranf,by_avg,bz_avg,q_p	!LAURA
      real(kind=CUSTOM_REAL):: rxe,rye,rze,fxe,fye,fze,dtxi,dtyi,dtzi,profile_factor   &	!LAURA
                        ,b_dipole_sqr,x_p,y_p,z_p,x_p_logical,y_p_logical        &
                        ,z_p_logical,r_c

      real(kind=CUSTOM_REAL):: ranval, pifac, x_dipole, y_dipole, z_dipole, &
                                curion, phibr, phicr,r_particle, vthfac, vcurr  !LAURA
      real(kind=CUSTOM_REAL):: xran,vmag,th,vxa,vya,vza,bx_avg,dtsav,field_subcycle !LAURA
      integer*8 :: j,is, ixe, iye, ize, nxa, nxcelb, nxb, nxc, k, ipb1, ipb2, isp, nk, ip  !LAURA

      one=1  !LAURA
 

      remake = 0

      dtxi = 1._CUSTOM_REAL/meshX%dt
      dtyi = 1._CUSTOM_REAL/meshY%dt
      dtzi = 1._CUSTOM_REAL/meshZ%dt

 
      random_init=5000000*myid+1000
      IDUM=0
      do i=1,random_init
        my_random_init=myranf()
      enddo
 
!  initialize the random number generator with a different
!  seed for different processors
 
      ranval=myranf()
      pi=dacos(-1.E+00_CUSTOM_REAL)	!LAURA
      pifac=180._CUSTOM_REAL/pi
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
      hxi=1._CUSTOM_REAL/hx
      hyi=1._CUSTOM_REAL/hy
      hzi=1._CUSTOM_REAL/hz
 
 
 
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

        xb        = 0._CUSTOM_REAL
        xe        = xmax
        xb_logical=MESH_UNMAP(meshX,xb)
        xe_logical=MESH_UNMAP(meshX,xe)
        yb_logical=MESH_UNMAP(meshY,yb)
        ye_logical=MESH_UNMAP(meshY,ye)
        zb_logical=MESH_UNMAP(meshZ,zb)
        ze_logical=MESH_UNMAP(meshZ,ze)


        do is=1,nspec
          npm=npx(is)*npy(is)*npz(is)*npes
          dfac(is)=dble(ny*nz*nx)/dble(npm)
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
 
      curion=0.5_CUSTOM_REAL
      phibr=phib/pifac
      phicr = 0.5_CUSTOM_REAL *(180.0_CUSTOM_REAL-phib)/ pifac
      do is=1,nspec
        ninj(is)=0
        ninj_global(is)=0
        npart(is)=0
        tx0(is)=btspec(is)/(2._CUSTOM_REAL*wpiwci**2)
        x0(is)=0._CUSTOM_REAL
        x1(is)=xmax
        if (is == 1) x1(is)=fxsho*xmax
        if (is == 2) x0(is)=fxsho*xmax
        call MPI_ALLREDUCE(npart(is),npart_global(is),1,MPI_INTEGER8,&
                         MPI_SUM,MPI_COMM_WORLD,IERR)
      enddo
      te0=bete/(2._CUSTOM_REAL*wpiwci**2)
      vbal=1._CUSTOM_REAL

      if (b_in_xy) then
        bxc = dcos(theta/pifac)/wpiwci
        byc = dsin(theta/pifac)/wpiwci
        bzc = 0._CUSTOM_REAL
      else
        bxc = dcos(theta/pifac)/wpiwci
        byc = 0._CUSTOM_REAL
        bzc = dsin(theta/pifac)/wpiwci
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
        dfac(is)=dble(ny*nz)*(x1(is)-x0(is))/(hx*dble(npm))
        vpar(is)=dsqrt(btspec(is)/(wspec(is)*wpiwci*wpiwci))
        vper(is)=vpar(is)*dsqrt(anisot(is))

        do ip=ipb1,ipb2
            
          if (ipb2 == ipb1) write(6,*) "myid = , # particles = ",myid,ipb1,ipb2

          if (uniform_loading_in_logical_grid) then
            X_P_LOGICAL  = XB_LOGICAL+(XE_LOGICAL-XB_LOGICAL)*myranf()
            Y_P_LOGICAL  = YB_LOGICAL+(YE_LOGICAL-YB_LOGICAL)*myranf()
            Z_P_LOGICAL  = ZB_LOGICAL+(ZE_LOGICAL-ZB_LOGICAL)*myranf()
            X_P          = MESH_MAP(meshX,X_P_LOGICAL)
            Y_P          = MESH_MAP(meshY,Y_P_LOGICAL)
            Z_P          = MESH_MAP(meshZ,Z_P_LOGICAL)
            IXE          = dtxi*X_P_LOGICAL+1.50000000000E+00_CUSTOM_REAL	!LAURA
            IYE          = dtyi*Y_P_LOGICAL+1.50000000000E+00_CUSTOM_REAL	!LAURA
            IZE          = dtzi*Z_P_LOGICAL+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
            r_particle=dsqrt((x_p-x_dipole)**2_CUSTOM_REAL+(y_p-y_dipole)**2_CUSTOM_REAL+(z_p-z_dipole)**2_CUSTOM_REAL)
          else
            r_particle=dsqrt((x_p-x_dipole)**2_CUSTOM_REAL+(y_p-y_dipole)**2_CUSTOM_REAL)
          endif

          if (r_particle < dipole_sphere_radius) goto 10
 
          np     = ipstore
          x(np)  = x_p
          y(np)  = y_p
          z(np)  = z_p
          qp(np) = q_p


!         Nonuniform mesh - using MESH_UNMAP
          rxe=dtxi*MESH_UNMAP(meshX,x(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
          rye=dtyi*MESH_UNMAP(meshY,y(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
          rze=dtzi*MESH_UNMAP(meshZ,z(np))+1.50000000000E+00_CUSTOM_REAL	!LAURA
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
              profile_factor = 1._CUSTOM_REAL
            else

              if (ndim /= 1) then
                r_c = dsqrt((meshX%xc(ixe)-x_dipole)**2_CUSTOM_REAL+(meshY%xc(iye+1)-y_dipole)**2_CUSTOM_REAL+(meshZ%xc(ize+1)-z_dipole)**2_CUSTOM_REAL)
              else
                r_c = dsqrt((meshX%xc(ixe)-x_dipole)**2_CUSTOM_REAL+(meshY%xc(iye+1)-y_dipole)**2_CUSTOM_REAL)
              endif

              if (r_c > dipole_sphere_radius+moat_zone) then
                b_dipole_sqr=(bdipole_x(ixe,iye+1,ize+1)**2_CUSTOM_REAL+bdipole_y(ixe+1,iye+1,ize+1)**2_CUSTOM_REAL+bdipole_z(ixe+1,iye+1,ize+1)**2_CUSTOM_REAL)    &
                             *wpiwci**2
                profile_factor = max(denmin,min(1.0_CUSTOM_REAL,(2.0_CUSTOM_REAL/dsqrt(b_dipole_sqr))**profile_power))
              else
                profile_factor = 1._CUSTOM_REAL
              endif

            endif

          else

            profile_factor = 1._CUSTOM_REAL

          endif
          qp(np) = qp(np)*profile_factor
!
 
          vthfac = 1._CUSTOM_REAL/dsqrt(profile_factor)
          vcurr = 0.0_CUSTOM_REAL
 
 
          xran=myranf()
          !LAURA vmag=dsqrt(-alog(1._CUSTOM_REAL-.999999_CUSTOM_REAL*xran))
          vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*xran))  !LAURA
          xran=myranf()
          th=2._CUSTOM_REAL*pi*xran
          vxa=vpar(is)*vmag*dcos(th)*vthfac
 

          xran=myranf()
          !LAURA vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*xran))
          vmag=dsqrt(-log(1._CUSTOM_REAL-.999999_CUSTOM_REAL*xran)) !LAURA
          xran=myranf()
          th=2._CUSTOM_REAL*pi*xran
          vya=vper(is)*vmag*dcos(th)*vthfac
          vza=vper(is)*vmag*dsin(th)*vthfac
 
 
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
      ex=0._CUSTOM_REAL;ey=0._CUSTOM_REAL;ez=0._CUSTOM_REAL
 
 
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
                bz_avg= 0.125_CUSTOM_REAL*( bz       (i+1,j+1,k  ) &
                               +bz       (i  ,j+1,k  ) &
                               +bz       (i+1,j  ,k  ) &
                               +bz       (i  ,j  ,k  ) &
                               +bz       (i+1,j+1,k+1) &
                               +bz       (i  ,j+1,k+1) &
                               +bz       (i+1,j  ,k+1) &
                               +bz       (i  ,j  ,k+1) &
                              )                        &
                       +0.125_CUSTOM_REAL*( bdipole_z(i+1,j+1,k  ) &
                               +bdipole_z(i  ,j+1,k  ) &
                               +bdipole_z(i+1,j  ,k  ) &
                               +bdipole_z(i  ,j  ,k  ) &
                               +bdipole_z(i+1,j+1,k+1) &
                               +bdipole_z(i  ,j+1,k+1) &
                               +bdipole_z(i+1,j  ,k+1) &
                               +bdipole_z(i  ,j  ,k+1) &
                              )
                by_avg= 0.125_CUSTOM_REAL*( by       (i+1,j+1,k  ) &
                               +by       (i  ,j+1,k  ) &
                               +by       (i+1,j  ,k  ) &
                               +by       (i  ,j  ,k  ) &
                               +by       (i+1,j+1,k+1) &
                               +by       (i  ,j+1,k+1) &
                               +by       (i+1,j  ,k+1) &
                               +by       (i  ,j  ,k+1) &
                              )                        &
                       +0.125_CUSTOM_REAL*( bdipole_y(i+1,j+1,k  ) &
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
                 bz_avg= 0.250_CUSTOM_REAL*( bz       (i+1,j+1,k+1) &
                                +bz       (i  ,j+1,k+1) &
                                +bz       (i+1,j  ,k+1) &
                                +bz       (i  ,j  ,k+1) &
                               )                        &
                        +0.250_CUSTOM_REAL*( bdipole_z(i+1,j+1,k+1) &
                                +bdipole_z(i  ,j+1,k+1) &
                                +bdipole_z(i+1,j  ,k+1) &
                                +bdipole_z(i  ,j  ,k+1) &
                               )                        
                 by_avg= 0.250_CUSTOM_REAL*( by       (i+1,j+1,k+1) &
                                +by       (i  ,j+1,k+1) &
                                +by       (i+1,j  ,k+1) &
                                +by       (i  ,j  ,k+1) &
                               )                        &
                        +0.250_CUSTOM_REAL*( bdipole_y(i+1,j+1,k+1) &
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
        call xrealbcc(ey,one,nx,ny,nz)  !LAURA
        call xrealbcc(ez,one,nx,ny,nz)  !LAURA
      else
        call xrealbcc_pack_e_2d(ex,ey,ez,one,nx,ny,nz)  !LAURA
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
            fox(i,j,k)=0._CUSTOM_REAL
            foy(i,j,k)=0._CUSTOM_REAL
            foz(i,j,k)=0._CUSTOM_REAL
            eta(i,j,k)=resis
          enddo
        enddo
      enddo
      dtsav=dt
      dt=0._CUSTOM_REAL
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
      implicit none !LAURA
      real(kind=CUSTOM_REAL):: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi	!LAURA
      integer*8 ix,iy,iz,ixp1,iyp1,izp1
      real(kind=CUSTOM_REAL):: rfrac, wmult, h, hh, w1, w2, w3, w4, w5, &
             w6, w7, w8, dns1, dns2, dnst, vxavg1, vxavg2, vxavg, vyavg1, vyavg2, &
             vyavg, vzavg1,vzavg2,vzavg,vxa,vya,vza,bxa,bya,bza,btota,wpar, wperp2 !LAURA
      integer*8 :: is, iize, iiye, iixe, l

      call date_and_time(values=time_begin_array(:,23))

      dtxi = 1._CUSTOM_REAL/meshX%dt
      dtyi = 1._CUSTOM_REAL/meshY%dt
      dtzi = 1._CUSTOM_REAL/meshZ%dt


      tpar=0._CUSTOM_REAL
      tperp=0._CUSTOM_REAL
      if (nspec >= 2) then
         rfrac = frac(2)/frac(1)
      else
         rfrac = 0._CUSTOM_REAL
      endif
 
      call date_and_time(values=time_begin_array(:,26))
      DO IS=1,NSPEC
        wmult=wspec(is)
        h=dt*qspec(is)/wmult
        hh=.5_CUSTOM_REAL*h
        dpedx = 0._CUSTOM_REAL
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
!                ry=hyi*y(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                rz=hzi*z(l)+0.5000000000000001E+00_CUSTOM_REAL	!LAURA
!                ix=rx
!                iy=ry
!                iz=rz
!                fx=rx-ix
!                fy=ry-iy
!                fz=rz-iz

!               Nonuniform mesh - using MESH_UNMAP
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

                w1=(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)*(1._CUSTOM_REAL-fz)
                w2=    fx *(1._CUSTOM_REAL-fy)*(1._CUSTOM_REAL-fz)
                w3=(1._CUSTOM_REAL-fx)*    fy *(1._CUSTOM_REAL-fz)
                w4=    fx*     fy *(1._CUSTOM_REAL-fz)
                w5=(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)*    fz
                w6=    fx *(1._CUSTOM_REAL-fy)*    fz
                w7=(1._CUSTOM_REAL-fx)*    fy*     fz
                w8=    fx*     fy*     fz

                dns1= dns(ix  ,iy  ,iz  ,1)*w1+dns(ixp1,iy  ,iz  ,1)*w2  &
                +     dns(ix  ,iyp1,iz  ,1)*w3+dns(ixp1,iyp1,iz  ,1)*w4  &
                +     dns(ix  ,iy  ,izp1,1)*w5+dns(ixp1,iy  ,izp1,1)*w6  &
                +     dns(ix  ,iyp1,izp1,1)*w7+dns(ixp1,iyp1,izp1,1)*w8

                dns2= 0._CUSTOM_REAL

                dnst = dns1 + dns2
      
                vxavg1=vxs(ix  ,iy  ,iz  ,1)*w1+vxs(ixp1,iy  ,iz  ,1)*w2  &
                +      vxs(ix  ,iyp1,iz  ,1)*w3+vxs(ixp1,iyp1,iz  ,1)*w4  &
                +      vxs(ix  ,iy  ,izp1,1)*w5+vxs(ixp1,iy  ,izp1,1)*w6  &
                +      vxs(ix  ,iyp1,izp1,1)*w7+vxs(ixp1,iyp1,izp1,1)*w8

                vxavg2= 0._CUSTOM_REAL

                vxavg = (dns1*vxavg1 + dns2*vxavg2)/dnst

                vyavg1=vys(ix  ,iy  ,iz  ,1)*w1+vys(ixp1,iy  ,iz  ,1)*w2  &
                +      vys(ix  ,iyp1,iz  ,1)*w3+vys(ixp1,iyp1,iz  ,1)*w4  &
                +      vys(ix  ,iy  ,izp1,1)*w5+vys(ixp1,iy  ,izp1,1)*w6  &
                +      vys(ix  ,iyp1,izp1,1)*w7+vys(ixp1,iyp1,izp1,1)*w8

                vyavg2=0._CUSTOM_REAL
  
                vyavg = (dns1*vyavg1 + dns2*vyavg2)/dnst

                vzavg1=vzs(ix  ,iy  ,iz  ,1)*w1+vzs(ixp1,iy  ,iz  ,1)*w2  &
                +      vzs(ix  ,iyp1,iz  ,1)*w3+vzs(ixp1,iyp1,iz  ,1)*w4  &
                +      vzs(ix  ,iy  ,izp1,1)*w5+vzs(ixp1,iy  ,izp1,1)*w6  &
                +      vzs(ix  ,iyp1,izp1,1)*w7+vzs(ixp1,iyp1,izp1,1)*w8

                vzavg2=0._CUSTOM_REAL

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

                btota=dsqrt(bxa**2_CUSTOM_REAL+bya**2_CUSTOM_REAL+bza**2_CUSTOM_REAL)
                if (btota < 1.e-20_CUSTOM_REAL) btota=1.e-20_CUSTOM_REAL
                wpar=(vxa*bxa+vya*bya+vza*bza)/btota
                wperp2=vxa**2_CUSTOM_REAL+vya**2_CUSTOM_REAL+vza**2_CUSTOM_REAL-wpar**2_CUSTOM_REAL

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
 
 
                np=link(np)
              ENDDO
            ENDDO
          ENDDO
        ENDDO


!        call XREAL(tpar (1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL_4D(tpar (1,jb-1,kb-1,1),is,NX,NY,NZ,NSPEC)
!        call XREAL(tperp(1,jb-1,kb-1,is),NX,NY,NZ)
        call XREAL_4D(tperp(1,jb-1,kb-1,1),is,NX,NY,NZ,NSPEC)
        call XREAL(dpedx(1,jb-1,kb-1   ),NX,NY,NZ)

        DO IZ = KB-1,KE
          DO IY = JB-1,JE
            DO IX = 1, NX1
              if (dpedx(ix,iy,iz) /= 0.) then
                tpar (ix,iy,iz,is) = tpar (ix,iy,iz,is)/(   tx0(is)*dpedx(ix,iy,iz))
                tperp(ix,iy,iz,is) = tperp(ix,iy,iz,is)/(2._CUSTOM_REAL*tx0(is)*dpedx(ix,iy,iz))
              endif
            ENDDO
          ENDDO
        ENDDO


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
!                dns2=0._CUSTOM_REAL
!                denum=dns1+rfrac*dns2
!              else
!                denum=dns(i,j,k,is)/(dfac(is)*frac(is))
!              endif
!              if (denum < denmin)  then
!               tpar(i,j,k,is)=1.e-5_CUSTOM_REAL
!               tperp(i,j,k,is)=1.e-5_CUSTOM_REAL
!              else
!               denum=denum*tx0(is)
!               tpar(i,j,k,is)=tpar(i,j,k,is)*wspec(is)/denum
!               tperp(i,j,k,is)=0.5_CUSTOM_REAL*tperp(i,j,k,is)*wspec(is)/denum
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
use PRECISION    !LAURA
      implicit none
      integer*8 ibnd,i,j,nx1m,ny1m,nz1m,k
      real(kind=CUSTOM_REAL) a_x(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
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
   
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),CUSTOM_MPI_TYPE,nbrtop ,0,&	!LAURA
                        packed_data_xy_recv,size(packed_data_xy_recv),CUSTOM_MPI_TYPE,nbrbot ,0,&	!LAURA
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
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),CUSTOM_MPI_TYPE,nbrbot ,1,&	!LAURA
                        packed_data_xy_recv,size(packed_data_xy_recv),CUSTOM_MPI_TYPE,nbrtop ,1,&	!LAURA
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
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
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
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
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
      implicit none !LAURA
      integer*4 :: maxchar  !LAURA
 
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
use PRECISION    !LAURA
      implicit none
      integer*8 ibnd,i,j,nx1m,ny1m,nz1m,k
      real(kind=CUSTOM_REAL) a_x(nxmax,jb-1:je+1,kb-1:ke+1)&	!LAURA
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
   
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),CUSTOM_MPI_TYPE,nbrtop ,0,&	!LAURA
                        packed_data_xy_recv,size(packed_data_xy_recv),CUSTOM_MPI_TYPE,nbrbot ,0,&	!LAURA
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
      call MPI_SENDRECV(packed_data_xy_send,size(packed_data_xy_send),CUSTOM_MPI_TYPE,nbrbot ,1,&	!LAURA
                        packed_data_xy_recv,size(packed_data_xy_recv),CUSTOM_MPI_TYPE,nbrtop ,1,&	!LAURA
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
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
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
      call MPI_SENDRECV(packed_data_xz_send,size(packed_data_xz_send),CUSTOM_MPI_TYPE,nbrleft,0,&	!LAURA
                        packed_data_xz_recv,size(packed_data_xz_recv),CUSTOM_MPI_TYPE,nbrrite,0,&	!LAURA
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
      implicit none !LAURA

      real(kind=CUSTOM_REAL):: rx,ry,rz,fx,fy,fz,dtxi,dtyi,dtzi,w1,w2,w3,w4,w5,w6,w7,w8	!LAURA
      integer*8:: ix,iy,iz,ixp1,iyp1,izp1,i,j,k,jmin,jmax,kmin,kmax
      real(kind=CUSTOM_REAL),dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(in) :: nonuniform_mesh	!LAURA
      real(kind=CUSTOM_REAL),dimension(nxmax,jb-1:je+1,kb-1:ke+1), intent(out):: uniform_mesh	!LAURA
      real(kind=CUSTOM_REAL),dimension(nxmax,0:ny+1,0:nz+1), intent(out):: nonuniform_mesh_global	!LAURA
      real(kind=CUSTOM_REAL),dimension(nxmax,0:ny+1,0:nz+1):: nonuniform_mesh_local	!LAURA
      real(kind=CUSTOM_REAL):: xc_uniform_pos,yc_uniform_pos,zc_uniform_pos	!LAURA

      dtxi = 1._CUSTOM_REAL/meshX%dt
      dtyi = 1._CUSTOM_REAL/meshY%dt
      dtzi = 1._CUSTOM_REAL/meshZ%dt

      uniform_mesh          = 0._CUSTOM_REAL
      nonuniform_mesh_local = 0._CUSTOM_REAL

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
                         ,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,IERR)	!LAURA

      do i=2,nx+1
        xc_uniform_pos = (dble(i)-1.5_CUSTOM_REAL)*hx
        rx   = dtxi*MESH_UNMAP(meshX,xc_uniform_pos)+1.50000000000E+00_CUSTOM_REAL	!LAURA
        ix   = rx
        fx   = rx-ix
        ixp1 = ix+1
        do j=jb,je
          yc_uniform_pos = (dble(j)-0.5_CUSTOM_REAL)*hy
          ry   = dtyi*MESH_UNMAP(meshY,yc_uniform_pos)+1.50000000000E+00_CUSTOM_REAL	!LAURA
          iy   = ry
          fy   = ry-iy
          iy   = iy-1             ! integer index in y direction starts at 0
          iyp1 = iy+1
          do k=kb,ke
            zc_uniform_pos = (dble(k)-0.5_CUSTOM_REAL)*hz
            rz   = dtzi*MESH_UNMAP(meshZ,zc_uniform_pos)+1.50000000000E+00_CUSTOM_REAL	!LAURA
            iz   = rz
            fz   = rz-iz
            iz   = iz-1             ! integer index in z direction starts at 0
            izp1 = iz+1
       
            w1=(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)*(1._CUSTOM_REAL-fz)
            w2=fx     *(1._CUSTOM_REAL-fy)*(1._CUSTOM_REAL-fz)
            w3=(1._CUSTOM_REAL-fx)*fy     *(1._CUSTOM_REAL-fz)
            w4=fx     *fy     *(1._CUSTOM_REAL-fz)
            w5=(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)*fz
            w6=fx     *(1._CUSTOM_REAL-fy)*fz
            w7=(1._CUSTOM_REAL-fx)*fy     *fz
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
