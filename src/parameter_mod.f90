!
!     Declare global parameters and global arrays in this module
!
      module parameter_mod
      implicit double precision (a-h,o-z),integer*8(i-n)
      save
      include "mpif.h"
      integer:: my_short_int,i_source,i_destination,i_tag,i_length,i_i
      integer  ,dimension(8,128):: time_begin_array,time_end_array
      double precision,dimension(128):: time_elapsed
      integer*8 :: nxmax,nymax,nzmax,nspecm,npes,nvar,nylmax,nzlmax,npm,npes_over_60

      integer :: numprocs,ndim,dims(2),ierr,comm2d,myid,req(8)                                                      &
                ,nbrtop,nbrbot,nbrritetop,nbrlefttop,nbrritebot,nbrleftbot                                          &
                ,nbrleft,nbrrite,ipe,stridery,striderz,iseed(1),coords(2)
      integer:: status (mpi_status_size),status1(mpi_status_size),status2(mpi_status_size),status_array(mpi_status_size,8)

      double precision:: zb,ze,yb,ye,teti,dt_field,volume_fraction,cell_volume_ratio,zb_logical,ze_logical,yb_logical,ye_logical       &
                        ,xb_logical,xe_logical,xb,xe,dipole_ramp_time,dipole_moment_max,time_reverse_B,reverse_B_ramp_time,smooth_coef
      double precision, dimension(:), allocatable:: zbglobal,zeglobal,ybglobal,yeglobal,xc_uniform,yc_uniform,zc_uniform   &
                                                                                       ,xv_uniform,yv_uniform,zv_uniform
      integer*8, dimension(:), allocatable:: kbglobal,keglobal,jbglobal,jeglobal,nsendp,nrecvp,ixc_2_c_map,iyc_2_c_map,izc_2_c_map     &
                                                                                              ,ixc_2_v_map,iyc_2_v_map,izc_2_v_map     &
                                                                                              ,ixv_2_c_map,iyv_2_c_map,izv_2_c_map     &
                                                                                              ,ixv_2_v_map,iyv_2_v_map,izv_2_v_map     
      double precision, dimension(:,:,:), allocatable:: ex,ey,ez,bx,by,bz,fox,foy,foz,eta,curlex,curley,curlez,bxs, &
                                                        bys,bzs,den,deno,dpedx,dpedy,dpedz,vix,viy,viz,vixo,viyo,   &
                                                        vizo,pe,curlbx,curlby,curlbz,bdipole_x,bdipole_y,bdipole_z, &
                                                        dipole_sphere_ex,dipole_sphere_ey,dipole_sphere_ez,         &
                                                        eta_times_b_dot_j
      double precision, dimension(:,:,:,:), allocatable:: dns, vxs, vys, vzs, tpar, tperp,qp_cell
      double precision, dimension(:,:,:,:), allocatable:: p_xx,p_xy,p_xz,p_yy,p_yz,p_zz

      double precision, dimension(:,:), allocatable ::ainjxz,ainjzx,deavxz,deavzx,vxavxz,vyavxz,vzavxz,vxavzx,      &
                                                      vyavzx,vzavzx,vxcaxz,vycaxz,vzcaxz,vxcazx,vycazx,vzcazx,      &
                                                      ainjyz,ainjzy,deavyz,deavzy,vxavyz,vyavyz,vzavyz,vxavzy,      &
                                                      vyavzy,vzavzy,vxcayz,vycayz,vzcayz,vxcazy,vycazy,vzcazy,      &
                                                      ainjxy,ainjyx,deavxy,deavyx,vxavxy,vyavxy,vzavxy,vxavyx,      &
                                                      vyavyx,vzavyx,vxcaxy,vycaxy,vzcaxy,vxcayx,vycayx,vzcayx
      double precision, dimension(:), allocatable :: x, y, z, vx, vy, vz, qp
      integer*8, dimension(:), allocatable :: link,porder
      integer*8 :: nplmax,ipstore,np,n_subcycles,field_subcycles
      integer*8, dimension(:,:,:,:), allocatable:: iphead, iptemp
      integer*8, dimension(:), allocatable ::  ninj,ninj_global,nescape,nescape_global,npart,npart_global
      integer*8, dimension(:),allocatable:: nescape_yz,nescape_zy,nescape_xy                 &
                             ,nescape_yx,nescape_xz,nescape_zx                 &
                             ,nescape_yz_global,nescape_zy_global              &
                             ,nescape_xy_global,nescape_yx_global              &
                             ,nescape_xz_global,nescape_zx_global
      double precision, dimension(:), allocatable :: qleft,qrite
      double precision, dimension(:), allocatable:: x0,x1,tx0,vpar,vper,bbal
      double precision, dimension(:,:) ,allocatable:: vbal
      double precision, dimension(5):: rcorr
      integer*8, dimension(5):: ishape
      double precision, dimension(5):: btspec,qspec,wspec,frac,vxbar,vybar,vzbar,anisot
      double precision :: denmin, resis, wpiwci, bete, fxsho,ave1,ave2,phib,theta,demin2,xmax,ymax,zmax,dt,theta_max&
                         ,bxc,byc,bzc,gama,dtwci,dipole_moment,dipole_sphere_radius,R_MP,R_obstacle_to_MP,bz_IMF    &
                         ,QUOTA,wall_clock_elapsed,maximum_simulation_time,buffer_zone,exc,eyc,ezc                  &
                         ,xaa,xbb,yaa,ybb,zaa,zbb,t_stopped,theta_old
      integer*8::        nax,nbx,nay,nby,naz,nbz
      integer*8, dimension(8) :: wall_clock_begin,wall_clock_end
      integer*8:: i_dipole,j_dipole,k_dipole,eta_par,nparbuf
      integer*8, dimension(5) ::  npx,npy,npz
      integer*8 :: iterb,norbskip,restrt_write,nxcel,netax,netay,netaz,nspec,nx,ny,nz,nprint,nwrtdata &
                ,nskipx,nskipy,nskipz
      double precision:: etamin,etamax,moat_zone
      integer*8:: ieta,profile_power
      logical :: testorbt,restart,setup_mesh,uniform_loading_in_logical_grid,MPI_IO_format,image_dipole,smoothing
      double precision ::  hx,hy,hz,hxi,hyi,hzi,pi,efld,bfld,time,te0
      logical ::prntinfo, wrtdat
      integer*8 :: nsteps0,itfin,iwt,nx1,nx2,ny1,ny2,nz1,nz2,it,iopen,file_unit(20),file_unit_time,notime             &
                ,file_unit_tmp,file_unit_read(20),nptot,npleaving,npentering,iclock_speed,clock_time_init             &
                ,clock_time_old,clock_time,clock_time1,nptotp
      double precision, dimension(:) ,allocatable:: dfac
      integer*8, dimension(:) ,allocatable:: nskip,ipleft,iprite,ipsendleft,ipsendrite,iprecv,ipsendtop,ipsendbot     &
                                          ,ipsendlefttop,ipsendleftbot,ipsendritetop,ipsendritebot,ipsend
      logical global
      integer*8:: idum
      logical:: harris
      integer*8, dimension(:)  , allocatable:: idmap
      integer*8, dimension(:,:), allocatable:: idmap_yz
      integer*8:: kb,ke,jb,je,nsendtotp,nrecvtotp,nbrid(8),nsendtot,nrecvtot
      integer*8, dimension(:), allocatable:: idfft,kvec,jvec,myid_stop
      integer*8:: ihstb,ihste,isendid(4),irecvid(4,4)
      real*4:: single_prec
      double precision:: double_prec, xtmp1m, xtmp2m,initial_time,final_time    &
                        ,xbox_l,xbox_r,ybox_l,ybox_r,zbox_l,zbox_r,t_begin,t_end
      integer*8::   recl_for_real,recl_for_double_precision
      logical:: periods(2),reorder,Yee,b_in_xy,absorbing_dipole,post_process
      double precision:: u_array_injection(1001,-1000:1000)
      integer*8:: restart_index
      character(len=2):: restart_index_suffix(2)
      character(len=160):: data_directory,cycle_ascii_new,myid_char,restart_directory
      character(len=160):: cycle_ascii,cleanup_status
!----------------------------------------------------------------------
!
! insitu Patrick O'Leary 2/21/2013
!
!  Declaration
!----------------------------------------------------------------------
      logical :: insitu
!----------------------------------------------------------------------
! end Declaration
!----------------------------------------------------------------------
!
!

!
!
      contains
!
!     Set global parameters and allocate global arrays
!
        subroutine set_parameters(numprocs)
          integer:: numprocs
          double_prec=0.
          single_prec=0.
          inquire (IOLENGTH=recl_for_double_precision) double_prec
          inquire (IOLENGTH=recl_for_real) single_prec
          nparbuf=nxmax*(nylmax+2)*(nzlmax+2)
          npes=numprocs
          npes_over_60 = npes / 60
!         count particle storage requirement
          nplmax=0
          do is=1,nspec
            nplmax=nplmax+npx(is)*npy(is)*npz(is)
          enddo
          nplmax=5*nplmax       ! pad storage requirement by factor of 2 
 
          allocate (zbglobal(0:npes-1),zeglobal(0:npes-1),ybglobal(0:npes-1),yeglobal(0:npes-1)                     &
                   ,kbglobal(0:npes-1),keglobal(0:npes-1),jbglobal(0:npes-1),jeglobal(0:npes-1)                     &
                   ,nsendp(0:npes-1),nrecvp(0:npes-1),myid_stop(0:npes-1))
          allocate ( x(nplmax),y(nplmax),z(nplmax),vx(nplmax),vy(nplmax),vz(nplmax),link(nplmax),porder(nplmax)     &
                    ,qp(nplmax))
          allocate ( ninj(nspecm), ninj_global(nspecm),nescape(nspecm),nescape_global(nspecm),npart(nspecm)         &
                    ,npart_global(nspecm),qleft(nspecm),qrite(nspecm))
          allocate ( nescape_xy(nspecm),nescape_yx(nspecm),nescape_xz(nspecm),nescape_zx(nspecm)                    &
                    ,nescape_yz(nspecm),nescape_zy(nspecm)                                                          &
                    ,nescape_xy_global(nspecm),nescape_yx_global(nspecm),nescape_xz_global(nspecm)                  &
                    ,nescape_zx_global(nspecm),nescape_yz_global(nspecm),nescape_zy_global(nspecm))
          allocate ( x0(nspecm),x1(nspecm),tx0(nspecm),vpar(nspecm),vper(nspecm),vbal(nxmax,nspecm),bbal(nxmax))
          allocate ( dfac(nspecm),nskip(nspecm),ipleft(nspecm),iprite(nspecm),ipsendleft(nspecm),ipsendrite(nspecm) &
                    ,iprecv(nspecm),ipsendtop(nspecm),ipsendbot(nspecm),ipsendlefttop(nspecm),ipsendleftbot(nspecm) &
                    ,ipsendritetop(nspecm),ipsendritebot(nspecm),ipsend(nspecm))
         allocate ( idmap_yz(0:nymax+1,0:nzmax+1),idmap(0:nzmax),idfft(nzmax),kvec(nzlmax),jvec(nylmax))
        end subroutine set_parameters
!
!       Allocate global arrays
!
        subroutine allocate_global_arrays
          allocate ( ex       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),ey       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,ez       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bx       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,by       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bz       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,fox      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),foy      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,foz      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),eta      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,curlex   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),curley   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,curlez   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bxs      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,bys      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bzs      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,den      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),deno     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,dpedx    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),dpedy    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,dpedz    (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),vix      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,viy      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),viz      (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,vixo     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),viyo     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,vizo     (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),pe       (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,curlbx   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),curlby   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,curlbz   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bdipole_x(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,bdipole_y(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax),bdipole_z(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)   &
                    ,dipole_sphere_ex(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)                                           &
                    ,dipole_sphere_ey(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)                                           &
                    ,dipole_sphere_ez(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax)                                           &
                    ,eta_times_b_dot_j(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax))
          allocate ( dns(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),vxs   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
                    ,vys(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),vzs   (nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
                    ,tpar(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),tperp(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) &
                    ,qp_cell(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm)) 
          allocate ( p_xx(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),p_xy(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm)&
                    ,p_xz(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),p_yy(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm)&
                    ,p_yz(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),p_zz(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm) )
          allocate ( ainjxz(nxmax,kb-1:kb+nzlmax),ainjzx(nxmax,kb-1:kb+nzlmax),deavxz(nxmax,kb-1:kb+nzlmax)          &
                    ,deavzx(nxmax,kb-1:kb+nzlmax),vxavxz(nxmax,kb-1:kb+nzlmax),vyavxz(nxmax,kb-1:kb+nzlmax)          &
                    ,vzavxz(nxmax,kb-1:kb+nzlmax),vxavzx(nxmax,kb-1:kb+nzlmax),vyavzx(nxmax,kb-1:kb+nzlmax)          &
                    ,vzavzx(nxmax,kb-1:kb+nzlmax),vxcaxz(nxmax,kb-1:kb+nzlmax),vycaxz(nxmax,kb-1:kb+nzlmax)          &
                    ,vzcaxz(nxmax,kb-1:kb+nzlmax),vxcazx(nxmax,kb-1:kb+nzlmax),vycazx(nxmax,kb-1:kb+nzlmax)          &
                    ,vzcazx(nxmax,kb-1:kb+nzlmax))
          allocate ( ainjyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),ainjzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,deavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),deavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,vxavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vyavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,vzavyz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vxavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,vyavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax),vzavzy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,vxcayz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vycayz(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,vzcayz(jb-1:jb+nylmax,kb-1:kb+nzlmax),vxcazy(jb-1:jb+nylmax,kb-1:kb+nzlmax)                     &
                    ,vycazy(jb-1:jb+nylmax,kb-1:kb+nzlmax),vzcazy(jb-1:jb+nylmax,kb-1:kb+nzlmax))
          allocate ( ainjxy(nxmax,jb-1:jb+nylmax),ainjyx(nxmax,jb-1:jb+nylmax),deavxy(nxmax,jb-1:jb+nylmax)          &
                    ,deavyx(nxmax,jb-1:jb+nylmax),vxavxy(nxmax,jb-1:jb+nylmax),vyavxy(nxmax,jb-1:jb+nylmax)          &
                    ,vzavxy(nxmax,jb-1:jb+nylmax),vxavyx(nxmax,jb-1:jb+nylmax),vyavyx(nxmax,jb-1:jb+nylmax)          &
                    ,vzavyx(nxmax,jb-1:jb+nylmax),vxcaxy(nxmax,jb-1:jb+nylmax),vycaxy(nxmax,jb-1:jb+nylmax)          &
                    ,vzcaxy(nxmax,jb-1:jb+nylmax),vxcayx(nxmax,jb-1:jb+nylmax),vycayx(nxmax,jb-1:jb+nylmax)          &
                    ,vzcayx(nxmax,jb-1:jb+nylmax))
          allocate (iphead(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm),iptemp(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax,nspecm))
          allocate (xc_uniform(nxmax),yc_uniform(nymax),zc_uniform(nzmax),xv_uniform(nxmax),yv_uniform(nymax),zv_uniform(nzmax))
          allocate (ixc_2_c_map(nx+1),iyc_2_c_map(ny+1),izc_2_c_map(nz+1))
          allocate (ixc_2_v_map(nx+1),iyc_2_v_map(ny+1),izc_2_v_map(nz+1))
          allocate (ixv_2_c_map(nx+1),iyv_2_c_map(ny+1),izv_2_c_map(nz+1))
          allocate (ixv_2_v_map(nx+1),iyv_2_v_map(ny+1),izv_2_v_map(nz+1))
        end subroutine allocate_global_arrays
!
!
      end module parameter_mod

