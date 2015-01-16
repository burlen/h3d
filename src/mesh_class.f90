!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! by Yuri Omelchenko, Scibernet Inc., 2003
! This module constructs a 1-D mesh by mapping a logical
! node-centered grid, t=[0,1] to a physical node-centered grid x=[0,xl] 
! surrounded by ghost grid points. 
! The cell-centered grid points are also computed.
! The following mapping is used (dt=1/nl):
! 0 <= t <= ta(=dt*na):
! x = xa - xa*[exp(alph1*(ta-t))-1]/[exp(alph1*ta)-1]
! ta <= t <= tb(=dt*nb):
! x = xa + dx*(t-ta)/dt, where dx = (xb-xa)/(nb-na)
! tb <= t <= 1:
! x = xb + (xl-xb)*[exp(alph2*(t-tb))-1]/[exp(alph2*(1-tb))-1]
! We define eps == exp(alph*dt)
! The first steps on both nonuniform grids are matched to be equal to dx.
! We solve nonlinear equations to find eps1 and eps2:
! eps1: (eps1-1)/(eps1**na-1) = dx/xa 
! eps2: (eps2-1)/(eps2**(nl-nb)-1) = dx/(xl-xb)
! NOTE: in order eps1>0 and eps2>0 we require:
! xa/na > dx and (xl-xb)/(nl-nb) > dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      module MESH_CLASS
use PRECISION    !LAURA
      implicit none

      type MESH
        real(kind=CUSTOM_REAL), dimension(:), pointer :: xn,xc,dxc,dxn	!LAURA
        integer*8 na,nb,nl
        real(kind=CUSTOM_REAL) xa,xb,xl	!LAURA
        real(kind=CUSTOM_REAL) dt,dx,dtdx,ta,tb,epsa,epsb,ca1,ca2,cb1,cb2	!LAURA
      end type MESH

      type MESHTYPE
        integer*8 type
      end type MESHTYPE

      type (MESHTYPE), parameter :: CELL = MESHTYPE(0)
      type (MESHTYPE), parameter :: NODE = MESHTYPE(1)
 
      interface MESH_INDEX
         module procedure MESH_INDEX_YURI,MESH_INDEX_HXV
      end interface

      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     initialize mesh attributes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      subroutine MESH_INIT(m,xa,xb,xl,na,nb,nl)
use PRECISION    !LAURA
      implicit none
      type(MESH), intent(out) :: m
      real(kind=CUSTOM_REAL), intent(in) :: xa,xb,xl	!LAURA
      integer*8, intent(in) :: na,nb,nl
      integer*8 i,nbb
      real(kind=CUSTOM_REAL) FINDEXP	!LAURA
     
      if((xa.ge.xb).or.(na.ge.nb)) then 
        call ERROR_ABORT('MESH_INIT(): bad parameters --- stop!')
      endif

      m%na=na ; m%nb=nb ; m%nl=nl
      m%xa=xa ; m%xb=xb ; m%xl=xl
      nbb = nl-nb

      allocate(m%xn(nl+3))  ! -1:nl+1
      allocate(m%xc(nl+2))  ! -1:nl
      allocate(m%dxn(nl+3)) ! -1:nl+1
      allocate(m%dxc(nl+2)) ! -1:nl

      m%dt = 1._CUSTOM_REAL/dble(nl)
      m%dx = (xb-xa)/(nb-na) 
      m%dtdx = m%dt/m%dx
      m%ta = m%dt*na
      m%tb = m%dt*nb
      if(na.gt.0)  then
         m%epsa = FINDEXP(m%dx/xa,na)
         m%ca1 = (m%epsa**na-1._CUSTOM_REAL)/xa
         m%ca2 = m%dt/dlog(m%epsa)
      else
         m%epsa = 1._CUSTOM_REAL
         m%ca1 = 0._CUSTOM_REAL
         m%ca2 = 0._CUSTOM_REAL
      endif
      if(nbb.gt.0) then
         m%epsb = FINDEXP(m%dx/(xl-xb),nbb)
         m%cb1 = (m%epsb**nbb-1._CUSTOM_REAL)/(xl-xb)
         m%cb2 = m%dt/dlog(m%epsb) 
      else
         m%epsb = 1._CUSTOM_REAL
         m%cb1 = 0._CUSTOM_REAL
         m%cb2 = 0._CUSTOM_REAL
      endif

      do i = 0,na-1
         m%xn(i+2) = xa - (m%epsa**(dble(na-i))-1._CUSTOM_REAL)/m%ca1 
         m%xc(i+2) = xa - (m%epsa**(dble(na-i)-0.5_CUSTOM_REAL)-1._CUSTOM_REAL)/m%ca1 
      enddo
      do i = na,nb
         m%xn(i+2) = xa + m%dx*(dble(i-na)) 
         m%xc(i+2) = xa + m%dx*(i+0.5_CUSTOM_REAL-na) 
      enddo
      do i = nb+1,nl
         m%xn(i+2) = xb + (m%epsb**(dble(i-nb))-1._CUSTOM_REAL)/m%cb1
         m%xc(i+2) = xb + (m%epsb**(dble(i+0.5_CUSTOM_REAL-nb))-1._CUSTOM_REAL)/m%cb1
      enddo

      m%xn(2) = 0. ! correct round-off errors

      m%xn(1)    = 2.*m%xn(2)-m%xn(3)
      m%xc(1)    = 2.*m%xn(2)-m%xc(2)
      m%xn(nl+3) = 2.*m%xn(nl+2)-m%xn(nl+1)
      m%xc(nl+2) = 2.*m%xn(nl+2)-m%xc(nl+1)
!     m%xc(1)    = 0.5*(m%xn(1)+m%xn(2))
!     m%xc(nl+2) = 0.5*(m%xn(nl+2)+m%xn(nl+3))

!     compute cell-based mesh sizes
      do i = 1,nl+2
         m%dxc(i) = m%xn(i+1)-m%xn(i)
      enddo
!     compute node-based mesh sizes
      do i = 2,nl+2
         m%dxn(i) = m%xc(i)-m%xc(i-1)
      enddo
      m%dxn(1) = m%dxn(2) ! fictitous
      m%dxn(nl+3) = m%dxn(nl+2) ! fictitous
      
      return
      end subroutine MESH_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     destroy mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      subroutine MESH_DESTRUCT(m)
use PRECISION    !LAURA
      implicit none
      type(MESH), intent(inout) :: m
      if(associated(m%xn)) then
        deallocate(m%xn)  
        deallocate(m%xc)  
        deallocate(m%dxc)  
        deallocate(m%dxn)  
        nullify(m%xn)
        nullify(m%xc)
        nullify(m%dxc)
        nullify(m%dxn)
      endif
      m%na = 0 ; m%nb = 0 ; m%nl = 0
      end subroutine MESH_DESTRUCT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     transform physical coordinate to logical space 
!     MESH_INIT() must be called prior to this call
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      function MESH_UNMAP(m,x) 
use PRECISION    !LAURA
      implicit none
      type(MESH), intent(in) :: m
      real(kind=CUSTOM_REAL), intent(in) :: x 	!LAURA
      real(kind=CUSTOM_REAL) MESH_UNMAP, t	!LAURA
      if(x.lt.m%xa) then
        t = m%ta - m%ca2*dlog( m%ca1*(m%xa-x)+1._CUSTOM_REAL )
      else if(x.le.m%xb) then
        t = m%ta + (x-m%xa)*m%dtdx 
      else 
        t = m%tb + m%cb2*dlog( m%cb1*(x-m%xb)+1._CUSTOM_REAL ) 
      endif
!     prevent mapping outiside of logical interval [0,1]
      if(t >= 1.) t = 1._CUSTOM_REAL-epsilon(dble(1)) 
      if(t <= 0.) t = epsilon(dble(1)) 
      MESH_UNMAP = t
      return
      end function MESH_UNMAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     transform physical coordinate to logical space 
!     MESH_INIT() must be called prior to this call
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      function MESH_MAP(m,t)
use PRECISION    !LAURA
      implicit none
      type(MESH), intent(in) :: m
      real(kind=CUSTOM_REAL), intent(in) :: t 	!LAURA
      real(kind=CUSTOM_REAL) MESH_MAP,x	!LAURA
      if(t.lt.m%ta) then
        x = m%xa-(dexp((m%ta-t)/m%ca2)-1._CUSTOM_REAL)/m%ca1
      else if(t.le.m%tb) then
        x = m%xa + (t-m%ta)/m%dtdx
      else 
        x = m%xb + (dexp((t-m%tb)/m%cb2)-1._CUSTOM_REAL)/m%cb1
      endif
!     prevent mapping outside of physical interval [0,xl]
      if(x >= m%xl) x = m%xl-epsilon(dble(1)) 
      if(x <= 0.) x = epsilon(dble(1)) 
      MESH_MAP = x
      return
      end function MESH_MAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     let xu(i=1:nnx)=[0,xl] be a uniform node-centered grid
!     then for each xu(i) find ix(i) such as 
!     m%xc(ix) <= xu(i) < m%xc(ix+1) if inode==CELL
!     m%xn(ix) <= xu(i) < m%xn(ix+1) if inode==NODE
!     MESH_INIT() must be called prior to this call
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      subroutine MESH_INDEX_HXV(m,inode,ix,inode_uniform)
use PRECISION    !LAURA
      implicit none
      type(MESH), intent(in) :: m
      type(MESHTYPE), intent(in) :: inode,inode_uniform
      integer*8, intent(out), dimension(:) :: ix
      real(kind=CUSTOM_REAL), dimension(:), pointer :: pp	!LAURA
      real(kind=CUSTOM_REAL) x,hx 	!LAURA
      integer*8 i,k,nnx

      ix(:)=0  
      nullify(pp)

      if(inode%type == CELL%type) then ! cell-centered
         pp => m%xc 
      else
         pp => m%xn 
      endif

      nnx = size(ix)
      hx = m%xl/(nnx-1) ! nnx == number of nodes
      do i=1,nnx
        if (inode_uniform%type == CELL%type) then     ! Interpolate locations of CELL uniform mesh
          x=(dble(i-1)-0.5_CUSTOM_REAL)*hx
        else                                      ! Interpolate locations of NODE uniform mesh
          x=(dble(i-1)    )*hx
        endif
        do k=2,size(pp)
           if(x.lt.pp(k)) then
             ix(i)=k-1
             exit
           endif
        enddo
      enddo

      nullify(pp)
      end subroutine MESH_INDEX_HXV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     let xu(i=1:nnx)=[0,xl] be a uniform node-centered grid
!     then for each xu(i) find ix(i) such as 
!     m%xc(ix) <= xu(i) < m%xc(ix+1) if inode==CELL
!     m%xn(ix) <= xu(i) < m%xn(ix+1) if inode==NODE
!     MESH_INIT() must be called prior to this call
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      subroutine MESH_INDEX_YURI(m,inode,ix)
use PRECISION    !LAURA
      implicit none
      type(MESH), intent(in) :: m
      type(MESHTYPE), intent(in) :: inode
      integer*8, intent(out), dimension(:) :: ix
      real(kind=CUSTOM_REAL), dimension(:), pointer :: pp	!LAURA
      real(kind=CUSTOM_REAL) x,hx 	!LAURA
      integer*8 i,k,nnx

      ix(:)=0  
      nullify(pp)

      if(inode%type == CELL%type) then ! cell-centered
         pp => m%xc 
      else
         pp => m%xn 
      endif

      nnx = size(ix)
      hx = m%xl/(nnx-1) ! nnx == number of nodes
      do i=1,nnx
        x=(i-1    )*hx
        do k=2,size(pp)
           if(x.lt.pp(k)) then
             ix(i)=k-1
             exit
           endif
        enddo
      enddo

      nullify(pp)
      end subroutine MESH_INDEX_YURI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!  interpolate (1:nx2,1:ny2) arrays from a nonuniform cell-centered grid 
!  to a uniform node-centered grid (1:nnx,1:nny) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      subroutine MESH_INTERPOLATE2D(inode,mx,my,a,ncx,ncy,ap,nnx,nny,lfirst)
use PRECISION    !LAURA
      implicit none
      type(MESHTYPE), intent(in) :: inode
      type(MESH), intent(in) :: mx,my 
      integer*8, intent(in) :: ncx,ncy,nnx,nny
      real(kind=CUSTOM_REAL), dimension(ncx,ncy), intent(in) :: a	!LAURA
      real(kind=CUSTOM_REAL), dimension(nnx,nny), intent(out) :: ap	!LAURA
      logical, intent(in) :: lfirst
!     integer*8, dimension(nnx), save :: ix  ! won't work under Windows
!     integer*8, dimension(nny), save :: iy  ! won't work under Windows
      integer*8, allocatable, save :: ix(:)
      integer*8, allocatable, save :: iy(:)

      real(kind=CUSTOM_REAL), save :: dx,dy	!LAURA
      real(kind=CUSTOM_REAL), dimension(:), pointer :: ppx,ppy	!LAURA
      integer*8 i,j,ii,jj
      real(kind=CUSTOM_REAL) rx,ry,fx,fy,w1,w2,w3,w4	!LAURA

      ap(:,:) = 0 ! zero output array 
      nullify(ppx)
      nullify(ppy)

      if(inode%type == CELL%type) then
        ppx => mx%xc
        ppy => my%xc
      else
        ppx => mx%xn
        ppy => my%xn
      endif

      if(lfirst) then
         if(allocated(ix)) then
           deallocate(ix)
           deallocate(iy)
         endif 
         allocate(ix(nnx))
         allocate(iy(nny))
         call MESH_INDEX(mx,inode,ix)
         call MESH_INDEX(my,inode,iy)
         dx = mx%xl/(nnx-1)
         dy = my%xl/(nny-1)
      endif

      do j=1,nny

         ry=(dble(j-1))*dy
         jj=iy(j)
         fy=( ry-ppy(jj) )/( ppy(jj+1)-ppy(jj) )

         do i=1,nnx

            rx=(i-1)*dx
            ii=ix(i)
            fx=( rx-ppx(ii) )/( ppx(ii+1)-ppx(ii) )
            w1=(1._CUSTOM_REAL-fx)*(1._CUSTOM_REAL-fy)
            w2=fx*(1._CUSTOM_REAL-fy)
            w3=(1._CUSTOM_REAL-fx)*fy
            w4=fx*fy
            ap(i,j)=w1*a(ii,jj)+w2*a(ii+1,jj)+w3*a(ii,jj+1)+w4*a(ii+1,jj+1)
        
         enddo

      enddo

      nullify(ppx)
      nullify(ppy)
      end subroutine MESH_INTERPOLATE2D

      end module MESH_CLASS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     Helper functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Solve for x: (x-1)/(x^N-1) - rhs = 0 
!     by using a simple bisection method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function FINDEXP(rhsi,ni)
use PRECISION    !LAURA
      implicit none
      real(kind=CUSTOM_REAL), intent(in) :: rhsi	!LAURA
      integer*8, intent(in) :: ni
      real(kind=CUSTOM_REAL) FINDEXP,tol,rhs,af,bf	!LAURA
      integer*8 n
      common /fparams/ rhs,n
      real(kind=CUSTOM_REAL) FUNC	!LAURA
      external FUNC
      tol = 10.*epsilon(dble(0))
!     These are common block parameters
      rhs=rhsi
      n = ni
!     These are common block parameters
      if(FUNC(1._CUSTOM_REAL).le.tol) then
        call ERROR_ABORT('FINDEXP(): dx_uniform too large --- stop!')
      endif
      af = 1._CUSTOM_REAL
      bf = af
      do while(FUNC(bf).ge.0.)
         bf = bf*2._CUSTOM_REAL
      enddo
      call BISECT(FUNC,af,bf,tol)
      FINDEXP = 0.5_CUSTOM_REAL*(af+bf)
      return
      end function FINDEXP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     F(t) = 0 to solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function FUNC(t) 
use PRECISION    !LAURA
      implicit none
      real(kind=CUSTOM_REAL), intent(in) :: t 	!LAURA
      real(kind=CUSTOM_REAL) FUNC,rhs	!LAURA
      integer*8 n
      common /fparams/ rhs,n
      if(t.le.1.) then
        FUNC = 1._CUSTOM_REAL/dble(n)-rhs
      else
        FUNC = (t-1._CUSTOM_REAL)/(t**n-1._CUSTOM_REAL)-rhs
      endif
      return
      end function FUNC 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     simple root finder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine BISECT(f,a,b,tol)
use PRECISION    !LAURA
      implicit none
      real(kind=CUSTOM_REAL), intent(in) :: tol	!LAURA
      real(kind=CUSTOM_REAL), intent(inout) :: a,b	!LAURA
      real(kind=CUSTOM_REAL) c,u,v,w	!LAURA
      real(kind=CUSTOM_REAL) f	!LAURA
      u = f(a) ; v = f(b)      
      if(u*v.gt.0.) then
        call ERROR_ABORT('BISECT(): bad initial interval --- stop!')
      else if(u.eq.0.) then
        b=a
        return
      else if(v.eq.0.) then
        a=b
        return
      endif
      do while((abs(f(a)).gt.tol).or.(abs(f(b)).gt.tol) ) 
         c = 0.5_CUSTOM_REAL*(a+b) 
         w = f(c)    
         if(w.eq.0.) then 
           a=c ; b=c
           return;
         endif
         if(w*u.lt.0.) then
           b=c ; v=w       
         else
           a=c ; u=w       
         endif
      enddo
      return      
      end subroutine BISECT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ERROR_ABORT(message)
      character(*), intent(in) :: message
      write(6,*) message
      stop
      end subroutine ERROR_ABORT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WARNING(message)
      character(*), intent(in) :: message
      write(6,*) message
      end subroutine WARNING
