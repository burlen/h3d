!
! DKV-1: parameter-dependent eta
!=======================================================================
!
      subroutine etacalc
!
!-----------------------------------------------------------------------
! - calculate a resistivity that depends on physical quantities
! - remember, eta and fo{x,y,z} are cell-centered quantities
! - ghost cells are only used in diagnostic output
! Here four/five choices, based on argument ieta (in input/block):
!       (1) constant eta
!       (2) 4th power of B-grad, 2nd power of density
!       (3) Arbitrary power of curent (here: 4)
!       (4) Homas method: 2nd derivative of current
!       (5) (2) and (4) combined, with (4) reduced by a factor (1/5)
!
!*************

      use parameter_mod
      use MESH2D
      real(kind=CUSTOM_REAL):: ajl(nxmax,jb-1:jb+nylmax,kb-1:kb+nzlmax) 	!LAURA
      real(kind=CUSTOM_REAL) :: dxa,dya,dza,dbxdy,dbxdz,dbydx,dbydz,dbzdx,dbzdy
      real(kind=CUSTOM_REAL) :: ba1, ba2, ba3, ba4, gb2, gb4, b2, cfront
      real(kind=CUSTOM_REAL) :: ajpx, ajmx, ajpy, ajmy, ajg, eps
      real(kind=CUSTOM_REAL) :: anetax, anetay, ainv4b, wpiwcigb, expo
      integer*8 :: k, j, i, ietg, ietn, ietb, ietj, l, itresis
      data eps / 1.e-25_CUSTOM_REAL/
 

        itresis=1000
 
 
      bx=bx+bdipole_x;by=by+bdipole_y;bz=bz+bdipole_z
 
 
      eta=0._CUSTOM_REAL
 
 
      if (ieta == 0)   then

        eta = resis
 
      else if (ieta == 1)   then

        anetax=dble(netax)
        anetay=dble(netay)
        do k=kb,ke
        do j=jb,je
          do i=1,nx2
            eta(i,j,k) = resis / ( dcosh((dble(i)-0.5_CUSTOM_REAL*dble(nx2+1))/anetax)         &
                                  *dcosh((dble(j)-0.5_CUSTOM_REAL*dble(ny2+1))/anetay))
          enddo
        enddo
        enddo
        eta=eta*dexp(-dble(it)/dble(itresis))	!LAURA
 
 
      else if ((ieta == 2).or.(ieta == 5) )   then
 
!                                use gradient of |B|, B, and n
!
!     good combination is (grad |B|)**4 / (n**4 * B**2)
!     set powers, here; and adjust cfront factor accordingly
!     e.g., cfront = 1e-3 for ietg=4, ietn=4, and ietb=2
!     this also uses an upper limit on eta, set desired value at top
 
        cfront = 8.0e-3_CUSTOM_REAL
        ietg = 4
        ietn = 4
        ietb = 2
        ainv4b = 1.0_CUSTOM_REAL/ ( dble(4_CUSTOM_REAL**ietb) )	!LAURA
        wpiwcigb = wpiwci**(dble(ietg) -dble(ietb))
        do k=kb,ke
          do j=jb,je
            do i=2,nx1
              if (den(i,j,k) .le. denmin)   then
                eta(i,j,k) = 0.0_CUSTOM_REAL
              else
                ba1 = ( bx(i+1,j  ,k)      &
                       +bx(i+1,j+1,k) )**2_CUSTOM_REAL &
                     +( by(i+1,j  ,k)      &
                       +by(i+1,j+1,k) )**2_CUSTOM_REAL &
                     +( bz(i+1,j  ,k)      &
                       +bz(i+1,j+1,k) )**2_CUSTOM_REAL
                ba2 = ( bx(i  ,j  ,k)      &
                       +bx(i  ,j+1,k) )**2_CUSTOM_REAL &
                     +( by(i  ,j  ,k)      &
                       +by(i  ,j+1,k) )**2_CUSTOM_REAL &
                     +( bz(i  ,j  ,k)      &
                       +bz(i  ,j+1,k) )**2_CUSTOM_REAL
                ba3 = ( bx(i  ,j+1,k)      &
                       +bx(i+1,j+1,k) )**2_CUSTOM_REAL &
                     +( by(i  ,j+1,k)      &
                       +by(i+1,j+1,k) )**2_CUSTOM_REAL &
                     +( bz(i  ,j+1,k)      &
                       +bz(i+1,j+1,k) )**2_CUSTOM_REAL
                ba4 = ( bx(i  ,j  ,k)      &
                       +bx(i+1,j  ,k) )**2_CUSTOM_REAL &
                     +( by(i  ,j  ,k)      &
                       +by(i+1,j  ,k) )**2_CUSTOM_REAL &
                     +( bz(i  ,j  ,k)      &
                       +bz(i+1,j  ,k) )**2_CUSTOM_REAL

!               Uniform mesh - Same as is in version 5.0
!                gb2 = 1.0 *(  ((sqrt(ba1) -sqrt(ba2))/ hx)**2    &
!                             +((sqrt(ba3) -sqrt(ba4))/ hy)**2  )

!               Nonuniform mesh
                gb2 = 1.0_CUSTOM_REAL *(  ((dsqrt(ba1) -dsqrt(ba2))/ meshX%dxc(i  ))**2_CUSTOM_REAL    &
                             +((dsqrt(ba3) - dsqrt(ba4))/ meshY%dxc(j+1))**2_CUSTOM_REAL  )


                gb4 = gb2**(dble(ietg)/2_CUSTOM_REAL)
                b2 = ainv4b   *(dsqrt(ba1) + dsqrt(ba2))**dble(ietb)
                eta(i,j,k) = cfront *wpiwcigb *resis *gb4/         &
                              (b2 *den(i,j,k)**dble(ietn))
                eta(i,j,k) = min(eta(i,j,k), etamax)
                if (eta(i,j,k) .lt. etamin)   eta(i,j,k) = 0.0_CUSTOM_REAL
              endif
            enddo
          enddo
        enddo
!
!-----------------------------------------------------------------------
!
      else if (ieta == 3)   then

!                                use arbitrary power of j
!
!     using the sum over all three components
!
!     uses arbitrary power ietj
!
!     adjust cfront factor accordingly
!     e.g., cfront = 1.0
!     this also uses an upper limit on eta, set desired value at top
 
        ietj = 4
        cfront = 2.0_CUSTOM_REAL *1.0e3_CUSTOM_REAL**dble(ietj)

!       Uniform mesh - Same as is in version 5.0
!        dxa=1.0/(2.*hx)
!        dya=1.0/(2.*hy)
        do k=kb,ke
          do j=jb,je

!           Nonuniform mesh
            dya=1.0_CUSTOM_REAL/(2.0_CUSTOM_REAL*meshY%dxc(j+1))

            do i=2,nx1

!             Nonuniform mesh
              dxa=1.0_CUSTOM_REAL/(2.0_CUSTOM_REAL*meshX%dxc(i  ))

              dbzdy= bz(i+1,j+1,k)+bz(i  ,j+1,k) &
                    -bz(i+1,j  ,k)-bz(i  ,j  ,k)
              dbzdx= bz(i+1,j+1,k)+bz(i+1,j  ,k) &
                    -bz(i  ,j+1,k)-bz(i  ,j  ,k)
              dbydx= by(i+1,j+1,k)+by(i+1,j  ,k) &
                    -by(i  ,j+1,k)-by(i  ,j  ,k)
              dbxdy= bx(i+1,j+1,k)+bx(i  ,j+1,k) &
                    -bx(i+1,j  ,k)-bx(i  ,j  ,k)
              ajl(i,j,k) = (dya*dbzdy          )**2_CUSTOM_REAL   &
                          +(         -dxa*dbzdx)**2_CUSTOM_REAL &
                          +(dxa*dbydx-dya*dbxdy)**2_CUSTOM_REAL
            enddo
          enddo
!
!                 --- end loop over squared components ---
!
!     use 1/2 power, because of square root of squared components
!     i.e., starting from |j| = sqrt(jx**2 +jy**2 +jz**2)
!
          expo = 0.5_CUSTOM_REAL *dble(ietj)	!LAURA
          do j=jb,je
            do i=2,nx1
              eta(i,j,k) = cfront *resis *ajl(i,j,k)**expo
              eta(i,j,k) = min(etamax, eta(i,j,k))
              if (eta(i,j,k) .lt. etamin)   eta(i,j,k) = 0.0_CUSTOM_REAL
            enddo
          enddo


        enddo
 
 
      else
 
        if ( (ieta .gt. 5).or.(ieta .lt. 0) )   then
          write(*,*) 'etacalc only accepts ieta 0, through 5, currently'
          stop
          return
        endif
 

      endif
 
                                                                                                       
!-----------------------------------------------------------------------
!
      if ( (ieta == 4) .or. (ieta == 5) )   then
!
!                                use 2nd gradient of j
!
!     this routine calculates something related to the
!     logarithmic derivative of grad j...
!     seems to be effective at inhibiting whistlers
!     protect against division by zero with "eps"
!
!     using the sum over all three components
!
!     uses arbitrary power ietj
!
!     adjust cfront factor accordingly
!     e.g., cfront = 1.0
!     this also uses an upper limit on eta, set desired value at top
!
      cfront = 0.5_CUSTOM_REAL
!
!     when adding this resitivity to the gradient method (case 2),
!     lower the contribution from this method:
!
        if (ieta == 5)  cfront = cfront/ 5.0_CUSTOM_REAL
          ietj = 4

!         Uniform mesh - Same as in version 5.0
!          dxa=1.0/(2.*hx)
!          dya=1.0/(2.*hy)
 
          do l = 1,3
            do k=kb,ke
              do j=jb,je

!               Nonuniform mesh
                dya=1.0_CUSTOM_REAL/(2.0_CUSTOM_REAL*meshY%dxc(j+1))

                do i=2,nx1

!                 Nonuniform mesh
                  dxa=1.0_CUSTOM_REAL/(2.0_CUSTOM_REAL*meshX%dxc(i  ))

                  if (l == 1)      then
                    dbzdy= bz(i+1,j+1,k)+bz(i  ,j+1,k) &
                          -bz(i+1,j  ,k)-bz(i  ,j  ,k)
                    ajl(i,j,k) =  dya*dbzdy
                    if (ieta == 4)  eta(i,j,k) = 0.0_CUSTOM_REAL
                  else if (l == 2)  then
                    dbzdx= bz(i+1,j+1,k)+bz(i+1,j  ,k) &
                          -bz(i  ,j+1,k)-bz(i  ,j  ,k)
                    ajl(i,j,k) = -dxa*dbzdx
                  else if (l == 3)  then
                    dbydx= by(i+1,j+1,k)+by(i+1,j  ,k) &
                          -by(i  ,j+1,k)-by(i  ,j  ,k)
                    dbxdy= bx(i+1,j+1,k)+bx(i  ,j+1,k) &
                          -bx(i+1,j  ,k)-bx(i  ,j  ,k)
                    ajl(i,j,k) =  dxa*dbydx -dya*dbxdy
                  endif
                enddo
              enddo
            enddo
!
!
!     set non-periodic boundary conditions
!
            do k=kb,ke
              do j=jb,je
                ajl(1  ,j,k) = ajl(2  ,j,k)
                ajl(nx2,j,k) = ajl(nx1,j,k)
              enddo
            enddo
 

            do k=kb,ke
              do i=1,nx2
                if (jb == 1) ajl(i,0   ,k) = ajl(i, 1,k)
                if (je == ny) ajl(i,ny1,k) = ajl(i,ny,k)
              enddo
            enddo
 
            do k=kb,ke
              do j=jb,je
                do i=2,nx1

!                 Uniform mesh - Same as is in version 5.0
                  ajpx = (ajl(i+1,j,k) -ajl(i,  j,k))/ hx
                  ajmx = (ajl(i,  j,k) -ajl(i-1,j,k))/ hx
                  ajpy = (ajl(i,j+1,k) -ajl(i,  j,k))/ hy
                  ajmy = (ajl(i,j  ,k) -ajl(i,j-1,k))/ hy
                  ajg = dabs( (ajpx - ajmx)/ hx/                  &
                       (dabs(ajpx) + dabs(ajmx) +eps) )**dble(ietj) &
                       *dabs( (ajpy - ajmy)/ hy/                  &
                       (dabs(ajpy) + dabs(ajmy) +eps) )**dble(ietj)

!                 Nonuniform mesh
!                  ajpx = (ajl(i+1,j,k) -ajl(i,  j,k))/ meshX%dxn(i+1)
!                  ajmx = (ajl(i,  j,k) -ajl(i-1,j,k))/ meshX%dxn(i+1)
!                  ajpy = (ajl(i,j+1,k) -ajl(i,  j,k))/ meshY%dxn(j+2)
!                  ajmy = (ajl(i,j  ,k) -ajl(i,j-1,k))/ meshY%dxn(j+2)
!                  ajg = abs( (ajpx - ajmx)/ meshX%dxc(i  )/                  &
!                       (abs(ajpx) +abs(ajmx) +eps) )**ietj                   &
!                       *abs( (ajpy - ajmy)/ meshY%dxc(j+1)/                  &
!                       (abs(ajpy) +abs(ajmy) +eps) )**ietj

                  eta(i,j,k) = eta(i,j,k) +cfront *resis *ajg
                enddo
              enddo
            enddo


          enddo


          do k=kb,ke
            do j=jb,je
              do i=2,nx1
                eta(i,j,k) = min(etamax, eta(i,j,k))
                if (eta(i,j,k) .lt. etamin)   eta(i,j,k) = 0.0_CUSTOM_REAL
              enddo
            enddo
          enddo
 

        endif
 
 
        bx=bx-bdipole_x;by=by-bdipole_y;bz=bz-bdipole_z
 
 
      return
      end
!
!

