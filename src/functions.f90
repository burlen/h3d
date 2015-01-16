!
!***********************************************************************
!
      module gammaln_mod
      contains
        FUNCTION gammln(xx)
use PRECISION    !LAURA
        implicit none
        real(kind=CUSTOM_REAL) gammln,xx	!LAURA
        INTEGER*8 j
        real(kind=CUSTOM_REAL) ser,stp,tmp,x,y,cof(6)	!LAURA
        SAVE cof,stp
        DATA cof,stp/76.18009172947146E+0_CUSTOM_REAL,-86.50532032941677E+0_CUSTOM_REAL,&	!LAURA
        24.01409824083091E+0_CUSTOM_REAL,-1.231739572450155E+0_CUSTOM_REAL,.1208650973866179E-2_CUSTOM_REAL,&	!LAURA
        -.5395239384953E-5_CUSTOM_REAL,2.5066282746310005E+0_CUSTOM_REAL/	!LAURA
        x=xx
        y=x
        tmp=x+5.5E+0_CUSTOM_REAL	!LAURA
        tmp=(x+0.5E+0_CUSTOM_REAL)*dlog(tmp)-tmp	!LAURA
        ser=1.000000000190015E+0_CUSTOM_REAL	!LAURA
        do 11 j=1,6
          y=y+1.E+0_CUSTOM_REAL	!LAURA
          ser=ser+cof(j)/y
11      continue
        gammln=tmp+dlog(stp*ser/x)
        END FUNCTION gammln
      end module gammaln_mod
!
!***********************************************************************
!
      module gcf_mod
      contains
        SUBROUTINE gcf(gammcf,a,x,gln)
        use gammaln_mod
use PRECISION    !LAURA
        implicit none
        INTEGER*8 ITMAX
        real(kind=CUSTOM_REAL) a,gammcf,gln,x,EPS,FPMIN	!LAURA
        PARAMETER (ITMAX=100,EPS=3.e-7_CUSTOM_REAL,FPMIN=1.e-30_CUSTOM_REAL)
!U      USES gammln
        INTEGER*8 i
        real(kind=CUSTOM_REAL) an,b,c,d,del,h	!LAURA
        gln=gammln(a)
        b=x+1._CUSTOM_REAL-a
        c=1._CUSTOM_REAL/FPMIN
        d=1._CUSTOM_REAL/b
        h=d
        do 11 i=1,ITMAX
          an=-dble(i)*(dble(i)-a)
          b=b+2._CUSTOM_REAL
          d=an*d+b
          if(dabs(d).lt.FPMIN)d=FPMIN
          c=b+an/c
          if(dabs(c).lt.FPMIN)c=FPMIN
          d=1._CUSTOM_REAL/d
          del=d*c
          h=h*del
          if(dabs(del-1._CUSTOM_REAL).lt.EPS)goto 1
11      continue
        pause 'a too large, ITMAX too small in gcf'
1       gammcf=dexp(-x+a*dlog(x)-gln)*h
        END SUBROUTINE gcf
      end module gcf_mod
!
!***********************************************************************
!
      module gser_mod
      contains
        SUBROUTINE gser(gamser,a,x,gln)
        use gammaln_mod
use PRECISION    !LAURA
        implicit none
        INTEGER*8 ITMAX
        real(kind=CUSTOM_REAL) a,gamser,gln,x,EPS	!LAURA
        PARAMETER (ITMAX=100,EPS=3.e-7_CUSTOM_REAL)
!U      USES gammln
        INTEGER*8 n
        real(kind=CUSTOM_REAL) ap,del,sum	!LAURA
        gln=gammln(a)
        if(x.le.0.)then
          if(x.lt.0.)pause 'x < 0 in gser'
          gamser=0._CUSTOM_REAL
          return
        endif
        ap=a
        sum=1._CUSTOM_REAL/a
        del=sum
        do 11 n=1,ITMAX
          ap=ap+1._CUSTOM_REAL
          del=del*x/ap
          sum=sum+del
          if(dabs(del).lt.dabs(sum)*EPS)goto 1
11      continue
        pause 'a too large, ITMAX too small in gser'
1       gamser=sum*dexp(-x+a*dlog(x)-gln)
        END SUBROUTINE gser
      end module gser_mod

!
!***********************************************************************
!
      module erf_mod
        contains
          function erf_double_precision(x)  !LAURA
!
!=======================================================================
!
use PRECISION    !LAURA
          implicit none
          integer*8:: n,i
          real(kind=SIZE_DOUBLE):: erf_double_precision,x,deltax,sqrtpi  !LAURA
!
!=======================================================================
!
          sqrtpi=dsqrt(dacos(-1.E+00_8))
          if (dabs(x).ge.3.0E+00_8) then
            erf_double_precision=sign(1.E+00_8,x)
          else
            n=max(1,int(abs(x)/0.001_8))
            deltax=dabs(x/dble(n))
            erf_double_precision=0.0_8
            do i=1,n
              erf_double_precision=erf_double_precision+dexp(-(deltax*(0.5E+00_8+dble(i-1)))**2_CUSTOM_REAL)
            enddo
            erf_double_precision=sign(1.E+00_8,x)*erf_double_precision*deltax*2._8/sqrtpi
          endif
!
!=======================================================================
!
        end function erf_double_precision  !LAURA

          function erf(x)  !LAURA
!
!=======================================================================
!
use PRECISION    !LAURA
          implicit none
          integer*8:: n,i
          real(kind=CUSTOM_REAL):: erf,x,deltax,sqrtpi  !LAURA
!
!=======================================================================
!
          sqrtpi=dsqrt(dacos(-1.E+00_CUSTOM_REAL))
          if (dabs(x).ge.3.0E+00_CUSTOM_REAL) then
            erf=sign(1.E+00_CUSTOM_REAL,x)
          else
            n=max(1,int(dabs(x)/0.001_CUSTOM_REAL))
            deltax=dabs(x/dble(n))
            erf=0.0_CUSTOM_REAL
            do i=1,n
              erf=erf+dexp(-(deltax*(0.5E+00_CUSTOM_REAL+real((i-1),KIND=CUSTOM_REAL)))**2_CUSTOM_REAL)
            enddo
            erf=sign(1.E+00_CUSTOM_REAL,x)*erf*deltax*2._CUSTOM_REAL/sqrtpi
          endif
!
!=======================================================================
!
        end function erf  !LAURA
          function erf_single_precision(x)  !LAURA
!
!=======================================================================
!
use PRECISION    !LAURA
          implicit none  !LAURA
          integer*8:: n,i  !LAURA
          real(kind=SIZE_REAL):: erf_single_precision,x,deltax,sqrtpi  !LAURA
!     
!=======================================================================
!
          sqrtpi=sqrt(acos(-1.E+00_4))  !LAURA
          if (abs(x).ge.3.0E+00_4) then  !LAURA
            erf_single_precision=sign(1.E+00_4,x)  !LAURA
          else  !LAURA
            n=max(1,int(abs(x)/0.001_4))  !LAURA
            deltax=abs(x/REAL(n))  !LAURA
            erf_single_precision=0.0_4  !LAURA
            do i=1,n  !LAURA
              erf_single_precision=erf_single_precision+exp(-(deltax*(0.5E+00_4+REAL(i-1)))**2)  !LAURA
            enddo  !LAURA
            erf_single_precision=sign(1.E+00_4,x)*erf_single_precision*deltax*2._4/sqrtpi  !LAURA
          endif  !LAURA
!
!=======================================================================
!
        end function erf_single_precision  !LAURA
      end module erf_mod
!
!***********************************************************************
!
      module sinc_mod
        contains
          function sinc(x)
!
!=======================================================================
!
use PRECISION    !LAURA
          implicit none
          real(kind=CUSTOM_REAL):: sinc,x	!LAURA
!
!=======================================================================
!
          if (x.eq.0.0) then
            sinc=1.
          else
            sinc=dsin(x)/x
          endif
!
!=======================================================================
!
        end function sinc
      end module sinc_mod
!
!***********************************************************************
!
      module zp_mod
        contains
          double complex function zp(u)
use PRECISION    !LAURA
          implicit none
          double complex:: u,z,u2,azp,azpold,usqm
          integer*8:: n,na
!
!=======================================================================
!
          na=10
          if(abs(u).ge.5.0) go to 3
          usqm=-u**2
          if (dble(usqm).lt.-100.) then
            usqm=cmplx(-100.E+00_SIZE_DOUBLE,0.E+00_SIZE_DOUBLE)+0.5_SIZE_DOUBLE*(usqm-conjg(usqm))	!LAURA
          else if (dble(usqm).gt.+100.) then
            usqm=cmplx(+100.E+00_SIZE_DOUBLE,0.E+00_SIZE_DOUBLE)+0.5_SIZE_DOUBLE*(usqm-conjg(usqm))	!LAURA
          endif
          zp=cmplx(0.E+00_SIZE_DOUBLE,1.E+00_SIZE_DOUBLE)*1.772453850905516027298167_SIZE_DOUBLE*exp(usqm)	!LAURA
          u2=-2._SIZE_DOUBLE*u**2_SIZE_DOUBLE
          azp=-2._SIZE_DOUBLE*u
          do 2 n=1,100
          zp=zp+azp
    2     azp=azp*u2/(2.*n+1.)
          zp=zp+azp
          go to 11
   3      z=1._SIZE_DOUBLE/u
          if(aimag(u).le.(0.)) go to 10
          zp=0.
          go to 20
   10     continue
          usqm=-u**2
          if (dble(usqm).lt.-100.) then
            usqm=cmplx(-100._SIZE_DOUBLE,0._SIZE_DOUBLE)+0.5_SIZE_DOUBLE*(usqm-conjg(usqm))
          else if (dble(usqm).gt.+100.) then
            usqm=cmplx(+100._SIZE_DOUBLE,0._SIZE_DOUBLE)+0.5_SIZE_DOUBLE*(usqm-conjg(usqm))
          endif
          zp=cmplx(0._SIZE_DOUBLE,1._SIZE_DOUBLE)*1.772453850905516027298167_SIZE_DOUBLE*exp(usqm)
          if(aimag(u).lt.0.) zp=2._SIZE_DOUBLE*zp
   20     azp=z
          u2=.5_SIZE_DOUBLE*z**2_SIZE_DOUBLE
          do 25 n=1,na
          zp=zp-azp
          azpold=azp
          azp=(2._SIZE_DOUBLE*n-1._SIZE_DOUBLE)*azp*u2
          if (abs(azp) .ge. abs(azpold)) go to 11
   25     continue
          zp=zp-azp
   11     continue
!
!=======================================================================
!
        end function zp
      end module zp_mod
!
!***********************************************************************
!
      module zprime_mod
        contains
          double complex function zprime(u)
          use zp_mod
use PRECISION    !LAURA
          implicit none
          double complex:: u
          zprime=-2._SIZE_DOUBLE*(1._SIZE_DOUBLE+u*zp(u))
        end function zprime
      end module zprime_mod
!
!***********************************************************************
!
      module z2prime_mod
        contains
          double complex function z2prime(u)
          use zp_mod
          use zprime_mod
use PRECISION    !LAURA
          implicit none
          double complex:: u
          z2prime=-2._SIZE_DOUBLE*(u*zprime(u)+zp(u))
        end function z2prime
      end module z2prime_mod
!
!***********************************************************************
!
      module f_mod
      contains
        double complex function f(omega,theta)
        use zprime_mod
use PRECISION    !LAURA
        implicit none
        real(kind=CUSTOM_REAL):: theta	!LAURA
        double complex:: omega
!
!=======================================================================
!
        f=1._CUSTOM_REAL-0.5_CUSTOM_REAL*theta*zprime(omega)
!
!=======================================================================
!
        end function f
      end module f_mod
!
!***********************************************************************
!
      module fprime_mod
      contains
        double complex function fprime(omega,theta)
        use z2prime_mod
use PRECISION    !LAURA
        implicit none
        real(kind=CUSTOM_REAL):: theta	!LAURA
        double complex:: omega
!
!=======================================================================
!
        fprime=-0.5_CUSTOM_REAL*theta*z2prime(omega)
!
!=======================================================================
!
        end function fprime
      end module fprime_mod
!
!*********************************************************************
!
      module flux_mod
      contains
        function flux(x,argum,vtherm,xran)
!
!=====================================================================
!
        use erf_mod
use PRECISION    !LAURA
        implicit none
        real(kind=CUSTOM_FLUX):: flux,argum,vtherm,x,x1,ss,xran,rat,sqrtpi	!LAURA
!
!=====================================================================
!
!  this function finds the F(v) whose zero will generate
!  the correct flux of particles from the left wall
!
!=====================================================================
!
        sqrtpi=dsqrt(dacos(-1.E+00_CUSTOM_FLUX))	!LAURA
        rat=argum/vtherm
        x1=(x-argum)/vtherm
!        if(CUSTOM_REAL == SIZE_REAL) then !LAURA
!            ss=vtherm*exp(-rat*rat)+argum*sqrtpi*(1.E+00_CUSTOM_REAL+erf_single_precision(rat))	!LAURA
!            flux=ss*(1.E+00_CUSTOM_REAL-xran) &	!LAURA
!                +argum*sqrtpi*(erf_single_precision(x1)-1.E+00_CUSTOM_REAL)-vtherm*exp(-x1*x1)	!LAURA
!        ELSE
            ss=vtherm*dexp(-rat*rat)+argum*sqrtpi*(1.E+00_CUSTOM_REAL+erf_double_precision(dble(rat)))
            flux=ss*(1.E+00_CUSTOM_FLUX-xran) &	!LAURA
               +argum*sqrtpi*(erf_double_precision(dble(x1))-1.E+00_CUSTOM_FLUX)-vtherm*dexp(-x1*x1)
!        END IF 								!LAURA
!
!=====================================================================
!
        end function flux
      end module flux_mod
!
!***********************************************************************
!
      module functions_f90
      contains
!
!***********************************************************************
!
        function sqrnoise(rkx,rky,netot,nitot &
        ,rkdesqr,rkdisqr,dx,dy)
!
!=======================================================================
!
        use sinc_mod
use PRECISION    !LAURA
        implicit none
        integer*8:: netot,nitot
        real(kind=CUSTOM_REAL):: sqrnoise,rksqr,tmpx,rkxhat,tmpy,rkyhat,rkhatsqr,s,ssqr &	!LAURA
        ,chiebar,chie,chii,epsbar,rkx,rky,rkdesqr,rkdisqr,dx,dy
!
!=======================================================================
!
        rksqr=rkx**2_CUSTOM_REAL+rky**2_CUSTOM_REAL
        tmpx=0.5_CUSTOM_REAL*rkx*dx
        rkxhat=rkx*sinc(tmpx)
        tmpy=0.5_CUSTOM_REAL*rky*dy
        rkyhat=rky*sinc(tmpy)
        rkhatsqr=rkxhat**2_CUSTOM_REAL+rkyhat**2_CUSTOM_REAL
        s=sinc(tmpx)**3_CUSTOM_REAL*sinc(tmpy)**3_CUSTOM_REAL
        ssqr=s**2_CUSTOM_REAL
        chiebar= rkdesqr/rkhatsqr
        chie   =(rkdesqr/rkhatsqr)*ssqr
        chii   =(rkdisqr/rkhatsqr)*ssqr
        epsbar =1._CUSTOM_REAL+chiebar+chii
        sqrnoise=ssqr*((rksqr/rkdesqr+1._CUSTOM_REAL/(1._CUSTOM_REAL+chie))/netot          &
        +(1._CUSTOM_REAL-ssqr)**2_CUSTOM_REAL*chiebar**2_CUSTOM_REAL/(nitot*(1._CUSTOM_REAL+chie)**2_CUSTOM_REAL*(1._CUSTOM_REAL+chiebar)  &
                                  *epsbar                        ) &
                                 )*chiebar**2_CUSTOM_REAL
!
!=======================================================================
!
      end function sqrnoise
!
!***********************************************************************
!
      function bessj0(x)
!
!=======================================================================
!
use PRECISION    !LAURA
      implicit none
      real(kind=CUSTOM_REAL):: ax,xx,z,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3 &	!LAURA
      ,r4,r5,r6,s1,s2,s3,s4,s5,s6,y,x,bessj0
      save p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4 &
      ,r5,r6,s1,s2,s3,s4,s5,s6
      data p1,p2,p3,p4,p5/1.e0_CUSTOM_REAL,-.1098628627e-2_CUSTOM_REAL,.2734510407e-4_CUSTOM_REAL, &
         -.2073370639e-5_CUSTOM_REAL,.2093887211e-6_CUSTOM_REAL/
      data  q1,q2,q3,q4,q5/-.1562499995e-1_CUSTOM_REAL, &
         .1430488765e-3_CUSTOM_REAL,-.6911147651E-5_CUSTOM_REAL, &
	.7621095161e-6_CUSTOM_REAL,-.934945152e-7_CUSTOM_REAL/	!LAURA
      data r1,r2,r3,r4,r5,r6/57568490574.e0_CUSTOM_REAL,-13362590354.e0_CUSTOM_REAL, &
          651619640.7e0_CUSTOM_REAL, &
         -11214424.18e0_CUSTOM_REAL,77392.33017e0_CUSTOM_REAL,-184.9052456e0_CUSTOM_REAL/
      data s1,s2,s3,s4,s5,s6/57568490411.e0_CUSTOM_REAL,1029532985.e0_CUSTOM_REAL, &
         9494680.718e0_CUSTOM_REAL,59272.64853e0_CUSTOM_REAL,267.8532712e0_CUSTOM_REAL,1.e0_CUSTOM_REAL/
!
!=======================================================================
!
      if (dabs(x).lt.8._CUSTOM_REAL) then
        y=x**2_CUSTOM_REAL
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
        /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=dabs(x)
        z=8._CUSTOM_REAL/ax
        y=z**2_CUSTOM_REAL
        xx=ax-.785398164_CUSTOM_REAL
        bessj0=dsqrt(.636619772_CUSTOM_REAL/ax)                &
        *dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))   &
        -z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5))))
      endif
!
!=======================================================================
!
      end function bessj0
!
!***********************************************************************
!
      double complex function shape(xix,xiy,f,rk0,xhs,yhs,xleft &
                                    ,plane_wave)
!
!======================================================================
!
use PRECISION    !LAURA
      implicit none
      logical:: plane_wave
      real(kind=CUSTOM_REAL):: xix,xiy,f,rk0,xhs,yhs,xleft	!LAURA
      double complex:: sigma2
!
!======================================================================
!
      if (plane_wave) then
        shape=1._CUSTOM_REAL
      else
        sigma2=(2._CUSTOM_REAL*f/rk0)**2_CUSTOM_REAL+dcmplx(0.E+00_CUSTOM_REAL,1.E+00_CUSTOM_REAL)*(xix-xhs)/(2._CUSTOM_REAL*rk0)	!LAURA
!
!     3D Gaussian beam
!
!      shape=exp(-xiy**2/(4.*sigma2)                   &
!                 +cmplx(0.E+00_CUSTOM_REAL,1.E+00_CUSTOM_REAL)*rk0*(xix-xhs)) &	!LAURA
!            *(2.*f/rk0)**2/sigma2
!
!     2D Gaussian beam
!
        shape=exp(-xiy**2_CUSTOM_REAL/(4._CUSTOM_REAL*sigma2)                   &
                   +cmplx(0.E+00_SIZE_DOUBLE,1.E+00_SIZE_DOUBLE)*rk0*(xix-xhs)) &	!LAURA
              *(2._CUSTOM_REAL*f/rk0)/sqrt(sigma2)
      endif
!
!======================================================================
!
      end function shape
!
!***********************************************************************
!
!      double precision FUNCTION myranf()
!      implicit none
!      INTEGER*8 idum
!      INTEGER*8 MBIG,MSEED,MZ
!      real(kind=CUSTOM_REAL) FAC	!LAURA
!      INTEGER*8 i,iff,ii,inext,inextp,k
!      INTEGER*8 mj,mk,ma(55)
!      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0 &
!                ,FAC=1.E-09_CUSTOM_REAL)	!LAURA
!      COMMON /MYRANDOM/ IDUM
!      SAVE iff,inext,inextp,ma
!      DATA iff /0/
!      if(idum.lt.0.or.iff.eq.0)then
!	iff=1
!        mj=MSEED-iabs(idum)
!	mj=mod(mj,MBIG)
!	ma(55)=mj
!	mk=1
!	do 11 i=1,54
!	  ii=mod(21*i,55)
!	  ma(ii)=mk
!	  mk=mj-mk
!	  if(mk.lt.MZ)mk=mk+MBIG
!	  mj=ma(ii)
!11      continue
!	do 13 k=1,4
!	  do 12 i=1,55
!	    ma(i)=ma(i)-ma(1+mod(i+30,55))
!	    if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
!12        continue
!13      continue
!	inext=0
!	inextp=31
!	idum=1
!      endif
!      inext=inext+1
!      if(inext.eq.56)inext=1
!      inextp=inextp+1
!      if(inextp.eq.56)inextp=1
!      mj=ma(inext)-ma(inextp)
!      if(mj.lt.MZ)mj=mj+MBIG
!      ma(inext)=mj
!      myranf=dble(mj)*FAC
!      END FUNCTION myranf
!
!***********************************************************************
!
      function fmaxwell(vx,vy,vthe)
!
!=======================================================================
!
use PRECISION    !LAURA
      implicit none
      real(kind=CUSTOM_REAL):: fmaxwell,vx,vy,vthe,twopi	!LAURA
!
!=======================================================================
!
      twopi=2._CUSTOM_REAL*dacos(-1.E+00_CUSTOM_REAL)	!LAURA
      fmaxwell=dexp(-(vx**2_CUSTOM_REAL+vy**2_CUSTOM_REAL)/(2._CUSTOM_REAL*vthe**2_CUSTOM_REAL))/(twopi*vthe**2_CUSTOM_REAL)
!
!=======================================================================
!
      end function fmaxwell
!
!***********************************************************************
!
      real(kind=8) function fmaxwell1d(vx,vy,vthe)
use PRECISION    !LAURA
      implicit none !LAURA
!
!=======================================================================
!
      real(kind=CUSTOM_REAL):: vx,vthe,vy,sqrttwopi	!LAURA
!
!=======================================================================
!
      sqrttwopi=dsqrt(dacos(-1.E+00_CUSTOM_REAL))	!LAURA
      fmaxwell1d=dexp(-vx**2_CUSTOM_REAL/(2._CUSTOM_REAL*vthe**2_CUSTOM_REAL))/(sqrttwopi*vthe)
!
!=======================================================================
!
      end function fmaxwell1d
!
!***********************************************************************
!
      FUNCTION gammp(a,x)
      use gser_mod
      use gcf_mod
use PRECISION    !LAURA
      implicit none !LAURA
      real(kind=CUSTOM_REAL) gammp,a,x	!LAURA
!U    USES gcf,gser
      real(kind=CUSTOM_REAL) gammcf,gamser,gln	!LAURA
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1._CUSTOM_REAL-gammcf
      endif
      END FUNCTION gammp
!
!***********************************************************************
!
      end module functions_f90
