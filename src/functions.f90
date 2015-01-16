!
!***********************************************************************
!
      module gammaln_mod
      contains
        DOUBLE PRECISION FUNCTION gammln(xx)
        implicit none
        DOUBLE PRECISION xx
        INTEGER*8 j
        DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
        SAVE cof,stp
        DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
        24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
        -.5395239384953d-5,2.5066282746310005d0/
        x=xx
        y=x
        tmp=x+5.5d0
        tmp=(x+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0
        do 11 j=1,6
          y=y+1.d0
          ser=ser+cof(j)/y
11      continue
        gammln=tmp+log(stp*ser/x)
        END FUNCTION gammln
      end module gammaln_mod
!
!***********************************************************************
!
      module gcf_mod
      contains
        SUBROUTINE gcf(gammcf,a,x,gln)
        use gammaln_mod
        implicit none
        INTEGER*8 ITMAX
        DOUBLE PRECISION a,gammcf,gln,x,EPS,FPMIN
        PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!U      USES gammln
        INTEGER*8 i
        DOUBLE PRECISION an,b,c,d,del,h
        gln=gammln(a)
        b=x+1.-a
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,ITMAX
          an=-i*(i-a)
          b=b+2.
          d=an*d+b
          if(abs(d).lt.FPMIN)d=FPMIN
          c=b+an/c
          if(abs(c).lt.FPMIN)c=FPMIN
          d=1./d
          del=d*c
          h=h*del
          if(abs(del-1.).lt.EPS)goto 1
11      continue
        pause 'a too large, ITMAX too small in gcf'
1       gammcf=exp(-x+a*log(x)-gln)*h
        END SUBROUTINE gcf
      end module gcf_mod
!
!***********************************************************************
!
      module gser_mod
      contains
        SUBROUTINE gser(gamser,a,x,gln)
        use gammaln_mod
        implicit none
        INTEGER*8 ITMAX
        DOUBLE PRECISION a,gamser,gln,x,EPS
        PARAMETER (ITMAX=100,EPS=3.e-7)
!U      USES gammln
        INTEGER*8 n
        DOUBLE PRECISION ap,del,sum
        gln=gammln(a)
        if(x.le.0.)then
          if(x.lt.0.)pause 'x < 0 in gser'
          gamser=0.
          return
        endif
        ap=a
        sum=1./a
        del=sum
        do 11 n=1,ITMAX
          ap=ap+1.
          del=del*x/ap
          sum=sum+del
          if(abs(del).lt.abs(sum)*EPS)goto 1
11      continue
        pause 'a too large, ITMAX too small in gser'
1       gamser=sum*exp(-x+a*log(x)-gln)
        END SUBROUTINE gser
      end module gser_mod

!
!***********************************************************************
!
      module erf_mod
        contains
          double precision function erf(x)
!
!=======================================================================
!
          implicit none
          integer*8:: n,i
          double precision:: x,deltax,sqrtpi
!
!=======================================================================
!
          sqrtpi=sqrt(acos(-1.d+00))
          if (abs(x).ge.3.0d+00) then
            erf=sign(1.d+00,x)
          else
            n=max(1,int(abs(x)/0.001))
            deltax=abs(x/dble(n))
            erf=0.0
            do i=1,n
              erf=erf+exp(-(deltax*(0.5d+00+dble(i-1)))**2)
            enddo
            erf=sign(1.d+00,x)*erf*deltax*2./sqrtpi
          endif
!
!=======================================================================
!
        end function erf
      end module erf_mod
!
!***********************************************************************
!
      module sinc_mod
        contains
          double precision function sinc(x)
!
!=======================================================================
!
          implicit none
          double precision:: x
!
!=======================================================================
!
          if (x.eq.0.0) then
            sinc=1.
          else
            sinc=sin(x)/x
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
          usqm=cmplx(-100.d+00,0.d+00)+0.5*(usqm-conjg(usqm))
          else if (dble(usqm).gt.+100.) then
          usqm=cmplx(+100.d+00,0.d+00)+0.5*(usqm-conjg(usqm))
          endif
          zp=cmplx(0.d+00,1.d+00)*1.772453850905516027298167*exp(usqm)
          u2=-2.*u**2
          azp=-2.*u
          do 2 n=1,100
          zp=zp+azp
    2     azp=azp*u2/(2.*n+1.)
          zp=zp+azp
          go to 11
   3      z=1./u
          if(aimag(u).le.(0.)) go to 10
          zp=0.
          go to 20
   10     continue
          usqm=-u**2
          if (dble(usqm).lt.-100.) then
          usqm=cmplx(-100.,0.)+0.5*(usqm-conjg(usqm))
          else if (dble(usqm).gt.+100.) then
          usqm=cmplx(+100.,0.)+0.5*(usqm-conjg(usqm))
          endif
          zp=cmplx(0.,1.)*1.772453850905516027298167*exp(usqm)
          if(aimag(u).lt.0.) zp=2.*zp
   20     azp=z
          u2=.5*z**2
          do 25 n=1,na
          zp=zp-azp
          azpold=azp
          azp=(2.*n-1.)*azp*u2
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
          implicit none
          double complex:: u
          zprime=-2.*(1.+u*zp(u))
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
          implicit none
          double complex:: u
          z2prime=-2.*(u*zprime(u)+zp(u))
        end function z2prime
      end module z2prime_mod
!
!***********************************************************************
!
      module f_mod
      contains
        double complex function f(omega,theta)
        use zprime_mod
        implicit none
        double precision:: theta
        double complex:: omega
!
!=======================================================================
!
        f=1.-0.5*theta*zprime(omega)
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
        implicit none
        double precision:: theta
        double complex:: omega
!
!=======================================================================
!
        fprime=-0.5*theta*z2prime(omega)
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
        double precision function flux(x,argum,vtherm,xran,ss)
!
!=====================================================================
!
        use erf_mod
        implicit none
        double precision:: argum,vtherm,x,x1,ss,xran,rat,sqrtpi
!
!=====================================================================
!
!  this function finds the F(v) whose zero will generate
!  the correct flux of particles from the left wall
!
!=====================================================================
!
        sqrtpi=sqrt(acos(-1.d+00))
        rat=argum/vtherm
        x1=(x-argum)/vtherm
        flux=ss*(1.d+00-xran) &
        +argum*sqrtpi*(erf(x1)-1.d+00)-vtherm*exp(-x1*x1)
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
        double precision function sqrnoise(rkx,rky,netot,nitot &
        ,rkdesqr,rkdisqr,dx,dy)
!
!=======================================================================
!
        use sinc_mod
        implicit none
        integer*8:: netot,nitot
        double precision:: rksqr,tmpx,rkxhat,tmpy,rkyhat,rkhatsqr,s,ssqr &
        ,chiebar,chie,chii,epsbar,rkx,rky,rkdesqr,rkdisqr,dx,dy
!
!=======================================================================
!
        rksqr=rkx**2+rky**2
        tmpx=0.5*rkx*dx
        rkxhat=rkx*sinc(tmpx)
        tmpy=0.5*rky*dy
        rkyhat=rky*sinc(tmpy)
        rkhatsqr=rkxhat**2+rkyhat**2
        s=sinc(tmpx)**3*sinc(tmpy)**3
        ssqr=s**2
        chiebar= rkdesqr/rkhatsqr
        chie   =(rkdesqr/rkhatsqr)*ssqr
        chii   =(rkdisqr/rkhatsqr)*ssqr
        epsbar =1.+chiebar+chii
        sqrnoise=ssqr*((rksqr/rkdesqr+1./(1.+chie))/netot          &
        +(1.-ssqr)**2*chiebar**2/(nitot*(1.+chie)**2*(1.+chiebar)  &
                                  *epsbar                        ) &
                                 )*chiebar**2
!
!=======================================================================
!
      end function sqrnoise
!
!***********************************************************************
!
      double precision function bessj0(x)
!
!=======================================================================
!
      implicit none
      double precision:: ax,xx,z,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3 &
      ,r4,r5,r6,s1,s2,s3,s4,s5,s6,y,x
      save p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4 &
      ,r5,r6,s1,s2,s3,s4,s5,s6
      data p1,p2,p3,p4,p5/1.e0,-.1098628627e-2,.2734510407e-4, &
         -.2073370639e-5,.2093887211e-6/
      data  q1,q2,q3,q4,q5/-.1562499995e-1, &
         .1430488765e-3,-.6911147651d-5,.7621095161e-6,-.934945152e-7/
      data r1,r2,r3,r4,r5,r6/57568490574.e0,-13362590354.e0, &
          651619640.7e0, &
         -11214424.18e0,77392.33017e0,-184.9052456e0/
      data s1,s2,s3,s4,s5,s6/57568490411.e0,1029532985.e0, &
         9494680.718e0,59272.64853e0,267.8532712e0,1.e0/
!
!=======================================================================
!
      if (abs(x).lt.8.) then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
        /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)                &
        *cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))   &
        -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5))))
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
      implicit none
      logical:: plane_wave
      double precision:: xix,xiy,f,rk0,xhs,yhs,xleft
      double complex:: sigma2
!
!======================================================================
!
      if (plane_wave) then
        shape=1.
      else
        sigma2=(2.*f/rk0)**2+cmplx(0.d+00,1.d+00)*(xix-xhs)/(2.*rk0)
!
!     3D Gaussian beam
!
!      shape=exp(-xiy**2/(4.*sigma2)                   &
!                 +cmplx(0.d+00,1.d+00)*rk0*(xix-xhs)) &
!            *(2.*f/rk0)**2/sigma2
!
!     2D Gaussian beam
!
        shape=exp(-xiy**2/(4.*sigma2)                   &
                   +cmplx(0.d+00,1.d+00)*rk0*(xix-xhs)) &
              *(2.*f/rk0)/sqrt(sigma2)
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
!      double precision FAC
!      INTEGER*8 i,iff,ii,inext,inextp,k
!      INTEGER*8 mj,mk,ma(55)
!      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0 &
!                ,FAC=1.d-09)
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
      double precision function fmaxwell(vx,vy,vthe)
!
!=======================================================================
!
      implicit none
      double precision:: vx,vy,vthe,twopi
!
!=======================================================================
!
      twopi=2.*acos(-1.d+00)
      fmaxwell=exp(-(vx**2+vy**2)/(2.*vthe**2))/(twopi*vthe**2)
!
!=======================================================================
!
      end function fmaxwell
!
!***********************************************************************
!
      double precision function fmaxwell1d(vx,vy,vthe)
!
!=======================================================================
!
      double precision:: vx,vthe,vy,sqrttwopi
!
!=======================================================================
!
      sqrttwopi=sqrt(acos(-1.d+00))
      fmaxwell1d=exp(-vx**2/(2.*vthe**2))/(sqrttwopi*vthe)
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
      DOUBLE PRECISION a,gammp,x
!U    USES gcf,gser
      DOUBLE PRECISION gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      END FUNCTION gammp
!
!***********************************************************************
!
      end module functions_f90
