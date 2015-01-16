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

