pro gimage,z,x,y,nxmin,nxmax,nymin,nymax,color_table,xtitle,ytitle,ptitle,panel_position,erase,mmax,mmin,clevel_1,clevel_2
!p.title='!6'
;jpeg
;!p.charsize=2.0
;postscript
!p.charsize=1.0
;
!x.thick=3.
!y.thick=3.
;mmax=max(z(nxmin:nxmax,nymin:nymax))
;
; Inverted color scale
;
;icolor=long(bytscl(mmax-z(nxmin:nxmax,nymin:nymax),/NaN))
;
;
; Regular color scale
;
icolor=long(bytscl(z(nxmin:nxmax,nymin:nymax),max=mmax,min=mmin,/NaN))
;
;
b=size(z)
;!z.range(0)=min(z(nxmin:nxmax,nymin:nymax))
;!z.range(1)=max(z(nxmin:nxmax,nymin:nymax))
!z.range(0)=mmin
!z.range(1)=mmax
!x.margin=[15,15]
!y.margin=[9,9]
kx_min=x(nxmin)
kx_max=x(nxmax)
ky_min=y(nymin)
ky_max=y(nymax)
!x.range=[kx_min,kx_max]
!y.range=[ky_min,ky_max]
!x.style=1
!y.style=1
xsc=fltarr(2)
ysc=fltarr(2)
xsc(0)=kx_min
xsc(1)=kx_max
ysc(0)=ky_min
ysc(1)=ky_max
if (!p.charsize eq 0) then !p.charsize=1.
!x.title=xtitle
!y.title=ytitle
if (erase eq 0) then begin
;  plot,xsc,ysc,/nodata,position=panel_position
  plot,xsc,ysc,/nodata
endif else begin
;  plot,xsc,ysc,/nodata,position=panel_position,/noerase
  plot,xsc,ysc,/nodata,/noerase
endelse
loadct,color_table
;
; Reverse color table
;
TVLCT,R,G,B,/GET
RR=REVERSE(R)
GG=REVERSE(G)
BB=REVERSE(B)
TVLCT,RR,GG,BB
;
;
xv=fltarr(4)
yv=fltarr(4)
for i=nxmin,nxmax-1 do begin
for j=nymin,nymax-1 do begin
xv(0)=x(i)
yv(0)=y(j)
xv(1)=x(i+1)
yv(1)=y(j)
xv(2)=x(i+1)
yv(2)=y(j+1)
xv(3)=x(i)
yv(3)=y(j+1)
if (z(i,j) ge mmin) then begin
polyfill,xv,yv,color=icolor(i-nxmin,j-nymin)
endif
endfor
endfor

height=0.025*(ysc(1)-ysc(0))*!p.charsize
;plot,xsc,ysc,/nodata,/noerase,position=panel_position
plot,xsc,ysc,/nodata,/noerase
loadct,color_table
;
; Reverse color table
;
TVLCT,R,G,B,/GET
RR=REVERSE(R)
GG=REVERSE(G)
BB=REVERSE(B)
TVLCT,RR,GG,BB

xscale=fltarr(257,2)
yscale=fltarr(257,2)
xscale(*,0)=xsc(0)+(xsc(1)-xsc(0))*findgen(257)/256.
yscale(*,0)=ysc(1)+1.5*height
xscale(*,1)=xscale(*,0)
yscale(*,1)=ysc(1)+2.5*height
;
; Inverted color scale
;
;icolor=255-indgen(256)
;
;
; Regular color scale
;
icolor=indgen(256)
;
;
for i=0,255 do begin
xv(0)=xscale(i,0)
yv(0)=yscale(i,0)
xv(1)=xscale(i+1,0)
yv(1)=yscale(i+1,0)
xv(2)=xscale(i+1,1)
yv(2)=yscale(i+1,1)
xv(3)=xscale(i,1)
yv(3)=yscale(i,1)
polyfill,xv,yv,color=icolor(i)
endfor
text='!6           '
strput,text,strcompress(string(!z.range(0),/print,format='(e10.2)'),/remove_all),3
xyouts,xsc(0),yscale(0,0),text,size=0.75*!p.charsize,alignment=1.0
text='!6           '
strput,text,strcompress(string(!z.range(1),/print,format='(e10.2)'),/remove_all),3
xyouts,xsc(1),yscale(0,0),text,size=0.75*!p.charsize,alignment=0.0
;text='!10#!6n!deS!n(k!dx!n,k!dy!n)!10#!6!u2!n'
text=ptitle
xyouts,0.5*(xsc(0)+xsc(1)),yscale(0,0)+3.*height,text,size=!p.charsize,alignment=0.5

xoffset = x + 0.5*(x(1)-x(0))
yoffset = y + 0.5*(y(1)-y(0))
!p.thick=2
clevel_1=clevel_1*max(z)
clevel_2=clevel_2*max(z)
;if (clevel_1 ge 0.) then contour,z,xoffset,yoffset,levels=clevel_1,position=panel_position,/noerase,c_color=255
if (clevel_1 ge 0.) then contour,z,xoffset,yoffset,levels=clevel_1,/noerase,c_color=255
!p.linestyle=1
;if (clevel_2 ge 0.) then contour,z,xoffset,yoffset,levels=clevel_2,position=panel_position,/noerase,c_color=255
if (clevel_2 ge 0.) then contour,z,xoffset,yoffset,levels=clevel_2,/noerase,c_color=255
!p.linestyle=0
end
