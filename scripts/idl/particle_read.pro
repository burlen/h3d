v_scale=6000.
nv = 100
nx = 100


device,decompose=0,true=24

;set_plot,'ps'
;device,/in,xsize=5.,ysize=5.,/color,file='particle.ps'

!p.charsize=1.2
!y.charsize=1.2
!z.charsize=1.2
;openr,1,'/nobackupp8/hkarimab/hkarimab_nbp1/run8b/data/particle_400000.bin'
openr,1,'/nobackupp8/hkarimab/hkarimab_nbp1/run8i_res/data/particle_110000.bin'
n_total=0l
n_count=0l

readu,1,n_total
n_count_down = n_total
print,"line_1, n_total = ",n_total
a=fltarr(9)
readu,1,a
xp=fltarr(n_total)
yp=fltarr(n_total)
zp=fltarr(n_total)
up=fltarr(n_total)
vp=fltarr(n_total)
wp=fltarr(n_total)
qp=fltarr(n_total)
vpar=fltarr(n_total)
vperp_1=fltarr(n_total)
vperp_2=fltarr(n_total)
ip_count = -1l

repeat begin
readu,1,n_count
print," line_2, n_count = ",n_count
a=fltarr(9)
readu,1,a
b=fltarr(10)
for ip=1,n_count do begin
  ip_count = ip_count + 1l
  readu,1,b
  xp(ip_count)      = b(0)
  yp(ip_count)      = b(1)
  zp(ip_count)      = b(2)
  up(ip_count)      = b(3) * v_scale
  vp(ip_count)      = b(4) * v_scale
  wp(ip_count)      = b(5) * v_scale
  qp(ip_count)      = b(6)
  vpar(ip_count)    = b(7) * v_scale
  vperp_1(ip_count) = b(8) * v_scale
  vperp_2(ip_count) = b(9) * v_scale
endfor
n_count_down = n_count_down - n_count
print,"line_3, n_count_down = ",n_count_down
endrep until n_count_down eq 0
close,1
;
; Reverse color table
;
TVLCT,R,G,B,/GET
RR=REVERSE(R)
GG=REVERSE(G)
BB=REVERSE(B)
TVLCT,RR,GG,BB
;
; The particle arrays are xp,yp,zp,up,vp,wp,qp.  Each array is of size n_total
;plot,vpar,sqrt(vperp_1^2+vperp_2^2),psym=4,symsize=0.1
;
; f(vpar)
;
vpar_min= min(vpar)
vpar_max= max(vpar)
nv_par  = nv 
binsize_vpar=(vpar_max-vpar_min)/(nv_par-1)
x_min   = min(xp)
x_max   = max(xp)
dx      = (x_max - x_min)/(nx-1)
x       = x_min + dx*findgen(nx)
vpar_arr=vpar_min+(findgen(nv_par) + 1.0)*binsize_vpar
f_vpar=histogram(vpar,binsize=binsize_vpar)
f_vpar=f_vpar/(total(f_vpar)*binsize_vpar)

f_vpar=fltarr(nx,nv_par)
for ip=0l,n_total-1l do begin
index_x   =(xp(ip)-x_min)/dx
index_vpar=(vpar(ip)-vpar_min)/binsize_vpar
if (index_x ge 0 and index_x le nx-1l and index_vpar ge 0 and index_vpar le nv_par-1l) then begin
f_vpar(index_x,index_vpar) = f_vpar(index_x,index_vpar) + 1l
endif
endfor
f_vpar = smooth(f_vpar,3)
f_vpar = f_vpar / (total(f_vpar)*binsize_vpar*dx)

;
; f(vperp_1)
;
vperp_1_min=min(vperp_1)
vperp_1_max=max(vperp_1)
nv_perp_1  = nv 
binsize_vperp_1=(vperp_1_max-vperp_1_min)/(nv_perp_1-1)
vperp_1_arr=vperp_1_min+(findgen(nv_perp_1) + 1.0)*binsize_vperp_1
f_vperp_1=histogram(vperp_1,binsize=binsize_vperp_1)
f_vperp_1=f_vperp_1/(total(f_vperp_1)*binsize_vperp_1)
!p.title='!6'
!x.title='!6v!dperp 1!n/v!dA!n'
!y.title='!6f(v!dperp 1!n/v!dA!n)'

window,0,xsize=400,ysize=400,xpos=0,ypos=500
wset,0

plot,vperp_1_arr,f_vperp_1
;
; f(vperp_2)
;
vperp_2_min=min(vperp_2)
vperp_2_max=max(vperp_2)
nv_perp_2  = nv 
binsize_vperp_2=(vperp_2_max-vperp_2_min)/(nv_perp_2-1)
vperp_2_arr=vperp_2_min+(findgen(nv_perp_2)+1.0)*binsize_vperp_2
f_vperp_2=histogram(vperp_2,binsize=binsize_vperp_2)
f_vperp_2=f_vperp_2/(total(f_vperp_2)*binsize_vperp_2)
!p.title='!6'
!x.title='!6v!dperp 2!n/v!dA!n'
!y.title='!6f(v!dperp 2!n/v!dA!n)'

window,1,xsize=400,ysize=400,xpos=500,ypos=500
wset,1

plot,vperp_2_arr,f_vperp_2
;
; f(vperp)
;
vperp = sqrt (vperp_1^2 + vperp_2^2)
vperp_min=min(vperp)
vperp_max=max(vperp)
nv_perp  = nv 
binsize_vperp=(vperp_max-vperp_min)/(nv_perp-1)
vperp_arr=vperp_min+(findgen(nv_perp)+1.0)*binsize_vperp
f_vperp=histogram(vperp,binsize=binsize_vperp)
f_vperp=f_vperp/(total(f_vperp)*binsize_vperp)
!p.title='!6'
!x.title='!10"!6v!dperp!n/v!dA!n!10"!6'
!y.title='!6f(!10"!6v!dperp!n/v!dA!n!10"!6)'

window,2,xsize=400,ysize=400,xpos=1000,ypos=500
wset,2

plot,vperp_arr,f_vperp
;
; f(vabs)
;
vabs  = sqrt (vperp^2 + vpar^2)
vabs_min=min(vabs)
vabs_max=max(vabs)
nv_abs  = nv 
binsize_vabs=(vabs_max-vabs_min)/(nv_abs-1)
vabs_arr=vabs_min+(findgen(nv_abs)+1.0)*binsize_vabs
f_vabs=histogram(vabs,binsize=binsize_vabs)
f_vabs=f_vabs/(total(f_vabs)*binsize_vabs)
!p.title='!6'
!x.title='!10"!6v/v!dA!n!10"!6'
!y.title='!6f(!10"!6v/v!dA!n!10"!6)'

window,3,xsize=400,ysize=400,xpos=0,ypos=0
wset,3

plot,vabs_arr,f_vabs
;
;
; compute 2D distribution function f(vperp_1,vperp_2)
;
f_2v=fltarr(nv_perp_1,nv_perp_2)
for ip=0l,n_total-1l do begin
index_1=(vperp_1(ip)-vperp_1_min)/binsize_vperp_1
index_2=(vperp_2(ip)-vperp_2_min)/binsize_vperp_2
if (index_1 ge 0 and index_1 le nv_perp_1-1l and index_2 ge 0 and index_2 le nv_perp_2-1l) then begin
f_2v(index_1,index_2) = f_2v(index_1,index_2) + 1l
endif
endfor
f_2v = smooth(f_2v,3)
f_2v = f_2v / (total(f_2v)*binsize_vperp_1*binsize_vperp_2)
panel_position=[0.15,0.15,0.80,0.80]
xtitle='!6v!dperp 1!n/v!dA!n'
ytitle='!6v!dperp 2!n/v!dA!n'
clevel_1 = 0.5
clevel_2 = 0.1
erase=0
!x.margin=[15,15]
!y.margin=[9,9]
mmax=max(f_2v)
mmin=0.
ptitle='!6f(v!dperp 1!n/v!dA!n,!6v!dperp 2!n/v!dA!n)'
color_table=5

window,4,xsize=400,ysize=400,xpos=500,ypos=0
wset,4

gimage,f_2v,vperp_1_arr,vperp_2_arr,0,nv_perp_1-1,0,nv_perp_2-1,color_table,xtitle,ytitle,ptitle,panel_position,erase,mmax,mmin,clevel_1,clevel_2
;
;
;
panel_position=[0.15,0.15,0.80,0.80]
xtitle='!6x/(c/!7x!6!dpi!n)'
ytitle='!6v!dpar!n/v!dA!n'
clevel_1 = -0.5
clevel_2 = -0.1
erase=0
!x.margin=[15,15]
!y.margin=[9,9]
mmax=max(f_vpar)
mmin=0.
ptitle='!6f(x/(c/!7x!6!dpi!n),!6v!dpar!n/v!dA!n)'
color_table=5
window,5,xsize=400,ysize=400,xpos=1000,ypos=0
wset,5

gimage,f_vpar,x,vpar_arr,1,nx-2,0,nv_par-1,color_table,xtitle,ytitle,ptitle,panel_position,erase,mmax,mmin,clevel_1,clevel_2

;device,/close



end
