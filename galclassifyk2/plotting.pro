; various plotting stuff

pro plotplx,data,model,ug,q,sig,absmag,absmage,hkcol,jhcol,col1s,col2s,bvcol,jkcol,cols,col5s,appmage,jhk

if not (jhk) then begin
    mabsm=model.mabs_v 
    obscol=bvcol[q]
    obscole=col5s
    xn='B-V'
    yn='M_V'
endif else begin
    mabsm=model.mabs_j
    obscol=jkcol[q]
    obscole=cols
    xn='J-K'
    yn='M_J'
endelse

ran=cgrandomindices(long(n_elements(ug)),long(n_elements(ug)*0.1))
ug=ug[ran]

!p.multi=0.
wset,1
!p.charsize=1.0
!p.multi=[0,2,1]

if (jhk) then begin
    plot,[0],[0],psym=3,xrange=[-0.2,1.5],yrange=[15,-10],/xs,$
    xtitle=xn,ytitle=yn		
    
    u=where(model[ug].logg gt 4.0)
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,mabsm[ug[u]],psym=4,symsize=0.2,color=1
    u=where(model[ug].logg gt 3.5 and model[ug].logg lt 4.0)
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,mabsm[ug[u]],psym=4,symsize=0.2,color=2
    u=where(model[ug].logg gt 2.0 and model[ug].logg lt 3.5)
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,mabsm[ug[u]],psym=4,symsize=0.2,color=3
    u=where(model[ug].logg lt 2.0)
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,mabsm[ug[u]],psym=4,symsize=0.2,color=4
    oploterror,[obscol],[absmag],[obscole],[absmage],thick=3,errthick=3

    plot,[0],[0],psym=3,xrange=[-0.2,0.5],yrange=[0.1,1.2],/xs,$
    xtitle='H-K',ytitle='J-H'
    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
    oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
    psym=4,symsize=0.2,color=1
    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
    oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
    psym=4,color=3,symsize=0.2
    oploterror,[hkcol[q]],[jhcol[q]],[col1s],[col2s],thick=3,errthick=3;,color=6,errcolor=6
 endif else begin
    !p.multi=0
    plot,[0],[0],psym=3,xrange=[-0.2,2.5],yrange=[15,-10],/xs,$
    xtitle=xn,ytitle=yn		
    u=where(model[ug].logg gt 4.0)
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,mabsm[ug[u]],psym=4,symsize=0.2,color=1
    u=where(model[ug].logg gt 3.5 and model[ug].logg lt 4.0)
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,mabsm[ug[u]],psym=4,symsize=0.2,color=2
    u=where(model[ug].logg gt 2.0 and model[ug].logg lt 3.5)
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,mabsm[ug[u]],psym=4,symsize=0.2,color=3
    u=where(model[ug].logg lt 2.0)
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,mabsm[ug[u]],psym=4,symsize=0.2,color=4
    oploterror,[obscol],[absmag],[obscole],[absmage],thick=3,errthick=3
endelse

!p.multi=0				
an=''
read,an

end


pro plotspec,model,teff_in,logg_in,teffe_in,logge_in

ug=cgrandomindices(long(n_elements(model.teff)),1e5)

!p.multi=0
wset,1
plot,[0],[0],psym=3,xrange=[8000,3000],yrange=[6,0],/xs,$
xtitle='Teff',ytitle='logg'
u=where(model[ug].logg gt 4.0)
oplot,model[ug[u]].teff,model[ug[u]].logg,psym=4,color=1,symsize=0.2
u=where(model[ug].logg lt 4.0 and model[ug].logg gt 3.5)
oplot,model[ug[u]].teff,model[ug[u]].logg,psym=4,color=2,symsize=0.2
u=where(model[ug].logg lt 3.5 and model[ug].logg gt 2.0)
oplot,model[ug[u]].teff,model[ug[u]].logg,psym=4,color=3,symsize=0.2
u=where(model[ug].logg lt 2.0)
oplot,model[ug[u]].teff,model[ug[u]].logg,psym=4,color=4,symsize=0.2

oploterror,[teff_in],[logg_in],[teffe_in],[logge_in],thick=3,$
errthick=3

an=''
read,an

end


pro plotrpm,data,model,ug,q,redpm,jkcol,hkcol,jhcol,grcol,ricol,cols,col1s,$
col2s,col3s,col4s,sig,appmage,rpm,rpme

!p.multi=0.

wset,1
!p.charsize=2
if (rpm[q] ne 0.) then begin
	!p.multi=[0,1,3]
	plot,[0],[0],psym=3,xrange=[-0.2,1.2],yrange=[10,-10],/xs,$
	xtitle='J-K',ytitle='RPM_J'		

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg gt 4.0)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,redpm[ug[u]],psym=4,symsize=0.2,color=1

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg gt 3.5 and model[ug].logg lt 4.0)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,redpm[ug[u]],psym=4,symsize=0.2,color=2

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg gt 2.0 and model[ug].logg lt 3.5)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,redpm[ug[u]],psym=4,symsize=0.2,color=3

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg lt 2.0)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,redpm[ug[u]],psym=4,symsize=0.2,color=4

	oploterror,[jkcol[q]],[rpm[q]],[cols],[rpme[q]],thick=3,errthick=3
endif else !p.multi=[0,1,2]

plot,[0],[0],psym=3,xrange=[-0.2,0.8],yrange=[0.1,1.5],/xs,$
xtitle='H-K',ytitle='J-H'
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
	psym=4,symsize=0.2,color=1	
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
	psym=4,color=3,symsize=0.2
oploterror,[hkcol[q]],[jhcol[q]],[col1s],[col2s],thick=3,errthick=3;,color=6,errcolor=6

plot,[0],[0],psym=3,/xs,xrange=[0,3.0],yrange=[0,1.5],/ys,$
xtitle='g-r',ytitle='r-i'
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
oplot,model[ug[u]].mapp_g-model[ug[u]].mapp_r,model[ug[u]].mapp_r-model[ug[u]].mapp_i,$
psym=4,symsize=0.2,color=1
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
if (u[0] ne -1) then oplot,model[ug[u]].mapp_g-model[ug[u]].mapp_r,model[ug[u]].mapp_r-model[ug[u]].mapp_i,$
psym=4,color=3,symsize=0.2
oploterror,[grcol[q]],[ricol[q]],[col3s],[col4s],thick=3,errthick=3;,color=6,errcolor=6
!p.multi=0		

an=''
read,an

end


pro plotpm,data,model,ug,q,redpm,jkcol,hkcol,jhcol,grcol,ricol,cols,col1s,$
col2s,col3s,col4s,sig,appmage,pm,pme

!p.multi=0.

wset,1
!p.charsize=2
if (pm ne 0.) then begin
	!p.multi=[0,1,3]
	plot,[0],[0],psym=3,xrange=[-0.2,2.0],yrange=[1000,0.01],/xs,$
	xtitle='J-K',ytitle='Proper Motion',/ylog		

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg gt 4.0)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,model[ug[u]].propmotion,psym=4,symsize=0.2,color=1

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg gt 3.5 and model[ug].logg lt 4.0)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,model[ug[u]].propmotion,psym=4,symsize=0.2,color=2

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg gt 2.0 and model[ug].logg lt 3.5)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,model[ug[u]].propmotion,psym=4,symsize=0.2,color=3

	u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
	data[q].jmag-sig*appmage and model[ug].logg lt 2.0)
	if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,model[ug[u]].propmotion,psym=4,symsize=0.2,color=4

	oploterror,[jkcol[q]],[pm],[cols],[pme],thick=3,errthick=3
endif else !p.multi=[0,1,2]

plot,[0],[0],psym=3,xrange=[-0.2,0.8],yrange=[0.1,1.5],/xs,$
xtitle='H-K',ytitle='J-H'
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
	psym=4,symsize=0.2,color=1	
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
	psym=4,color=3,symsize=0.2
oploterror,[hkcol[q]],[jhcol[q]],[col1s],[col2s],thick=3,errthick=3;,color=6,errcolor=6

plot,[0],[0],psym=3,/xs,xrange=[0,3.0],yrange=[0,1.5],/ys,$
xtitle='g-r',ytitle='r-i'
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
oplot,model[ug[u]].mapp_g-model[ug[u]].mapp_r,model[ug[u]].mapp_r-model[ug[u]].mapp_i,$
psym=4,symsize=0.2,color=1
u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
if (u[0] ne -1) then oplot,model[ug[u]].mapp_g-model[ug[u]].mapp_r,model[ug[u]].mapp_r-model[ug[u]].mapp_i,$
psym=4,color=3,symsize=0.2
oploterror,[grcol[q]],[ricol[q]],[col3s],[col4s],thick=3,errthick=3;,color=6,errcolor=6
!p.multi=0		

an=''
read,an

end


pro plotplx2,data,model,ug,q,redpm,jkcol,hkcol,jhcol,grcol,ricol,cols,col1s,$
col2s,col3s,col4s,sig,appmage,plx,plxe,jhk,bvcol,col5s

!p.multi=0.

wset,1
!p.charsize=2

if (jhk) then begin
    !p.multi=[0,1,3] 

    plot,[0],[0],psym=3,xrange=[-0.2,2.0],yrange=[1000,0.01],/xs,$
    xtitle='J-K',ytitle='Parallax',/ylog		

    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg gt 4.0)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=1

    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg gt 3.5 and model[ug].logg lt 4.0)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=2

    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg gt 2.0 and model[ug].logg lt 3.5)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=3

    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg lt 2.0)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_j-model[ug[u]].mapp_k,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=4
    oploterror,[jkcol[q]],[plx],[cols],[plxe],thick=3,errthick=3

    plot,[0],[0],psym=3,xrange=[-0.2,0.8],yrange=[0.1,1.5],/xs,$
    xtitle='H-K',ytitle='J-H'
    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
    if (u[0] ne -1 and n_elements(u) gt 1) then $
            oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
            psym=4,symsize=0.2,color=1	
    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
    if (u[0] ne -1 and n_elements(u) gt 1) then $
            oplot,model[ug[u]].mapp_h-model[ug[u]].mapp_k,model[ug[u]].mapp_j-model[ug[u]].mapp_h,$
            psym=4,color=3,symsize=0.2
    oploterror,[hkcol[q]],[jhcol[q]],[col1s],[col2s],thick=3,errthick=3;,color=6,errcolor=6

    plot,[0],[0],psym=3,/xs,xrange=[0,3.0],yrange=[0,1.5],/ys,$
    xtitle='g-r',ytitle='r-i'
    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg gt 3.5)			
    oplot,model[ug[u]].mapp_g-model[ug[u]].mapp_r,model[ug[u]].mapp_r-model[ug[u]].mapp_i,$
    psym=4,symsize=0.2,color=1
    u=where(model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt $
    data[q].jmag-sig*appmage and model[ug].logg lt 3.5)		
    if (u[0] ne -1) then oplot,model[ug[u]].mapp_g-model[ug[u]].mapp_r,model[ug[u]].mapp_r-model[ug[u]].mapp_i,$
    psym=4,color=3,symsize=0.2
    oploterror,[grcol[q]],[ricol[q]],[col3s],[col4s],thick=3,errthick=3;,color=6,errcolor=6
    !p.multi=0		

 endif else begin

    !p.multi=0
    plot,[0],[0],psym=3,xrange=[-0.2,2.5],yrange=[1000,0.01],/xs,$
    xtitle='B-V',ytitle='Parallax',/ylog		

    u=where(model[ug].mapp_v lt data[q].vmag+sig*appmage and model[ug].mapp_v gt $
    data[q].vmag-sig*appmage and model[ug].logg gt 4.0)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=1

    u=where(model[ug].mapp_v lt data[q].vmag+sig*appmage and model[ug].mapp_v gt $
    data[q].vmag-sig*appmage and model[ug].logg gt 3.5 and model[ug].logg lt 4.0)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=2

    u=where(model[ug].mapp_v lt data[q].vmag+sig*appmage and model[ug].mapp_v gt $
    data[q].vmag-sig*appmage and model[ug].logg gt 2.0 and model[ug].logg lt 3.5)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=3

    u=where(model[ug].mapp_v lt data[q].vmag+sig*appmage and model[ug].mapp_v gt $
    data[q].vmag-sig*appmage and model[ug].logg lt 2.0)
    if (u[0] ne -1 and n_elements(u) gt 1) then $
    oplot,model[ug[u]].mapp_b-model[ug[u]].mapp_v,1000D/model[ug[u]].dis,psym=4,symsize=0.2,color=4
    oploterror,[bvcol[q]],[plx],[col5s],[plxe],thick=3,errthick=3

endelse

!p.multi=0

an=''
read,an

end
