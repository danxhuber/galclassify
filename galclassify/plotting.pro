; various plotting stuff
pro plotconstraints,model,hkcol,jhcol,hkcole,jhcole,grcol,ricol,grcole,ricole,$
bvcol,ubcol,bvcole,ubcole,jkcol,jkcole,pm,pme,plx,plxe,teff_in,logg_in,$
teffe_in,logge_in

!p.multi=[0,2,3]
!p.charsize=2.5
; JHK
plot,[0],[0],psym=3,xrange=[-0.2,0.8],yrange=[0.1,1.5],/xs,$
xtitle='H-K',ytitle='J-H',title='JHK'
u=where(model.logg gt 3.5)			
if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[u].mapp_h-model[u].mapp_k,model[u].mapp_j-model[u].mapp_h,$
	psym=4,symsize=0.2,color=1	
u=where(model.logg lt 3.5)		
if (u[0] ne -1 and n_elements(u) gt 1) then $
	oplot,model[u].mapp_h-model[u].mapp_k,model[u].mapp_j-model[u].mapp_h,$
	psym=4,color=3,symsize=0.2
oploterror,[hkcol],[jhcol],[hkcole],[jhcole],thick=3,errthick=3;,color=6,errcolor=6


; gri
plot,[0],[0],psym=3,/xs,xrange=[0,3.0],yrange=[0,1.5],/ys,$
xtitle='g-r',ytitle='r-i',title='gri'
u=where(model.logg gt 3.5)			
oplot,model[u].mapp_g-model[u].mapp_r,model[u].mapp_r-model[u].mapp_i,$
psym=4,symsize=0.2,color=1
u=where(model.logg lt 3.5)		
if (u[0] ne -1) then oplot,model[u].mapp_g-model[u].mapp_r,model[u].mapp_r-model[u].mapp_i,$
psym=4,color=3,symsize=0.2
oploterror,[grcol],[ricol],[grcole],[ricole],thick=3,errthick=3;,color=6,errcolor=6


; UBV
plot,[0],[0],psym=3,/xs,xrange=[0,2.0],yrange=[-0.5,2.5],/ys,$
xtitle='B-V',ytitle='U-B',title='UBV'
u=where(model.logg gt 3.5)			
oplot,model[u].mapp_b-model[u].mapp_v,model[u].mapp_u-model[u].mapp_b,$
psym=4,symsize=0.2,color=1
u=where(model.logg lt 3.5)		
if (u[0] ne -1) then oplot,model[u].mapp_b-model[u].mapp_v,model[u].mapp_u-model[u].mapp_b,$
psym=4,color=3,symsize=0.2
oploterror,[bvcol],[ubcol],[bvcole],[ubcole],thick=3,errthick=3;,color=6,errcolor=6

; proper motion
plot,[0],[0],psym=3,xrange=[-0.2,1.2],yrange=[1000,0.01],/xs,$
	xtitle='J-K',ytitle='Proper Motion',/ylog,title='Proper Motion'		
u=where(model.logg gt 4.0)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,model[u].propmotion,psym=4,symsize=0.2,color=1
u=where(model.logg gt 3.5 and model.logg lt 4.0)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,model[u].propmotion,psym=4,symsize=0.2,color=2
u=where(model.logg gt 2.0 and model.logg lt 3.5)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,model[u].propmotion,psym=4,symsize=0.2,color=3
u=where(model.logg lt 2.0)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,model[u].propmotion,psym=4,symsize=0.2,color=4
oploterror,[jkcol],[pm],[jkcole],[pme],thick=3,errthick=3


; proper motion
plot,[0],[0],psym=3,xrange=[-0.2,1.2],yrange=[1000,0.01],/xs,$
	xtitle='J-K',ytitle='Parallax',/ylog,title='Parallax'			
u=where(model.logg gt 4.0)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,1000D/model[u].dis,psym=4,symsize=0.2,color=1
u=where(model.logg gt 3.5 and model.logg lt 4.0)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,1000D/model[u].dis,psym=4,symsize=0.2,color=2
u=where(model.logg gt 2.0 and model.logg lt 3.5)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,1000D/model[u].dis,psym=4,symsize=0.2,color=3
u=where(model.logg lt 2.0)
if (u[0] ne -1 and n_elements(u) gt 1) then $
oplot,model[u].mapp_j-model[u].mapp_k,1000D/model[u].dis,psym=4,symsize=0.2,color=4
oploterror,[jkcol],[plx],[jkcole],[plxe],thick=3,errthick=3


; spectra
plot,[0],[0],psym=3,xrange=[8000,3000],yrange=[6,0],/xs,$
xtitle='Teff',ytitle='logg',title='Spectroscopy'
u=where(model.logg gt 4.0)
oplot,model[u].teff,model[u].logg,psym=4,color=1,symsize=0.2
u=where(model.logg lt 4.0 and model.logg gt 3.5)
oplot,model[u].teff,model[u].logg,psym=4,color=2,symsize=0.2
u=where(model.logg lt 3.5 and model.logg gt 2.0)
oplot,model[u].teff,model[u].logg,psym=4,color=3,symsize=0.2
u=where(model.logg lt 2.0)
oplot,model[u].teff,model[u].logg,psym=4,color=4,symsize=0.2

oploterror,[teff_in],[logg_in],[teffe_in],[logge_in],thick=3,$
errthick=3


end
