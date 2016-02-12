;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; take raw classification output and perform sanity checks
;
; required input:
;   cam=cam 	    	    ... K2 campaign number
;   
; optional input:
;   pl=pl   	    	    ... active plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro consolidate,cam=cam,pl=pl

name='c'+strtrim(string(cam),2);+'_'+strtrim(string(maglim[0],format='(f6.1)'),2)+'_'$
	;+strtrim(string(maglim[1],format='(f6.1)'),2)

; default paths
readcol,'paths.txt',paths,format='A'
pathtocat=paths[0]
pathtoepic=paths[1]
outpath=paths[2]

; load raw classification results
restore,outpath+'k2'+name+'_stparas_raw.sav'

; load Galaxia model
ebf_read,'galaxia/'+name+'_model.ebf','/model',model

u=where(model.mapp_j lt max(pars.jmag)+0.1); and model.popid lt 8)
model=model[u]
pos=where((10.^model.age)*1e-9 lt 0.5 and model.mass lt 0.15,complement=u)
model=model[u]

; kick out targets that weren't classified
u=where(pars.teffm ne 0.)
pars=pars[u]

pars2=pars

; this fixed a bug that existed in classifyk2 after paper acceptance and 
; prevented distance uncertainties from parallaxes to be adopted
u=where(strmatch(pars2.stpropflag,'plx'))
if (u[0] ne -1) then begin
    for q=0.,n_elements(u)-1 do begin
    	if (pars2[u[q]].dise1 ne 0. and pars2[u[q]].dise2 ne 0.) then begin
	    pars2[u[q]].disme1=pars2[u[q]].dise1
    	    pars2[u[q]].disme2=pars2[u[q]].dise2
	endif
    endfor
endif

pars = replicate($
	{epic:0L, stpropflag:'', cor:0, $
	jmag:0D, hmag:0D, kmag:0D, kpmag:0D, $
	gmag:0D, rmag:0D, imag:0D, $
	teff:0D, teffe1:0D, teffe2:0D, $
	logg:0D, logge1:0D, logge2:0D, $
	feh:0D, fehe1:0D, fehe2:0D, $
	rad:0D, rade1:0D, rade2:0D, $
	mass:0D, masse1:0D, masse2:0D, $
	rho:0D, rhoe1:0D, rhoe2:0D, $
	dis:0D, dise1:0D, dise2:0D, $
        ebv:0D, ebve1:0D, ebve2:0D},n_elements(pars))

pars.epic=pars2.epic
pars.stpropflag=pars2.stpropflag
pars.jmag=pars2.jmag
pars.hmag=pars2.hmag
pars.kmag=pars2.kmag
pars.kpmag=pars2.kpmag
;pars.bmag=pars2.bmag
;pars.vmag=pars2.vmag
pars.gmag=pars2.gmag
pars.rmag=pars2.rmag
pars.imag=pars2.imag

ranm=cgrandomindices(long(n_elements(model.teff)),long(1e5))

;;; correct stars that lie in a crazy parameter space because the
;;; median of a bimodal distribution is off
ranm=cgrandomindices(long(n_elements(model.teff)),long(1e5),seed=10)

pos=intarr(n_elements(pars2.teff))
dt=100.
dg=0.1

xax=[6000.,5700.,4800.,3000]
yax=[5.,4.,2.5,2.0]

;goto,skip2
for q=0.,n_elements(pos)-1 do begin
   temp=interpol(yax,xax,pars2[q].teffm)
   if (pars2[q].loggm lt temp or pars2[q].teffm lt 3500.) then continue
   u=where( model[ranm].teff gt (pars2[q].teffm)-dt and model[ranm].teff lt (pars2[q].teffm)+dt $
            and model[ranm].logg gt (pars2[q].loggm)-dg and model[ranm].logg lt (pars2[q].loggm)+dg)
   if u[0] eq -1 then pos[q]=-1
endfor
;skip2:

;pos=match_2d([alog10(pars2.teffm)],[pars2.loggm],alog10(model[ranm].teff),$
;	model[ranm].logg,[0.02,0.02])

if keyword_set(pl) then begin
    !p.charsize=1.5
    !p.multi=[0,1,2]
    plot,model[ranm].teff,model[ranm].logg,psym=3,xrange=[12000,1000],yrange=[6,-1],/xs,/ys,$
    xtitle='teff',ytitle='logg',title='original'
    oplot,xax,yax,color=3
    oplot,pars2.teffm,pars2.loggm,psym=4,color=1
    u=where(pos eq -1 )
    if n_elements(u) gt 2 then oplot,pars2[u].teffm,pars2[u].loggm,psym=4,color=2
 endif

print,'correcting ',n_elements(u),' stars'

print,pars2[u].epic

for q=0.,n_elements(pars)-1 do begin
	
	; first, if the posterior Teff+logg mode are far off the grid, we assume that 
	; the modes correspond to different peaks of the PDFs were chosen. in this case, adopt 
	; the best fit value and scale uncertainties according to posteriors
	if (pos[q] eq -1) then begin
		pars[q].teff=pars2[q].teffb
		pars[q].teffe1=pars2[q].teffbe1;pars2[q].teffb*(pars2[q].teffme1/pars2[q].teffm)
		pars[q].teffe2=pars2[q].teffbe2;pars2[q].teffb*(pars2[q].teffme2/pars2[q].teffm)
		pars[q].logg=pars2[q].loggb
		pars[q].logge1=pars2[q].loggbe1
		pars[q].logge2=pars2[q].loggbe2
		pars[q].feh=pars2[q].fehb
		pars[q].fehe1=pars2[q].fehbe1
		pars[q].fehe2=pars2[q].fehbe2
		pars[q].rad=pars2[q].radb
		pars[q].rade1=pars2[q].radbe1;*(pars2[q].radme1/pars2[q].radm)
		pars[q].rade2=pars2[q].radbe2;*(pars2[q].radme2/pars2[q].radm)
		pars[q].mass=pars2[q].massb
		pars[q].masse1=pars2[q].massbe1;*(pars2[q].massme1/pars2[q].massm)
		pars[q].masse2=pars2[q].massbe2;*(pars2[q].massme2/pars2[q].massm)
		pars[q].rho=pars2[q].rhob
		pars[q].rhoe1=pars2[q].rhobe1;*(pars2[q].rhome1/pars2[q].rhom)
		pars[q].rhoe2=pars2[q].rhobe2;*(pars2[q].rhome2/pars2[q].rhom)
		pars[q].dis=pars2[q].disb
		pars[q].dise1=pars2[q].disbe1;*(pars2[q].disme1/pars2[q].dism)
		pars[q].dise2=pars2[q].disbe2;*(pars2[q].disme2/pars2[q].dism)
                pars[q].ebv=pars2[q].ebvb
		pars[q].ebve1=pars2[q].ebvbe1
		pars[q].ebve2=pars2[q].ebvbe2
	
	endif else begin
	; otherwise, adopt posterior median
		
		pars[q].teff=pars2[q].teffm
		pars[q].teffe1=pars2[q].teffme1
		pars[q].teffe2=pars2[q].teffme2
		pars[q].logg=pars2[q].loggm
		pars[q].logge1=pars2[q].loggme1
		pars[q].logge2=pars2[q].loggme2
		pars[q].feh=pars2[q].fehm
		pars[q].fehe1=pars2[q].fehme1
		pars[q].fehe2=pars2[q].fehme2
		pars[q].dis=pars2[q].dism
		pars[q].dise1=pars2[q].disme1
		pars[q].dise2=pars2[q].disme2
                pars[q].mass=pars2[q].massm
		pars[q].masse1=pars2[q].massme1
		pars[q].masse2=pars2[q].massme2
                pars[q].rad=pars2[q].radm
		pars[q].rade1=pars2[q].radme1
		pars[q].rade2=pars2[q].radme2
                pars[q].rho=pars2[q].rhom
		pars[q].rhoe1=pars2[q].rhome1
		pars[q].rhoe2=pars2[q].rhome2
		pars[q].ebv=pars2[q].ebvm
		pars[q].ebve1=pars2[q].ebvme1
		pars[q].ebve2=pars2[q].ebvme2
	endelse
endfor

if keyword_set(pl) then begin
    plot,model[ranm].teff,model[ranm].logg,psym=3,xrange=[12000,2000],yrange=[6,-1],/xs,/ys,$
    xtitle='teff',ytitle='logg',title='corrected'
    oplot,pars.teff,pars.logg,psym=4,color=1
    u=where(pos eq -1)
    if n_elements(u) gt 2 then oplot,pars[u].teff,pars[u].logg,psym=4,color=2
    !p.multi=0
    an=''
    read,an
endif 


; set a floor on uncertainties since they don't include model systematics
; assuming 1% in Teff + 2% in mass and radius yields 0.02 dex log(g), 6% in density, 
; and 3% in distance 
; NB these are not actually rigorous limits on model systematics (which should be 
; a function of spectral type/ev state), but just broadly avoid crazy small uncertainties
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Teff
tefflimup=10.0
tefflimlow=0.01
u=where(pars.teffe1/pars.teff gt tefflimup)
if (u[0] ne -1) then pars[u].teffe1=pars[u].teff*tefflimup
u=where(pars.teffe2/pars.teff gt tefflimup)
if (u[0] ne -1) then pars[u].teffe2=pars[u].teff*tefflimup

u=where(pars.teffe1/pars.teff lt tefflimlow or pars.teffe1 eq 0.)
if (u[0] ne -1) then pars[u].teffe1=pars[u].teff*tefflimlow
u=where(pars.teffe2/pars.teff lt tefflimlow or pars.teffe2 eq 0.)
if (u[0] ne -1) then pars[u].teffe2=pars[u].teff*tefflimlow

; [Fe/H]
fehlimup=5.0
fehlimlow=0.05
u=where(pars.fehe1 gt fehlimup)
if (u[0] ne -1) then pars[u].fehe1=fehlimup
u=where(pars.fehe2 gt fehlimup)
if (u[0] ne -1) then pars[u].fehe2=fehlimup

u=where(pars.fehe1 lt fehlimlow or pars.fehe1 eq 0.)
if (u[0] ne -1) then pars[u].fehe1=fehlimlow
u=where(pars.fehe2 lt fehlimlow or pars.fehe2 eq 0.)
if (u[0] ne -1) then pars[u].fehe2=fehlimlow

; Mass
masslimup=10.0
masslimlow=0.02
u=where(pars.masse1/pars.mass gt masslimup)
if (u[0] ne -1) then pars[u].masse1=pars[u].mass*masslimup
u=where(pars.masse2/pars.mass gt masslimup)
if (u[0] ne -1) then pars[u].masse2=pars[u].mass*masslimup

u=where(pars.masse1/pars.mass lt masslimlow or pars.masse1 eq 0.)
if (u[0] ne -1) then pars[u].masse1=pars[u].mass*masslimlow
u=where(pars.masse2/pars.mass lt masslimlow or pars.masse2 eq 0.)
if (u[0] ne -1) then pars[u].masse2=pars[u].mass*masslimlow


; for logg, radius, rho and distance we don't set upper bounds,
; since dwarfs can be giants and vice versa
logglimup=100.0
logglimlow=0.02
u=where(pars.logge1 gt logglimup)
if (u[0] ne -1) then pars[u].logge1=logglimup
u=where(pars.logge2 gt logglimup)
if (u[0] ne -1) then pars[u].logge2=logglimup

u=where(pars.logge1 lt logglimlow or pars.logge1 eq 0.)
if (u[0] ne -1) then pars[u].logge1=logglimlow
u=where(pars.logge2 lt logglimlow or pars.logge2 eq 0.)
if (u[0] ne -1) then pars[u].logge2=logglimlow

radlimup=100.0
radlimlow=0.02
u=where(pars.rade1/pars.rad lt radlimlow or pars.rade1 eq 0.)
if (u[0] ne -1) then pars[u].rade1=pars[u].rad*radlimlow
u=where(pars.rade2/pars.rad lt radlimlow or pars.rade2 eq 0.)
if (u[0] ne -1) then pars[u].rade2=pars[u].rad*radlimlow

rholimup=100.0
rholimlow=0.06
u=where(pars.rhoe1/pars.rho lt rholimlow or pars.rhoe1 eq 0.)
if (u[0] ne -1) then pars[u].rhoe1=pars[u].rho*rholimlow
u=where(pars.rhoe2/pars.rho lt rholimlow or pars.rhoe2 eq 0.)
if (u[0] ne -1) then pars[u].rhoe2=pars[u].rho*rholimlow

dislimup=100.0
dislimlow=0.03
u=where(pars.dise1/pars.dis lt dislimlow or pars.dise1 eq 0.)
if (u[0] ne -1) then pars[u].dise1=pars[u].dis*dislimlow
u=where(pars.dise2/pars.dis lt dislimlow or pars.dise2 eq 0.)
if (u[0] ne -1) then pars[u].dise2=pars[u].dis*dislimlow

; set some reasonable floor on all non-zero E(B-V) error bars
u=where(pars.ebv ne 0.)
pars[u].ebve1 = pars[u].ebve1>0.001
pars[u].ebve2 = pars[u].ebve2>0.001
u=where(pars.ebv eq 0.)
pars[u].ebve1 = 0.
pars[u].ebve2 = 0.

;goto,skip
; save the final results
save,file=outpath+'k2'+name+'_stparas.sav',pars

openw,1,outpath+'k2'+name+'_stparas.txt'
printf,1,'# EPIC|Teff|sTeff+|sTeff-|logg|slogg+|slogg-|FeH|sFeH+|sFeH-|Rad|sRad+|sRad-|Mass|sMass+|sMass-|rho|srho+|srho-|dis|dis+|dis-|ebv|ebv+|ebv-|prov'
for q=0.,n_elements(pars)-1 do $
	if (pars[q].teff ne 0.) then $
	printf,1,pars[q].epic,'|',$
	pars[q].teff,'|',pars[q].teffe1,'|',pars[q].teffe2,'|',$
	pars[q].logg,'|',pars[q].logge1,'|',pars[q].logge2,'|',$
	pars[q].feh,'|',pars[q].fehe1,'|',pars[q].fehe2,'|',$
	pars[q].rad,'|',pars[q].rade1,'|',pars[q].rade2,'|',$
	pars[q].mass,'|',pars[q].masse1,'|',pars[q].masse2,'|',$
	pars[q].rho,'|',pars[q].rhoe1,'|',pars[q].rhoe2,'|',$
	pars[q].dis,'|',pars[q].dise1,'|',pars[q].dise2,'|',$
        pars[q].ebv,'|',pars[q].ebve1,'|',pars[q].ebve2,'|',$
	pars[q].stpropflag,$
	format='(I12,A2,I5,A2,I5,A2,I5,A2,d6.3,A2,d6.3,A2,d6.3,A2,d7.3,A2,d6.3,A2,d6.3,A2,d8.3,A2,d8.3,A2,d8.3,A2,d8.3,A2,d8.3,A2,d8.3,A2,e10.3,A2,e10.3,A2,e10.3,A2,e10.3,A2,e10.3,A2,e10.3,A2,d7.4,A2,d7.4,A2,d7.4,A2,A6)'
close,1
;skip:

; do some stats and sanity checks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; target list
readcol,'targetlists/K2Campaign'+strtrim(string(cam),2)+'targets.csv',epicc,dum,dum,dum,propid,$
format='L,D,D,D,A',/silent

print,'------------------------------------------'
print,'campaign',cam
print,'total # of stars:',n_elements(epicc)
print,'total # of stars going into classification:',n_elements(pars)
u=where(pars.teff ne 0.)
print,'total # of stars classified:',n_elements(u)
u=where(pars.teff ne 0. and strmatch(pars.stpropflag,'*plx*'))
print,'total # of stars classified with plx:',n_elements(u)
u=where(pars.teff ne 0. and (strmatch(pars.stpropflag,'*lam*') or $
strmatch(pars.stpropflag,'*apo*') or strmatch(pars.stpropflag,'*seg*') or $
strmatch(pars.stpropflag,'*rav*')))
print,'total # of stars classified with spectra:',n_elements(u)
u=where(pars.teff ne 0. and strmatch(pars.stpropflag,'*rpm*'))
print,'total # of stars classified with rpm:',n_elements(u)
u=where(pars.teff ne 0. and strmatch(pars.stpropflag,'*col*'))
print,'total # of stars classified with col:',n_elements(u)
print,'------------------------------------------'


if keyword_set(pl) then begin

!p.multi=[0,1,2]

plot,pars.teff,pars.logg,psym=4,xrange=[max(pars.teff),3000],yrange=[5.5,0],/xs,/ys,$
xtitle='teff',ytitle='logg',title='corrected',/nodata
u=where(strmatch(pars.stpropflag,'*rpm*'))
if (u[0] ne -1) then oplot,pars[u].teff,pars[u].logg,psym=4
u=where(strmatch(pars.stpropflag,'*col*'))
if (u[0] ne -1) then oplot,pars[u].teff,pars[u].logg,psym=4,color=1
u=where(strmatch(pars.stpropflag,'*apo*'))
if (u[0] ne -1) then oplot,pars[u].teff,pars[u].logg,psym=4,color=2
u=where(strmatch(pars.stpropflag,'*rav*'))
if (u[0] ne -1) then oplot,pars[u].teff,pars[u].logg,psym=4,color=3
u=where(strmatch(pars.stpropflag,'*lam*'))
if (u[0] ne -1) then oplot,pars[u].teff,pars[u].logg,psym=4,color=4
u=where(strmatch(pars.stpropflag,'*plx*'))
if (u[0] ne -1) then oplot,pars[u].teff,pars[u].logg,psym=4,color=5

legend,['rpm','col','apo','rav','lam','plx'],psym=[4,4,4,4,4,4],colors=[255,1,2,3,4,5],/top,/left

plot,pars.teff,pars.logg,xrange=[max(pars.teff),3000],yrange=[6,1],psym=3,/nodata
oploterror,pars.teff,pars.logg,pars.teffe1,pars.logge1,/hibar,color=11,errcolor=11,psym=3
oploterror,pars.teff,pars.logg,pars.teffe2,pars.logge2,/lobar,color=11,errcolor=11,psym=3
oplot,pars.teff,pars.logg,psym=4,thick=2

!p.multi=0

an=''
read,an


; plot sample and highlight different science cases according to
; proposals (red giant galactic archeology, M dwarfs, FGK dwarfs)
match,pars.epic,epicc,ix,iy
parst=pars[ix]
propid=propid[iy]

!p.multi=[0,2,2]
!p.charsize=1.7

plot,[0],[0],psym=4,xrange=[9000,2000],yrange=[6,0],/ys,/xs,xtitle='Teff',ytitle='logg'
u=where(parst.logg ne 0.)
oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=11
if cam eq 1 then u=where(parst.logg ne 0. and strmatch(propid,'*1059*'))
if cam eq 2 then u=where(parst.logg ne 0. and strmatch(propid,'*2051*'))
if cam eq 3 then u=where(parst.logg ne 0. and strmatch(propid,'*3051*'))
if cam eq 4 then u=where(parst.logg ne 0. and strmatch(propid,'*4020*'))
if cam eq 5 then u=where(parst.logg ne 0. and strmatch(propid,'*5020*'))
if cam eq 6 then u=where(parst.logg ne 0. and strmatch(propid,'*6032*'))
if cam eq 7 then u=where(parst.logg ne 0. and strmatch(propid,'*7032*'))
if cam eq 8 then u=where(parst.logg ne 0. and strmatch(propid,'*8042*'))
if (n_elements(u) gt 2) then oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=1
xyouts,8000,1,'GAP',color=1,charsize=1.5

plot,[0],[0],psym=4,xrange=[9000,2000],yrange=[6,0],/ys,/xs
u=where(parst.logg ne 0.)
oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=11
if cam eq 1 then u=where(parst.logg ne 0. and strmatch(propid,'*1053*'))
if cam eq 2 then u=where(parst.logg ne 0. and strmatch(propid,'*2069*'))
if cam eq 3 then u=where(parst.logg ne 0. and strmatch(propid,'*3069*'))
if cam eq 4 then u=where(parst.logg ne 0. and strmatch(propid,'*4097*'))
if cam eq 5 then u=where(parst.logg ne 0. and strmatch(propid,'*5097*'))
if cam eq 6 then u=where(parst.logg ne 0. and strmatch(propid,'*6084*'))
if cam eq 7 then u=where(parst.logg ne 0. and strmatch(propid,'*7087*'))
;; rr lyraes
;if cam eq 7 then u=where(parst.logg ne 0. and strmatch(propid,'*7082*'))
if cam eq 8 then u=where(parst.logg ne 0. and strmatch(propid,'*8056*'))
if (n_elements(u) gt 2) then oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=2
xyouts,8000,1,'M dwarfs planets I',color=2,charsize=1.5

plot,[0],[0],psym=4,xrange=[9000,2000],yrange=[6,0],/ys,/xs
u=where(parst.logg ne 0.)
oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=11
if cam eq 1 then u=where(parst.logg ne 0. and strmatch(propid,'*1054*'))
if cam eq 2 then u=where(parst.logg ne 0. and strmatch(propid,'*2104*'))
if cam eq 3 then u=where(parst.logg ne 0. and strmatch(propid,'*3054*'))
if cam eq 4 then u=where(parst.logg ne 0. and strmatch(propid,'*4033*'))
if cam eq 5 then u=where(parst.logg ne 0. and strmatch(propid,'*5033*'))
if cam eq 6 then u=where(parst.logg ne 0. and strmatch(propid,'*6030*'))
if cam eq 7 then u=where(parst.logg ne 0. and strmatch(propid,'*7030*'))
if cam eq 8 then u=where(parst.logg ne 0. and strmatch(propid,'*8077*'))
if (n_elements(u) gt 2) then oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=3
xyouts,8000,1,'FGK planets',color=3,charsize=1.5
;xyouts,8000,1,'Petigura',color=3,charsize=1.5

plot,[0],[0],psym=4,xrange=[9000,2000],yrange=[6,0],/ys,/xs
u=where(parst.logg ne 0.)
oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=11
if cam eq 1 then u=where(parst.logg ne 0. and strmatch(propid,'*1036*'))
if cam eq 2 then u=where(parst.logg ne 0. and strmatch(propid,'*2107*'))
if cam eq 3 then u=where(parst.logg ne 0. and strmatch(propid,'*3107*'))
if cam eq 4 then u=where(parst.logg ne 0. and strmatch(propid,'*4011*'))
if cam eq 5 then u=where(parst.logg ne 0. and strmatch(propid,'*5011*'))
if cam eq 6 then u=where(parst.logg ne 0. and strmatch(propid,'*6008*'))
if cam eq 7 then u=where(parst.logg ne 0. and strmatch(propid,'*7008*'))
if cam eq 8 then u=where(parst.logg ne 0. and strmatch(propid,'*8068*'))
if (n_elements(u) gt 2) then oplot,parst[u].teff,parst[u].logg,psym=4,symsize=0.4,color=4
xyouts,8000,1,'M dwarf / FGK planets II',color=4,charsize=1.5
!p.multi=0

an=''
read,an

; plot distance histogram
!p.multi=[0,1,2]
plothist,pars.dis,bin=5,xtitle='distance (pc)',ytitle='# of stars',/xlog
plothist,pars.ebv,bin=0.001,xtitle='E(B-V)',ytitle='# of stars'
!p.multi=0

print,'min distance:',min(pars.dis,pos)
print,pars[pos].epic
print,'median distance:',median(pars.dis)
print,'max distance:',max(pars.dis,pos)
print,pars[pos].epic
print,' '
print,'min ebv:',min(pars.ebv)
print,'median ebv:',median(pars.ebv)
print,'max ebv:',max(pars.ebv)

an=''
read,an

; plot relative uncertainty histograms
!p.charsize=2.5
!p.multi=[0,2,8]
plothist,pars.teffe1,bin=20,xrange=[0,1000],title='Teff+'
plothist,pars.teffe2,bin=20,xrange=[0,1000],title='Teff-';,/overplot,color=11
;print,'median Teffe+',median(pars.teffe1/pars.teff)
;print,'median Teffe-',median(pars.teffe2/pars.teff)

plothist,pars.logge1,bin=0.02,title='logg+',xrange=[0,1]
plothist,pars.logge2,bin=0.02,title='logg-',xrange=[0,1];/overplot,color=11
;print,'median logge+',median(pars.logge1)
;print,'median logge-',median(pars.logge2)

plothist,pars.fehe1,bin=0.03,title='feh+'
plothist,pars.fehe2,bin=0.03,title='feh-';/overplot,color=11
;print,'median fehe+',median(pars.fehe1)
;print,'median fehe-',median(pars.fehe2)

plothist,(pars.rade1/pars.rad)<2.,bin=0.02,title='rad+',xrange=[0,1]
plothist,(pars.rade2/pars.rad)<2.,bin=0.02,title='rad-',xrange=[0,1];,/overplot,color=11
;print,'median rade+',median(pars.rade1/pars.rad)
;print,'median rade-',median(pars.rade2/pars.rad)

plothist,pars.masse1/pars.mass,bin=0.02,title='mass+'
plothist,pars.masse2/pars.mass,bin=0.02,title='mass-';/overplot,color=11
;print,'median masse+',median(pars.masse1/pars.mass)
;print,'median masse-',median(pars.masse2/pars.mass)

plothist,(pars.rhoe1/pars.rho)<2.,bin=0.02,title='rho+',xrange=[0,1]
plothist,(pars.rhoe2/pars.rho)<2.,bin=0.02,title='rho-',xrange=[0,1]
;print,'median logrhoe+',median(pars.rhoe1)
;print,'median logrhoe-',median(pars.rhoe2)

plothist,(pars.dise1/pars.dis)<2.,bin=0.02,title='dis+',xrange=[0,1]
plothist,(pars.dise2/pars.dis)<2.,bin=0.02,title='dis-',xrange=[0,1];/overplot,color=11
;print,'median dise+',median(pars.dise1/pars.dis)
;print,'median dise-',median(pars.dise2/pars.dis)

plothist,pars.ebve1,bin=0.001,title='ebv+',xrange=[0,0.1]
plothist,pars.ebve2,bin=0.001,title='ebv-',xrange=[0,0.1];/overplot,color=11
;print,'median dise+',median(pars.dise1/pars.dis)
;print,'median dise-',median(pars.dise2/pars.dis)
!p.multi=0.

endif
stop
end
