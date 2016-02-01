;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; classify K2 targets using a Galaxia model
;
; required input:
;   cam=cam 	    	    ... K2 campaign number
;   pathtocat=pathtocat     ... path to spectroscopic input catalogs
;   pathtoepic=pathtoepic   ... path to EPIC (must be in output format of readepic.pro)
;   outpath=outpath 	    ... path for output
;   
; optional input:
;   sample=sample   	    ... output posteriors
;   pl=pl   	    	    ... enable plotting
;   epicin=epicin   	    ... only do this EPIC ID
;   outpath=outpath 	    ... path for output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@getpdf.pro
@addextinction.pro
@gsmooth.pro
@plotting.pro

pro classifyK2,cam=cam,$
sample=sample,epicin=epicin,$
; plotting
pl=pl,verbose=verbose

if not keyword_set(cam) and not keyword_set(epicin) then begin
    print,'no campaign number or EPIC ID specified, stopping.'
    stop
 endif

; get the campaign ID for a specified EPIC ID
if keyword_set(epicin) then begin
   spawn,'grep -r '+strtrim(string(epicin),2)+' targetlists/*',tmp
   cam=strmid(tmp,22,1)
   cam=long(cam[0])
endif

; default paths
readcol,'paths.txt',paths,format='A'
pathtocat=paths[0]
pathtoepic=paths[1]
outpath=paths[2]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; load model, data, and do various setups
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
; various parameters
minrpmerr=0.2		; minimum uncertainty for reduced proper motion
minpmerr=1.0		; minimum uncertainty for proper motion (mas)
mincolerr=0.03		; minimum uncertainty for colors (mag)

;plxlim=0.8              ; maximum relative error for parallax
sig=5.		        ; sigma limit outside which we don't consider models (to speed things up)
limm_bright=1000.	; minimum number of models within sigma box for a bright star
limm_faint=1000.        ; minimum number of models within sigma box for a faint star
colstep=0.01		; stepsize to increase color uncertainties
rpmstep=0.01		; stepsize to increase reduced proper motion uncertainties
pmstep=0.5		; stepsize to increase proper motion uncertainties
apstep=0.05		; stepsize to increase apparent mag uncertainties

;collim=1.0		; upper limit for color uncertainty
aplim=1.0		; upper limit for apparent mag uncertainty
;appmage=0.1		; default apparent magnitude error (not used)

teffe_in=150.		; input uncertainty for teff
logge_in=0.15		; input uncertainty for logg
fehe_in=0.15		; input uncertainty for [Fe/H]

;magcut=9.0  	    	; magnitude cut for using bright extension

itlim=1000.             ; maximum iterations to find increase uncertainties

range = 3./60./60.  	; position match radius for spectroscopic catalogs (degrees)

; stepsizes for PDF integration
; teff(%), logg, feh, log(rad), mass(%), log(rho), log(dist), E(B-V)
steps=[0.01,0.005,0.05,0.001,0.01,0.005,0.001,0.001]

;ocollim=collim
oaplim=aplim
;oappmage=appmage

; get campaigns coordinates
readcol,'coords.txt',campaign,raco,deco,format='I,D,D',/silent
u=where(campaign eq cam)
cra = raco[u[0]]
cdec = deco[u[0]]

; load EPIC and Galaxia model
name='c'+strtrim(string(cam),2);+'_'+strtrim(string(maglim[0],format='(f6.1)'),2)+'_'$
	;+strtrim(string(maglim[1],format='(f6.1)'),2)
outfile='k2'+name

print,'loading EPIC: ',pathtoepic+'c'+strtrim(string(cam),2)+'epic.ebf'
ebf_read,pathtoepic+'c'+strtrim(string(cam),2)+'epic.ebf','/data',data

print,'loading Galaxia model: ','galaxia/'+name+'_model.ebf'
ebf_read,'galaxia/'+name+'_model.ebf','/model',model

; add JHK corrections to model colors
model.mapp_j=model.mapp_j-0.046
model.mapp_h=model.mapp_h-0.020
model.mapp_k=model.mapp_k-0.030

; add extinction to the model colors
mapp_j_nored=model.mapp_j
addextinction,model

; figure out bright cut
ug=where(model.sid eq 1)
magcut=max(mapp_j_nored[ug])-0.5
;magcut=11.

; calculate reduced proper motions for the model
redpm = model.mapp_j+5.*alog10(model.propmotion/1000.)
model.mabs_j=model.mapp_j-5*alog10((1D/(1000D/model.dis))*100.0) 
model.mabs_v=model.mapp_v-5*alog10((1D/(1000D/model.dis))*100.0) 

; load target list
readcol,'targetlists/K2Campaign'+strtrim(string(cam),2)+'targets.csv',epicc,dum,dum,$
dum,propid,format='L,D,D,D,A',/silent
print,'total # of targets:',n_elements(epicc)

; make various cuts to the EPIC
; first, exclude all sources that aren't stars
u=where(strmatch(string(data.obj),'*STAR*'))
data=data[u]

; then, exclude everything that doesn't have valid 2MASS
; photometry, *unless* it has a parallax and valid BV photometry
; (i.e. it's bright and saturated in 2MASS)
u=where((data.jmag ne 0. and data.kmag ne 0. and data.hmag ne 0.) or $
(data.plx gt 0. and data.bmag gt 0. and data.vmag gt 0.))
data=data[u]

; match targets to EPIC
match,data.epic,epicc,ix,iy
data=data[ix]
epicc=epicc[iy]
propid=propid[iy]
print,'total # of targets that can be classified:',n_elements(epicc)

; drop models that are much fainter than targets
u=where(model.mapp_j lt max(data.jmag)+0.5); and model.popid le 7)
model=model[u]
mapp_j_nored=mapp_j_nored[u]
redpm=redpm[u]

; this next line drops pre-main sequence models for M dwarfs
pos=where((10.^model.age)*1e-9 lt 0.5 and model.mass lt 0.15,complement=u)
model=model[u]
mapp_j_nored=mapp_j_nored[u]
redpm=redpm[u]

; randomly sample model to match apparent mag distribution of data (for plotting only)
count=0
bins=0.2
ix=findgen(n_elements(model))
h=histogram(data.jmag,bin=bins,locations=loc)
for i=0.,n_elements(loc)-1 do begin
	sel=where(model.mapp_j ge loc[i] and model.mapp_j lt loc[i]+bins)
	;print,loc[i],n_elements(sel),h[i]
	if (h[i] eq 0.) then continue
	if (n_elements(sel) gt h[i]) then $
	r = CGRANDOMINDICES(long(n_elements(sel)),long(h[i]),seed=10)  $
	else r = findgen(n_elements(sel))
	if (count eq 0) then ran=ix[sel[r]] else ran=[ran,ix[sel[r]]]
        ;print,loc[i],h[i],n_elements(r)
        ;an=''
        ;read,an
	count++
endfor
mran=model[ran]
;save,file=outpath+outfile+'_random.sav',mran
;stop

; load and match spectroscopic catalogs
print,'loading and matching up spectroscopic catalogs ...'
apogee=mrdfits(pathtocat+'apogee_dr12.fits',1,h,/silent)
distance=sqrt( (apogee.RA-cra)^2D + (apogee.DEC-cdec)^2D )
u=where(distance lt 10. and apogee.TEFF gt 0.)
if (u[0] ne -1) then begin
	apogee=apogee[u] 
	match_apogee=MATCH_2D(data.ra,data.dec,apogee.ra,apogee.dec,range)
endif else match_apogee=replicate(-1,n_elements(epicc))

rave=mrdfits(pathtocat+'rave_dr4.fits',1,h,/silent)
distance=sqrt( (rave._RAJ2000-cra)^2D + (rave._DEJ2000-cdec)^2D )
u=where(distance lt 10. and rave.teffk ne 0.)
if (u[0] ne -1) then begin
	rave=rave[u] 
	match_rave=MATCH_2D(data.ra,data.dec,rave._RAJ2000,rave._DEJ2000,range)
endif else match_rave=replicate(-1,n_elements(epicc))

lamost=mrdfits(pathtocat+'lamost_dr1.fits',1,h,/silent)
distance=sqrt( (lamost.RA-cra)^2D + (lamost.DEC-cdec)^2D )
u=where(distance lt 10. and lamost.teff ne 0.)
if (u[0] ne -1) then begin
	lamost=lamost[u] 
	match_lamost=MATCH_2D(data.ra,data.dec,lamost.RA,lamost.DEC,range)
endif else match_lamost=replicate(-1,n_elements(epicc))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; some basic model + data sanity checks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(verbose) then begin

    ; show spectro catalog coverage in FOV
    ranm=cgrandomindices(long(n_elements(model.teff)),long(1e5))
    plot,model[ranm].ra,model[ranm].dec,psym=3
    oplot,data.ra,data.dec,psym=3,color=1
    oplot,apogee.RA,apogee.DEC,psym=4,color=2
    oplot,rave._RAJ2000,rave._DEJ2000,psym=4,color=3
    oplot,lamost.RA,lamost.DEC,psym=4,color=4
    legend,['APOGEE','RAVE','LAMOST'],colors=[2,3,4],psym=[4,4,4],/bottom,/left
    ;saveimage,'screenshots/'+outfile+'_spec.gif'
    an=''
    read,an

    ; compare the model apparent magnitude distribution of baseline and extended set
    plothist,mapp_j_nored,bin=0.5,/ylog
    u=where(model.sid eq 0.)
    if (u[0] ne -1) then plothist,mapp_j_nored[u],bin=0.5,/ylog,color=1,/overplot
    u=where(model.sid eq 1.)
    if (u[0] ne -1) then plothist,mapp_j_nored[u],bin=0.5,/ylog,color=2,/overplot
    u=where(model.sid eq 2.)
    if (u[0] ne -1) then plothist,mapp_j_nored[u],bin=0.5,/ylog,color=3,/overplot
    oplot,[magcut,magcut],[1e-6,1e10],color=4,thick=3
    oplot,[-10,20],[5000,5000],color=3,thick=3
    ;saveimage,'screenshots/'+outfile+'_mag.gif'
    an=''
    read,an

    ; compare RPM-color distribution of model and data
    !p.charsize=2
    !p.multi=[0,1,2]
    plot,[0],[0],psym=3,xrange=[-0.2,1.3],yrange=[10,-10],/xs,$
    xtitle='J-K',ytitle='RPM_J'
    u=where(model[ran].logg gt 4.0)
    if (u[0] ne -1) then oplot,model[ran[u]].mapp_j-model[ran[u]].mapp_k,redpm[ran[u]],psym=4,color=1,symsize=0.2
    u=where(model[ran].logg gt 3.5 and model[ran].logg lt 4.0)
    if (u[0] ne -1) then oplot,model[ran[u]].mapp_j-model[ran[u]].mapp_k,redpm[ran[u]],psym=4,color=2,symsize=0.2
    u=where(model[ran].logg lt 3.5 and model[ran].logg gt 2.0)
    if (u[0] ne -1) then oplot,model[ran[u]].mapp_j-model[ran[u]].mapp_k,redpm[ran[u]],psym=4,color=3,symsize=0.2
    u=where(model[ran].logg lt 2.0)
    if (u[0] ne -1) then oplot,model[ran[u]].mapp_j-model[ran[u]].mapp_k,redpm[ran[u]],psym=4,color=4,symsize=0.2

    plot,[0],[0],psym=3,xrange=[-0.2,1.3],yrange=[10,-10],/xs,$
    xtitle='J-K',ytitle='RPM_J'
    oplot,data.jmag-data.kmag,data.jmag+5.*alog10(data.pmt/1000.),psym=4,symsize=0.2
	!p.multi=0
    ;saveimage,'screenshots/'+outfile+'_rpm.gif'

    an=''
    read,an

    ; compare color distribution for mag limited sample; these should roughly match if 
    ; reddening is correct
	!p.charsize=2
    !p.multi=[0,2,2,0,1]
    u=where(model.mapp_j lt 13. and model.sid eq 0)
    if (u[0] ne -1) then plothist,model[u].mapp_j-model[u].mapp_k,bin=0.01,xrange=[-0.5,1.2],/xs,$
    xtitle='J-K',/norm1
    u=where(data.jmag lt 13. and data.jmag ne 0.)
    if (u[0] ne -1) then plothist,data[u].jmag-data[u].kmag,bin=0.01,color=1,/overplot,/norm1
    
    u=where(model.mapp_j lt 13. and model.sid eq 0)
    if (u[0] ne -1) then plothist,model[u].mapp_j-model[u].mapp_h,bin=0.01,xrange=[-0.5,1.2],/xs,$
    xtitle='J-H',/norm1
    u=where(data.jmag lt 13. and data.jmag ne 0.)
    if (u[0] ne -1) then plothist,data[u].jmag-data[u].hmag,bin=0.01,color=1,/overplot,/norm1
    
    u=where(model.mapp_g lt 14. and model.sid eq 0)
    if (u[0] ne -1) then plothist,model[u].mapp_g-model[u].mapp_r,bin=0.01,xrange=[-0.5,2.0],/xs,$
    xtitle='g-r',/norm1
    u=where(data.gmag lt 14. and data.gmag ne 0.)
    if (u[0] ne -1) then plothist,data[u].gmag-data[u].rmag,bin=0.01,color=1,/overplot,/norm1
    
    u=where(model.mapp_g lt 14. and model.sid eq 0)
    if (u[0] ne -1) then plothist,model[u].mapp_r-model[u].mapp_i,bin=0.01,xrange=[-0.5,1.0],/xs,$
    xtitle='r-i',/norm1
    u=where(data.gmag lt 14. and data.gmag ne 0. and data.rmag and data.imag)
    if (u[0] ne -1) then plothist,data[u].rmag-data[u].imag,bin=0.01,color=1,/overplot,/norm1
    !p.multi=0.

	an=''
	read,an

endif

!p.charsize=2.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set up arrays and loop over all stars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; results
; b = best fit
; m = posterior mode 
; no suffix = adopted value
pars = replicate($
	{epic:0L, stpropflag:'', $
	jmag:0D, hmag:0D, kmag:0D, kpmag:0D, $
	gmag:0D, rmag:0D, imag:0D, $
	teff:0D, teffe1:0D, teffe2:0D, $
	logg:0D, logge1:0D, logge2:0D, $
	feh:0D, fehe1:0D, fehe2:0D, $
	rad:0D, rade1:0D, rade2:0D, $
	mass:0D, masse1:0D, masse2:0D, $
	rho:0D, rhoe1:0D, rhoe2:0D, $
	dis:0D, dise1:0D, dise2:0D, $
	ebv:0D, ebve1:0D, ebve2:0D, $
	teffm:0D, teffme1:0D, teffme2:0D, $
	loggm:0D, loggme1:0D, loggme2:0D, $
	fehm:0D, fehme1:0D, fehme2:0D, $
	radm:0D, radme1:0D, radme2:0D, $
	massm:0D, massme1:0D, massme2:0D, $
	rhom:0D, rhome1:0D, rhome2:0D, $
	dism:0D, disme1:0D, disme2:0D, $
	ebvm:0D, ebvme1:0D, ebvme2:0D, $
	teffb:0D, loggb:0D, fehb:0D, radb:0D, massb:0D, rhob:0D, disb:0D, ebvb:0D},n_elements(epicc))
	pars.epic=epicc

; the Galaxia model has 2-3 groups: an oversampled bright end, and a standard (base) model, 
; and (for C2 and C7) an undersampled faint end 
; get indices of these subsamples, and precalculate colors
ug=where(model.sid)
mjh1=model[ug].mapp_j-model[ug].mapp_h
mhk1=model[ug].mapp_h-model[ug].mapp_k
map1=model[ug].mapp_j
mrpm1=redpm[ug]
mpm1=model[ug].propmotion

ug=where(model.sid eq 0)
mjh0=model[ug].mapp_j-model[ug].mapp_h
mhk0=model[ug].mapp_h-model[ug].mapp_k
map0=model[ug].mapp_j
mrpm0=redpm[ug]
mpm0=model[ug].propmotion

; check whether the model has an undersampled faint extension
ug=where(model.sid eq 2)
if (ug[0] ne -1) then begin
   faintcut=min(model[ug].mapp_j+0.5)
   mjh2=model[ug].mapp_j-model[ug].mapp_h
   mhk2=model[ug].mapp_h-model[ug].mapp_k
   map2=model[ug].mapp_j
   mrpm2=redpm[ug]
   mpm2=model[ug].propmotion
endif else faintcut=0

; store some stuff in arrays; not actually needed ...
modind=findgen(n_elements(model))
jkcol=dblarr(n_elements(epicc))
jhcol=dblarr(n_elements(epicc))
hkcol=dblarr(n_elements(epicc))
grcol=dblarr(n_elements(epicc))
ricol=dblarr(n_elements(epicc))
bvcol=dblarr(n_elements(epicc))
rpm=dblarr(n_elements(epicc))
rpme=dblarr(n_elements(epicc))

ct=0.
; loop over and classify stars
for q=0.,n_elements(epicc)-1 do begin

    	; reset parameters that lazily change in the loop
        ;collim=ocollim
        aplim=oaplim
        ;appmage=oappmage

        ;if (data[q].gmag eq 0. or data[q].plx eq 0.) then continue
        ;if (data[q].pmte le data[q].pmt) then continue

        ;if (data[q].plx le 0. or data[q].jmag ne 0.) then continue
        ;if (data[q].plx eq 0. or data[q].plxe/data[q].plx lt 0.8) then continue

	if keyword_set(epicin) then begin
	    u=where(epicc[q] eq epicin)
	    if (u[0] eq -1) then continue
	endif

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; print various input about the target
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;		
	print,'-------'	
	print,'EPIC',epicc[q]	
	
        ;if (data[q].jmag gt 6.0) then continue

	; calculate colors and reduced proper motions
	jkcol[q] = data[q].jmag-data[q].kmag
	jhcol[q] = data[q].jmag-data[q].hmag
	grcol[q] = data[q].gmag-data[q].rmag
	ricol[q] = data[q].rmag-data[q].imag
	hkcol[q] = data[q].hmag-data[q].kmag
        bvcol[q] = data[q].bmag-data[q].vmag
	
	if (data[q].pmt gt 0.) then begin
		rpm[q] = data[q].jmag+5.*alog10(data[q].pmt/1000.)  
		rpme[q] = sqrt((5D/((data[q].pmt/1000D)*alog(10.)))^2*(data[q].pmte/1000D)^2+data[q].jmage^2)
	endif else begin
		rpm[q]=0.
		rpme[q]=0.
        endelse

    	; uncertainties of observables
          if (data[q].pmt gt 0.) then rpme[q] = rpme[q]>minrpmerr
          col1s=(sqrt(data[q].jmage^2D + data[q].hmage^2D))>mincolerr
          col2s=(sqrt(data[q].hmage^2D + data[q].kmage^2D))>mincolerr
          col3s=(sqrt(data[q].gmage^2D + data[q].rmage^2D))>mincolerr
          col4s=(sqrt(data[q].rmage^2D + data[q].image^2D))>mincolerr
          col5s=(sqrt(data[q].bmage^2D + data[q].vmage^2D))>mincolerr
          cols=(sqrt(data[q].jmage^2D + data[q].kmage^2D))>mincolerr
          pms=data[q].pmte>minpmerr
          appmage=data[q].jmage

          print,'JHK:',data[q].jmag,data[q].jmage,data[q].hmag,data[q].hmage,$
	  data[q].kmag,data[q].kmage,format='(A5,d8.3,d8.3,d8.3,d8.3,d8.3,d8.3)'
          print,'BV:',data[q].bmag,data[q].bmage,data[q].vmag,data[q].vmage,format='(A5,d8.3,d8.3,d8.3,d8.3)'
          print,'gri:',data[q].gmag,data[q].gmage,data[q].rmag,data[q].rmage,data[q].imag,data[q].image,$
	  format='(A5,d8.3,d8.3,d8.3,d8.3,d8.3,d8.3)'
          print,'pm:',data[q].pmt,data[q].pmte,format='(A5,d10.5,d10.5)'      
          print,'rpm:',rpm[q],rpme[q],format='(A5,d10.5,d10.5)'
          print,'plx:',data[q].plx,data[q].plxe,format='(A5,d10.5,d10.5)'
          ct++

          pars[q].jmag=data[q].jmag
          pars[q].hmag=data[q].hmag
          pars[q].kmag=data[q].kmag
          pars[q].kpmag=data[q].kepmag

          ; set flag to indicate which Galaxia subsample to use, depending on app mag
          samp='base'
          if (data[q].jmag ne 0. and data[q].jmag lt magcut) then samp='bright'
          if (data[q].jmag eq 0. and data[q].vmag lt magcut) then samp='bright'
          if (faintcut gt 0. and data[q].jmag gt faintcut) then samp='faint'

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;  best case: parallax is available
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
          ;if (abs(data[q].plxe/data[q].plx) lt plxlim and data[q].plx gt 0.) then begin
          if (data[q].plx gt 0.) then begin   
             if (data[q].jmag ne 0.) then begin
                ap=data[q].jmag
                ape=data[q].jmage
             endif else begin
                ap=data[q].vmag
                ape=data[q].vmage
             endelse
	     ; distance, absolute magnitude + uncertainty
             absmag = -5D*alog10(1D/(data[q].plx/1000D))+ap+5D
             absmage = sqrt((-5D/(data[q].plx*alog(10)))^2*data[q].plxe^2+ape^2)
             dis = 1D/(data[q].plx/1000D)
             dise = sqrt((-1D/(data[q].plx/1000D)^2)^2D*(data[q].plxe/1000D)^2D)
          endif else absmag = 0.
          
          if (absmag ne 0.) then begin

                   if (strmatch(samp,'base')) then ug=where(model.sid eq 0)
                   if (strmatch(samp,'bright')) then ug=where(model.sid eq 1)
                   if (strmatch(samp,'faint')) then ug=where(model.sid eq 2)
                   print,'using ',samp,' model'

                  print,'using parallax ...'
		  pars[q].stpropflag='plx'

                  ; check whether JHK are valid; if not, use BV
                  if (data[q].jmag ne 0 and data[q].hmag ne 0. and data[q].kmag ne 0.) then $
                     jhk=1 else jhk=0

                  print,'absolute mag:',absmag,absmage
                  print,'distance:',dis,dise

                  ; first, constrain by apparent mag and parallax
                  nmodel=0.
                  nit=0.
                  limmodel=limm_bright
                  while (nmodel lt limmodel) do begin 

                     if (jhk) then begin
                        ut=where( (model[ug].mapp_j-model[ug].mapp_h gt jhcol[q]-sig*col1s and $
                           model[ug].mapp_j-model[ug].mapp_h lt jhcol[q]+sig*col1s) $
                           and (model[ug].mapp_h-model[ug].mapp_k gt hkcol[q]-sig*col2s and $
                                model[ug].mapp_h-model[ug].mapp_k lt hkcol[q]+sig*col2s) $
                           and (model[ug].mapp_j lt data[q].jmag+sig*appmage and $
                                model[ug].mapp_j gt data[q].jmag-sig*appmage) $
                           and (1000D/model[ug].dis lt data[q].plx+sig*data[q].plxe and $
                                1000D/model[ug].dis gt data[q].plx-sig*data[q].plxe))
                     endif else begin
                           ut=where( (model[ug].mapp_b-model[ug].mapp_v gt bvcol[q]-sig*col5s and $
                           model[ug].mapp_b-model[ug].mapp_v lt bvcol[q]+sig*col5s) $
                           and (model[ug].mapp_v lt data[q].vmag+sig*appmage and model[ug].mapp_v gt data[q].vmag-sig*appmage) $
                           and (1000D/model[ug].dis lt data[q].plx+sig*data[q].plxe and $
                                1000D/model[ug].dis gt data[q].plx-sig*data[q].plxe))         
                           ;and (model[ug].mabs_v lt absmag+sig*absmage and model[ug].mabs_v gt absmag-sig*absmage))
                     endelse

                     uam=modind[ug[ut]]
                     nmodel=n_elements(uam)
                     if (nmodel lt limmodel) then appmage+=apstep
                     if (appmage gt aplim and nmodel lt limmodel) then begin
                        appmage=aplim
                        col1s+=colstep		
                        col2s+=colstep
                        col5s+=colstep	
                     endif	
                     
                     if (nit gt itlim or (appmage gt aplim and col1s gt 1.0)) then break
                     nit++
                     
                  endwhile

                  ;;; if there are no models around forget about distance and use absolute mag 
                  ;;; this isn't ideal since the d will be wrong, but if a star is that close we
                  ;;; shouldn't be using population models to get its stellar parameters anyway
                  nodis=0
                  if (ut[0] eq -1 or n_elements(ut) lt limmodel) then begin
                          print,'fitting absolute mag'
                          ut=findgen(n_elements(ug))
                          nodis=1
 
                          col1s=(sqrt(data[q].jmage^2D + data[q].hmage^2D))>mincolerr
                          col2s=(sqrt(data[q].hmage^2D + data[q].kmage^2D))>mincolerr
                          col3s=(sqrt(data[q].gmage^2D + data[q].rmage^2D))>mincolerr
                          col4s=(sqrt(data[q].rmage^2D + data[q].image^2D))>mincolerr
                          col5s=(sqrt(data[q].bmage^2D + data[q].vmage^2D))>mincolerr

                          if (jhk) then begin
                             ut=where( (model[ug].mapp_j-model[ug].mapp_h gt jhcol[q]-sig*col1s and $
                             model[ug].mapp_j-model[ug].mapp_h lt jhcol[q]+sig*col1s) $
                             and (model[ug].mapp_h-model[ug].mapp_k gt hkcol[q]-sig*col2s and $
		             model[ug].mapp_h-model[ug].mapp_k lt hkcol[q]+sig*col2s) $
		             and (model[ug].mabs_j lt absmag+sig*absmage and model[ug].mabs_j gt absmag-sig*absmage))
                          endif else begin
                             ut=where( (model[ug].mapp_b-model[ug].mapp_v gt bvcol[q]-sig*col5s and $
                             model[ug].mapp_b-model[ug].mapp_v lt bvcol[q]+sig*col5s) $
		             and (model[ug].mabs_v lt absmag+sig*absmage and model[ug].mabs_v gt absmag-sig*absmage))
                          endelse
                  endif

                  ; if that still doesn't work, assume the star is screwed and bail
                  if (ut[0] eq -1) then continue

    	    	  ; models to be used
                  uam=modind[ug[ut]]                  
                  print,'using appmage:',appmage
                  print,'using cols:',col1s,col2s,col5s
                  print,n_elements(uam),'models within 5 sigma for parallax'

                  ; likelihoods			
                  if (jhk) then begin
                     print,'using JHK'
                     lh_col1 = exp( -(jhcol[q]-(model[uam].mapp_j-model[uam].mapp_h))^2D / (2D*col1s^2D) )
                     lh_col2 = exp( -(hkcol[q]-(model[uam].mapp_h-model[uam].mapp_k))^2D / (2D*col2s^2D) )
                     lh_am = exp( -(absmag-model[uam].mabs_j)^2D / (2D*absmage^2D) )
                     lh_ap = exp( -(data[q].jmag-model[uam].mapp_j)^2D / (2D*appmage^2D) )	
                  endif else begin
                     print,'using BV'
                     lh_col1 = exp( -(bvcol[q]-(model[uam].mapp_b-model[uam].mapp_v))^2D / (2D*col5s^2D) )
                     lh_col2 = 1D
                     lh_ap = exp( -(data[q].vmag-model[uam].mapp_v)^2D / (2D*appmage^2D) )
                     lh_am = exp( -(absmag-model[uam].mabs_v)^2D / (2D*absmage^2D) )
                  endelse
                  
                  lh_plx = exp( -(data[q].plx-(1000D/model[uam].dis))^2D / (2D*data[q].plxe^2D) )

                  if (data[q].gmag ne 0. and data[q].rmag ne 0. and data[q].imag ne 0.) then begin
                          print,'using gri'
                          lh_col3 = exp( -(grcol[q]-(model[uam].mapp_g-model[uam].mapp_r))^2D / (2D*col3s^2D) )
                          lh_col4 = exp( -(ricol[q]-(model[uam].mapp_r-model[uam].mapp_i))^2D / (2D*col4s^2D) )
                          checkcol,model[uam],grcol[q],ricol[q],col3s,col4s,lh_col3,lh_col4
                  endif else begin
                          print,'not using gri'
                          lh_col3=1.
                          lh_col4=1.
                  endelse
	
    	    	  ; joint likelihood
                  if (nodis eq 1) then p=lh_col1*lh_col2*lh_col3*lh_col4*lh_am else $
                  p=lh_col1*lh_col2*lh_col3*lh_col4*lh_ap*lh_plx
                  p=p/total(p)

    	    	  ; get PDFs
                  getsum,model[uam],p,pars,q,steps,pl=pl,outpath=outpath

                  pars[q].dism=dis
                  pars[q].dise1=dise
                  pars[q].dise2=dise
                  pars[q].disb=dis
                  
                  ; if we don't fit a distance, forget about the reddening constraint
                  if (nodis eq 1) then begin
                     pars[q].ebvm=0.
                     pars[q].ebvme1=0.
                     pars[q].ebvme2=0.
                     pars[q].ebvb=0.
                  endif
    	    	    
		  ; sample posterior
	    	  if keyword_set(sample) then begin
	    	    psamp,model[upm],p,sample,post	    
	    	    openw,1,outpath+'/'+strtrim(string(epicc[q]),2)+'_samp.txt'
    	    	    printf,1,'teff,logg,feh,rad,mass,rho,dis'
    	    	    for i=0.,n_elements(post)-1 do printf,1,post[i].teff,post[i].logg,post[i].feh,$
	    	    post[i].rad,post[i].mass,post[i].rho,post[i].dis,$
	    	    format='(d8.0,d8.3,d8.3,d15.8,d15.8,d15.8,d15.8)'
    	    	    close,1
    	    	  endif

                  if keyword_set(pl) then if (nodis eq 1) then $
                     plotplx,data,model,ug,q,sig,absmag,absmage,hkcol,jhcol,$
                             col1s,col2s,bvcol,jkcol,cols,col5s,appmage,jhk else $
                                plotplx2,data,model,ug,q,redpm,jkcol,hkcol,jhcol,grcol,ricol,$
                                         cols,col1s,col2s,col3s,col4s,sig,appmage,data[q].plx,$
                                         data[q].plxe,jhk,bvcol,col5s
 
                  continue

               endif

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;  2nd best case: spectra are available
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	  ; prioritize APOGEE>RAVE>LAMOST
          if (match_apogee[q] ne -1 or match_rave[q] ne -1 or match_lamost[q] ne -1) then begin

                  if (match_apogee[q] ne -1) then begin		
                          ;DR12
                          teff_in = apogee[match_apogee[q]].TEFF
                          logg_in = apogee[match_apogee[q]].LOGG
                          feh_in = apogee[match_apogee[q]].FE_H	
                          ; if apogee doesn't provide metallicity, assume solar
                          if (feh_in lt -100.) then feh_in=0.0
                          pars[q].stpropflag='apo'
                  endif 
                  if (match_rave[q] ne -1) then begin	
                          teff_in = rave[match_rave[q]].TEFFK
                          logg_in = rave[match_rave[q]].LOGGK
                          feh_in = rave[match_rave[q]].C_M_H_K
                          pars[q].stpropflag='rav'
                  endif
                  if (match_lamost[q] ne -1) then begin	
                          teff_in = lamost[match_lamost[q]].TEFF
                          logg_in = lamost[match_lamost[q]].LOGG
                          feh_in = lamost[match_lamost[q]].FEH
                          pars[q].stpropflag='lam'
                  endif

                  print,'input teff/logg/feh:',teff_in,logg_in,feh_in
                  
                  if (strmatch(samp,'base')) then ug=where(model.sid eq 0)
                  if (strmatch(samp,'bright')) then ug=where(model.sid eq 1)
                  if (strmatch(samp,'faint')) then ug=where(model.sid eq 2)
                  print,'using ',samp,' model'

    	    	  ; select models until minimum required number of models is reached
                  nmodel=0.
                  limmodel=limm_bright
                  while (nmodel lt limmodel) do begin 

                          ut=where( (model[ug].teff gt teff_in-sig*teffe_in and model[ug].teff lt teff_in+sig*teffe_in) $
                          and (model[ug].logg gt logg_in-sig*logge_in and model[ug].logg lt logg_in+sig*logge_in) $
                          and (model[ug].mapp_j lt data[q].jmag+sig*appmage and model[ug].mapp_j gt data[q].jmag-sig*appmage) $
                          and (model[ug].feh gt feh_in-sig*fehe_in and model[ug].feh lt feh_in+sig*fehe_in))

                          if (ut[0] eq -1) then begin
                                  usp=findgen(n_elements(ug))
                                  break
                          endif

                          usp=modind[ug[ut]]

                          nmodel=n_elements(usp)
                          if (nmodel lt limmodel) then appmage+=apstep
                          if (appmage gt 5.0) then break	
                  endwhile

                  print,'using appmage:',appmage
                  print,n_elements(usp),'within 5 sigma for SP ('+pars[q].stpropflag+')'

                  lh_ap = exp( -(data[q].jmag-model[usp].mapp_j)^2D / (2D*appmage^2D) )	
                  lh_teff = exp( -(teff_in-model[usp].teff)^2D / (2D*teffe_in^2D) )
                  lh_logg = exp( -(logg_in-model[usp].logg)^2D / (2D*logge_in^2D) )
                  lh_feh = exp( -(feh_in-model[usp].feh)^2D / (2D*fehe_in^2D) )
                  p=lh_teff*lh_logg*lh_feh*lh_ap
                  p=p/total(p)

                  getsum,model[usp],p,pars,q,steps,pl=pl,outpath=outpath

    	    	  ; sample posterior
	    	  if keyword_set(sample) then begin
	    	    psamp,model[upm],p,sample,post	    
	    	    openw,1,outpath+'/'+strtrim(string(epicc[q]),2)+'_samp.txt'
    	    	    printf,1,'teff,logg,feh,rad,mass,rho,dis'
    	    	    for i=0.,n_elements(post)-1 do printf,1,post[i].teff,post[i].logg,post[i].feh,$
	    	    post[i].rad,post[i].mass,post[i].rho,post[i].dis,$
	    	    format='(d8.0,d8.3,d8.3,d15.8,d15.8,d15.8,d15.8)'
    	    	    close,1
    	    	  endif
	
	    	  ; plot observables
    	    	  if keyword_set(pl) then plotspec,model,teff_in,logg_in,teffe_in,logge_in

                  continue
          endif


          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;  if nothing else is available, we're stuck with RPMs
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
	nmodel=0.
	upm=0.
        if strmatch(samp,'bright') then limmodel=limm_bright else limmodel=limm_faint
        if (strmatch(samp,'base')) then ug=where(model.sid eq 0)
        if (strmatch(samp,'bright')) then ug=where(model.sid eq 1)
        if (strmatch(samp,'faint')) then ug=where(model.sid eq 2)
        print,'using ',samp,' model'

        nit=0.
        while (nmodel lt limmodel) do begin 			

           if (data[q].pmt gt 0.) then begin

                if (strmatch(samp,'base')) then begin
                        ut=where((map0 lt data[q].jmag+sig*appmage and map0 gt data[q].jmag-sig*appmage) and $
                        (mjh0 gt jhcol[q]-sig*col1s and mjh0 lt jhcol[q]+sig*col1s) and $
                        (mhk0 gt hkcol[q]-sig*col2s and mhk0 lt hkcol[q]+sig*col2s) and $
                        ;(mrpm0 lt rpm[q]+sig*rpme[q] and mrpm0 gt rpm[q]-sig*rpme[q]))
                        (mpm0 lt data[q].pmt+sig*pms and mpm0 gt data[q].pmt-sig*pms))                     
                endif
                if (strmatch(samp,'faint')) then begin
                        ut=where((map2 lt data[q].jmag+sig*appmage and map2 gt data[q].jmag-sig*appmage) and $
                        (mjh2 gt jhcol[q]-sig*col1s and mjh2 lt jhcol[q]+sig*col1s) and $
                        (mhk2 gt hkcol[q]-sig*col2s and mhk2 lt hkcol[q]+sig*col2s) and $
                        ;(mrpm2 lt rpm[q]+sig*rpme[q] and mrpm2 gt rpm[q]-sig*rpme[q]))
                        (mpm2 lt data[q].pmt+sig*pms and mpm2 gt data[q].pmt-sig*pms))
                endif
                if (strmatch(samp,'bright')) then begin
                        ut=where((map1 lt data[q].jmag+sig*appmage and map1 gt data[q].jmag-sig*appmage) and $
                        (mjh1 gt jhcol[q]-sig*col1s and mjh1 lt jhcol[q]+sig*col1s) and $
                        (mhk1 gt hkcol[q]-sig*col2s and mhk1 lt hkcol[q]+sig*col2s) and $
                        ;(mrpm1 lt rpm[q]+sig*rpme[q] and mrpm1 gt rpm[q]-sig*rpme[q]))
                        (mpm1 lt data[q].pmt+sig*pms and mpm1 gt data[q].pmt-sig*pms))
                endif

             endif else begin

                if (strmatch(samp,'base')) then begin
                        ut=where((map0 lt data[q].jmag+sig*appmage and map0 gt data[q].jmag-sig*appmage) and $
                        (mjh0 gt jhcol[q]-sig*col1s and mjh0 lt jhcol[q]+sig*col1s) and $
                        (mhk0 gt hkcol[q]-sig*col2s and mhk0 lt hkcol[q]+sig*col2s)) 
                endif
                if (strmatch(samp,'faint')) then begin
                        ut=where((map2 lt data[q].jmag+sig*appmage and map2 gt data[q].jmag-sig*appmage) and $
                        (mjh2 gt jhcol[q]-sig*col1s and mjh2 lt jhcol[q]+sig*col1s) and $
                        (mhk2 gt hkcol[q]-sig*col2s and mhk2 lt hkcol[q]+sig*col2s))
                endif
                if (strmatch(samp,'bright')) then begin
                        ut=where((map1 lt data[q].jmag+sig*appmage and map1 gt data[q].jmag-sig*appmage) and $
                        (mjh1 gt jhcol[q]-sig*col1s and mjh1 lt jhcol[q]+sig*col1s) and $
                        (mhk1 gt hkcol[q]-sig*col2s and mhk1 lt hkcol[q]+sig*col2s))
                endif

             endelse

                if (ut[0] ne -1) then upm=modind[ug[ut]]
        ;	print,appmage
                nmodel=n_elements(upm)
                if (nmodel lt limmodel) then appmage+=apstep
                if (appmage gt aplim and nmodel lt limmodel) then begin
                        appmage=aplim
                        col1s+=colstep		
                        col2s+=colstep
                        col3s+=colstep	
                        col4s+=colstep
                        ;rpme[q]+=rpmstep
                        pms+=pmstep
                endif	

                if (nit gt itlim) then break
                nit++
                ;if (col1s gt collim or col2s gt collim) then break		

        endwhile

        print,'did ',nit,' iterations'
	print,'using appmage:',appmage
	print,'using cols:',col1s,col2s,col3s,col4s
	print,'using pme:',pms
;	print,'using rpme:',rpme[q]
	print,n_elements(upm),' models within 5 sigma for RPM'

	; calculate likelihoods: JHK, apparent mag, red. proper motion, and abs mag
	lh_col1 = exp( -(jhcol[q]-(model[upm].mapp_j-model[upm].mapp_h))^2D / (2D*col1s^2D) )
	lh_col2 = exp( -(hkcol[q]-(model[upm].mapp_h-model[upm].mapp_k))^2D / (2D*col2s^2D) )
	;lh_pm = exp( -(rpm[q]-redpm[upm])^2D / (2D*rpme[q]^2D) )
	lh_ap = exp( -(data[q].jmag-model[upm].mapp_j)^2D / (2D*appmage^2D) )
 	lh_pm = exp( -(data[q].pmt-model[upm].propmotion)^2D / (2D*pms^2D) )
			
	if (data[q].gmag ne 0. and data[q].rmag ne 0. and data[q].imag ne 0.) then begin
                print,'using gri'
		lh_col3 = exp( -(grcol[q]-(model[upm].mapp_g-model[upm].mapp_r))^2D / (2D*col3s^2D) )
		lh_col4 = exp( -(ricol[q]-(model[upm].mapp_r-model[upm].mapp_i))^2D / (2D*col4s^2D) )
                checkcol,model[upm],grcol[q],ricol[q],col3s,col4s,lh_col3,lh_col4
	endif else begin
                print,'not using gri'
		lh_col3=1.
		lh_col4=1.
	endelse

	; joint probability
	if (data[q].pmt gt 0.) then begin
		pars[q].stpropflag='rpm'
		p=lh_col1*lh_col2*lh_col3*lh_col4*lh_pm*lh_ap 
	endif else begin
		p=lh_col1*lh_col2*lh_col3*lh_col4
		pars[q].stpropflag='col'
	endelse
	p=p/total(p)
			        
        getsum,model[upm],p,pars,q,steps,epicc[q],sample=sample,pl=pl,outpath=outpath
    	
	; sample posterior
	if keyword_set(sample) then begin
	    psamp,model[upm],p,sample,post	    
	    openw,1,outpath+'/'+strtrim(string(epicc[q]),2)+'_samp.txt'
    	    printf,1,'teff,logg,feh,rad,mass,rho,dis'
    	    for i=0.,n_elements(post)-1 do printf,1,post[i].teff,post[i].logg,post[i].feh,$
	    post[i].rad,post[i].mass,post[i].rho,post[i].dis,$
	    format='(d8.0,d8.3,d8.3,d15.8,d15.8,d15.8,d15.8)'
    	    close,1
    	endif
	
	; plot observables
	if keyword_set(pl) then $
            plotpm,data,model,ug,q,redpm,jkcol,hkcol,jhcol,grcol,ricol,cols,col1s,$
    	    col2s,col3s,col4s,sig,appmage,data[q].pmt,pms
            ;plotrpm,data,model,ug,q,redpm,jkcol,hkcol,jhcol,grcol,ricol,cols,col1s,$
    	    ;col2s,col3s,col4s,sig,appmage,rpm[q],rpme[q]
	if (pars[q].teffm eq 0.) then stop

endfor

if keyword_set(epicin) then stop

save,file=outpath+outfile+'_stparas_raw.sav',pars

end
