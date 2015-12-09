;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Classify an arbitrary star given some input observables (colors, spectroscopy, 
;;; proper motions, parallaxes) using a synthetic Galaxia population
;;; This code uses the same method adopted for the K2 EPIC stellar 
;;; classifications, except for calculating a synthetic populations for each
;;; source locally
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@makegalaxiasim.pro
@getpdf.pro
@addextinction.pro
@gsmooth.pro
@plotting.pro
@readinput.pro
@loadcolors.pro

pro galclassify,input=input,sample=sample,pl=pl

if not keyword_set(input) then begin
    print,'no input file specified, stopping.'
    stop
endif

; read input
readinput,input,name,glon,glat,jmag,jmage,hmag,hmage,kmag,kmage,$
umag,umage,bmag,bmage,vmag,vmage,gmag,gmage,rmag,rmage,imag,image,$
pmra,pmrae,pmde,pmdee,plx,plxe,teff,teffe,logg,logge,feh,fehe,irad,nm

; calculate colors and proper motion
bvcol=bmag-vmag
bvcole=sqrt(bmage^2.+vmage^2.)
ubcol=umag-bmag
ubcole=sqrt(umage^2.+bmage^2.)
jhcol=jmag-hmag
jhcole=sqrt(jmage^2.+hmage^2.)
jkcol=jmag-kmag
jkcole=sqrt(jmage^2.+kmage^2.)
hkcol=hmag-kmag
hkcole=sqrt(hmage^2.+kmage^2.)
grcol=gmag-rmag
grcole=sqrt(gmage^2.+rmage^2.)
ricol=rmag-imag
ricole=sqrt(rmage^2.+image^2.)
pm = sqrt(pmde^2D + pmra^2D)
x = pmde^2D + pmra^2D
if (pm ne 0.) then $
pme = sqrt( (x^(-0.5D)*pmde*pmdee)^2D + (x^(-0.5D)*pmra*pmrae)^2D ) else pme = 0D

; calculate synthetic population using Galaxia
lowlim=jmag-0.5
uplim=jmag+0.5
area=!DPI*irad^2.
; don't override old simulations!
if not file_test('galaxia/'+name+'_model.ebf') then $
   makegalaxiasim,name,lowlim,uplim,glat,glon,area,nm,model else $
      ebf_read,'galaxia/'+name+'_model.ebf','/model',model

; output struct (not really needed)
pars = replicate($
	{teffm:0D, teffme1:0D, teffme2:0D, $
	loggm:0D, loggme1:0D, loggme2:0D, $
	fehm:0D, fehme1:0D, fehme2:0D, $
	radm:0D, radme1:0D, radme2:0D, $
	massm:0D, massme1:0D, massme2:0D, $
	rhom:0D, rhome1:0D, rhome2:0D, $
	dism:0D, disme1:0D, disme2:0D, $
	ebvm:0D, ebvme1:0D, ebvme2:0D, $
        teffb:0D, loggb:0D, fehb:0D, radb:0D, massb:0D, rhob:0D, disb:0D, ebvb:0D},1)

; add JHK corrections to model colors
model.mapp_j=model.mapp_j-0.046
model.mapp_h=model.mapp_h-0.020
model.mapp_k=model.mapp_k-0.030

; add extinction to the model colors
mapp_j_nored=model.mapp_j
addextinction,model


;;;; calculate likelihoods; for apparent mag use priority: J>V>g

if (gmag ne 0.) then lh_ap = exp( -(gmag-model.mapp_g)^2D / (2D*gmage^2D) ) 
if (vmag ne 0.) then lh_ap = exp( -(vmag-model.mapp_v)^2D / (2D*vmage^2D) )
if (jmag ne 0.) then lh_ap = exp( -(jmag-model.mapp_j)^2D / (2D*jmage^2D) )

; gri
if (gmag ne 0. and rmag ne 0. and imag ne 0.) then begin
    lh_col3 = exp( -(grcol-(model.mapp_g-model.mapp_r))^2D / (2D*grcole^2D) )   	 
    lh_col4 = exp( -(ricol-(model.mapp_r-model.mapp_i))^2D / (2D*ricole^2D) )   	 
endif else begin
    lh_col3=1.
    lh_col4=1.
endelse

; UBV
if (bmag ne 0. and vmag ne 0. and umag ne 0.) then begin
    lh_col5 = exp( -(bvcol-(model.mapp_b-model.mapp_v))^2D / (2D*bvcole^2D) )
    lh_col6 = exp( -(ubcol-(model.mapp_u-model.mapp_b))^2D / (2D*ubcole^2D) )
endif else begin
    lh_col5=1D
    lh_col6=1D
endelse

; JHK
if (jmag ne 0. and hmag ne 0. and kmag ne 0.) then begin
   lh_col1 = exp( -(jhcol-(model.mapp_j-model.mapp_h))^2D / (2D*jhcole^2D) )
   lh_col2 = exp( -(hkcol-(model.mapp_h-model.mapp_k))^2D / (2D*hkcole^2D) )
endif else begin
   lh_col1=1D
   lh_col2=1D
endelse

; proper motion and parallax
if (pm ne 0.) then lh_pm = exp( -(pm-model.propmotion)^2D / (2D*pme^2D) ) else lh_pm=1D
if (plx ne 0.) then lh_plx = exp( -(plx-(1000D/model.dis))^2D / (2D*plxe^2D) ) else lh_plx=1D

; atmospheric parameters
if (teff ne 0.) then begin
   lh_teff = exp( -(teff-model.teff)^2D / (2D*teffe^2D) )
   lh_logg = exp( -(logg-model.logg)^2D / (2D*logge^2D) )
   lh_feh = exp( -(feh-model.feh)^2D / (2D*fehe^2D) )
endif else begin
   lh_teff=1D
   lh_logg=1D
   lh_feh=1D
endelse

; stepsizes for PDF integration
; teff(%), logg, feh, log(rad), mass(%), log(rho), log(dist)
steps=[0.01,0.005,0.05,0.001,0.01,0.005,0.001,0.001]

if keyword_set(pl) then begin
	wset,1
	plotconstraints,model,hkcol,jhcol,hkcole,jhcole,grcol,ricol,grcole,ricole,$
	bvcol,ubcol,bvcole,ubcole,jkcol,jkcole,pm,pme,plx,plxe,teff,logg,$
	teffe,logge
	wset,0
endif

; joint probability
p=lh_col1*lh_col2*lh_col3*lh_col4*lh_col5*lh_plx*lh_ap*lh_pm*lh_teff*lh_logg*lh_feh
p=p/total(p)
print,'----------------------------'
print,'parameter median +1sig -1sig'
getsum,model,p,pars,0,steps,name,sample=sample,outpath='output/',pl=pl

if keyword_set(sample) then begin
    psamp,model,p,sample,post	    
    openw,1,'output/'+name+'_samp.txt'
    printf,1,'teff,logg,feh,rad,mass,rho,dis'
    for i=0.,n_elements(post)-1 do printf,1,post[i].teff,post[i].logg,post[i].feh,$
    post[i].rad,post[i].mass,post[i].rho,post[i].dis,$
    format='(d8.0,d8.3,d8.3,d15.8,d15.8,d15.8,d15.8)'
close,1
endif

stop

end

