;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; wrapper to make a Galaxia simulation of a Kepler/K2 field, including a bright end extension
;
; required input:
;   cam=cam 	    	    ... K2 campaign number
;
; optional input:
;   maglim=maglim           ... array containing the lower and upper J-band mag limit
;                           ... for the model; default: [0,17]
;   faintcut=faintcut       ... calculate undersampled model above this limit; default=17
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@abs2app.pro
@extcorr.pro

; main routine to calculate a Galaxia model with bright end extension
pro makegalaxiasim,cam=cam,maglim=maglim,faintcut=faintcut

if not keyword_set(maglim) then begin
	print,'no magnitude limit specified, assuming J<17.0'
        maglim=[0.,17.0]
endif

; sample the model to the specified cut +0.5 mag to account for uncertainties
magoff=0.5

; subsampling fraction for model below faintcut

subsamp=0.2 ; this was the default used for the paper
;subsamp=0.05	; for C11 (galactic center)

; adjust oversampling so we have at least this many stars in the brightest mag bin
nt=5000.

; bright mag limit for oversampled population (i.e. we're assuming
; everything brigher than this has a parallax)
lowlimit = 5.

; constants
gconst = double(6.6726e-8)
m_sun = double(1.9891e33)		; g
r_sun = double(6.9599e10)		; cm
r_jup = 71.398e8/r_sun
;nocor=1

; NB: in this entire subroutine glat=galactic longitude and
; glon=galactic latitude
;kepler field
 glat=76.321
 glon=13.496
 name='kepler'

if keyword_set(cam) then begin
	readcol,'coords.txt',campaign,raco,deco,format='I,D,D',/silent
	u=where(campaign eq cam)
	GLACTC, raco[u[0]], deco[u[0]], 2000., glat, glon, 1, /DEGREE
	name='c'+strtrim(string(cam),2)
        ;+'_'+strtrim(string(maglim[0],format='(f6.1)'),2)+'_'$
	;+strtrim(string(maglim[1],format='(f6.1)'),2)
endif

print,'calculating Galaxia model for:',name
dir=name+'_add'

lowlim=maglim[0]-magoff
if keyword_set(faintcut) then uplim=faintcut+magoff else uplim=maglim[1]+magoff

; make the base Galaxia simulation (no oversampling)
openw,1,'galaxia/input'
printf,1,'outputFile                          '+name
printf,1,'outputDir                           ./galaxia/'
printf,1,'photoSys                            DCMC'
printf,1,'magcolorNames                       J,J-Ks'
printf,1,'appMagLimits[0]                     '+strtrim(string(lowlim,format='(f10.2)'),2)
printf,1,'appMagLimits[1]                     '+strtrim(string(uplim,format='(f10.2)'),2)
printf,1,'absMagLimits[0]                     -1000'
printf,1,'absMagLimits[1]                     1000'
printf,1,'colorLimits[0]                      -1000'
printf,1,'colorLimits[1]                      1000'
printf,1,'geometryOption                      1'
printf,1,'longitude                           '+strtrim(string(glat,format='(f12.3)'),2)
printf,1,'latitude                            '+strtrim(string(glon,format='(f12.3)'),2)
printf,1,'surveyArea                          254.47'
printf,1,'fSample                             1'
printf,1,'popID                               -1'
printf,1,'warpFlareOn                         1'
printf,1,'seed                                12'
printf,1,'r_max                               1000'
printf,1,'starType                            0'
printf,1,'photoError                          0'
close,1
spawn,'galaxia -r galaxia/input'
spawn,'galaxia -a --psys=UBV '+'galaxia/'+name+'.ebf'
spawn,'galaxia -a --psys=SDSS '+'galaxia/'+name+'.ebf'
readgalaxia,'galaxia/'+name+'.ebf',name

spawn,'mkdir galaxia/'+dir

; if faint cut is set, calculate an undersampled faint end model
if keyword_set(faintcut) then begin

   lowlim=faintcut-magoff
   uplim=maglim[1]+magoff

   openw,1,'galaxia/input'
   printf,1,'outputFile                          usamp'
   printf,1,'outputDir                           ./galaxia/'+dir
   printf,1,'photoSys                            DCMC'
   printf,1,'magcolorNames                       J,J-Ks'
   printf,1,'appMagLimits[0]                     '+strtrim(string(lowlim,format='(f10.2)'),2)
   printf,1,'appMagLimits[1]                     '+strtrim(string(uplim,format='(f10.2)'),2)
   printf,1,'absMagLimits[0]                     -1000'
   printf,1,'absMagLimits[1]                     1000'
   printf,1,'colorLimits[0]                      -1000'
   printf,1,'colorLimits[1]                      1000'
   printf,1,'geometryOption                      1'
   printf,1,'longitude                           '+strtrim(string(glat,format='(f12.3)'),2)
   printf,1,'latitude                            '+strtrim(string(glon,format='(f12.3)'),2)
   printf,1,'surveyArea                          254.47'
   printf,1,'fSample                             '+strtrim(string(subsamp,format='(f12.3)'),2)
   printf,1,'popID                               -1'
   printf,1,'warpFlareOn                         1'
   printf,1,'seed                                12'
   printf,1,'r_max                               1000'
   printf,1,'starType                            0'
   printf,1,'photoError                          0'
   close,1
   spawn,'galaxia -r galaxia/input'
   spawn,'galaxia -a --psys=UBV galaxia/'+dir+'/usamp.ebf'
   spawn,'galaxia -a --psys=SDSS galaxia/'+dir+'/usamp.ebf'
   
 endif

; make oversampled Galaxia simulations for the brightest mag bins
ebf_read,'galaxia/'+name+'.ebf','/',data
mapp_j_nored=data.DCMC_J+5*alog10(data.rad*100.0) 

low=(maglim[0]-magoff)>lowlimit
if keyword_set(faintcut) then up=faintcut+magoff else up=maglim[1]+magoff
mags=low+findgen((up-low+0.5)/0.5)*0.5
;mags=(maglim[0]-magoff)+findgen(((maglim[1]-maglim[0])+2.*magoff)/0.5)*0.5

maxbin=0.
maxsample=0.
for q=0.,n_elements(mags)-2 do begin
	u=where(mapp_j_nored ge mags[q] and mapp_j_nored le mags[q+1])
	osample=nt/double(n_elements(u))
        if (osample gt maxsample) then maxsample=osample
	if (osample lt 1.) then continue
	maxbin=mags[q]
endfor

; add a couple of mags at the faint end to account for inflated
; uncertainties
mags=low+findgen(((maxbin-low)+1.0)/0.5)*0.5
maxbin=max(mags)

;for q=0.,n_elements(mags)-2 do begin
;	u=where(mapp_j_nored ge mags[q] and mapp_j_nored le mags[q+1])
;	osample=nt/double(n_elements(u))
;	if (osample lt 1.) then osample = 1.
	;print,q,mags[q],mags[q]+1,n_elements(u),osample
	
	openw,1,'galaxia/input'
	printf,1,'outputFile                          osamp';+strtrim(string(fix(q)),2)
	printf,1,'outputDir                           ./galaxia/'+dir
	printf,1,'photoSys                            DCMC'
	printf,1,'magcolorNames                       J,J-Ks'
	;printf,1,'appMagLimits[0]                     '+strtrim(string(mags[q],format='(f10.2)'),2)
	;printf,1,'appMagLimits[1]                     '+strtrim(string(mags[q+1],format='(f10.2)'),2)
	printf,1,'appMagLimits[0]                     '+strtrim(string(lowlimit,format='(f10.2)'),2)
	printf,1,'appMagLimits[1]                     '+strtrim(string(maxbin,format='(f10.2)'),2)
	printf,1,'absMagLimits[0]                     -1000'
	printf,1,'absMagLimits[1]                     1000'
	printf,1,'colorLimits[0]                      -1000'
	printf,1,'colorLimits[1]                      1000'
	printf,1,'geometryOption                      1'
	printf,1,'longitude                           '+strtrim(string(glat,format='(f12.3)'),2)
	printf,1,'latitude                            '+strtrim(string(glon,format='(f12.3)'),2)
	printf,1,'surveyArea                          254.47'
 	printf,1,'fSample                             '+strtrim(string(fix(maxsample)),2)       
	;printf,1,'fSample                             '+strtrim(string(fix(osample)),2)
	printf,1,'popID                               -1'
	printf,1,'warpFlareOn                         1'
	printf,1,'seed                                12'
	printf,1,'r_max                               1000'
	printf,1,'starType                            0'
	printf,1,'photoError                          0'
	close,1
	spawn,'galaxia -r galaxia/input'
	;spawn,'galaxia -a --psys=UBV galaxia/'+dir+'/osamp'+strtrim(string(fix(q)),2)+'.ebf'
	;spawn,'galaxia -a --psys=SDSS galaxia/'+dir+'/osamp'+strtrim(string(fix(q)),2)+'.ebf'
        spawn,'galaxia -a --psys=UBV galaxia/'+dir+'/osamp.ebf'
	spawn,'galaxia -a --psys=SDSS galaxia/'+dir+'/osamp.ebf'
;endfor

; append all the different subsimulations
spawn,'find galaxia/'+dir+'/ > galaxia/files'
readcol,'galaxia/files',files,format='A'

if (n_elements(files) gt 1) then begin

    files=files[1:n_elements(files)-1]

    for i=0.,n_elements(files)-1 do begin

	    ebf_read,files[i],'/',data

	    abs2app,data,/NOEXT
	    
	    mapp_j=data.DCMC_J;+5*alog10(data.rad*100.0) 
	    mapp_h=data.DCMC_H;+5*alog10(data.rad*100.0) 
	    mapp_k=data.DCMC_Ks;+5*alog10(data.rad*100.0) 
	    mapp_g=data.SDSS_g;+5*alog10(data.rad*100.0) 
	    mapp_r=data.SDSS_r;+5*alog10(data.rad*100.0) 
	    mapp_i=data.SDSS_i;+5*alog10(data.rad*100.0) 
	    mapp_b=data.UBV_B;+5*alog10(data.rad*100.0) 
	    mapp_v=data.UBV_V;+5*alog10(data.rad*100.0) 

	    ; convert to ra and dec
	    GLACTC, ra, dec, 2000D, data.glon, data.glat, 2, /DEGREE

	    propmotion=dblarr(n_elements(ra))
	    for q=0.,n_elements(ra)-1 do propmotion[q] = convtopm(ra[q],dec[q],(data.rad[q])*1000.,$
	    data.vx[q],data.vy[q],data.vz[q]) 

	    par = 1D/(data.rad)
	    radius = sqrt(gconst*data.smass*m_sun/10^data.grav)/r_sun
	    ;numax = 3100.*(10^data.grav/10^4.44 * (10^data.teff/5777.)^(-0.5))
	    logg=data.grav
	    teff=10D^data.teff
	    feh=data.feh
	    rad=radius
	    propmotion=propmotion
	    ra=ra
	    dec=dec
	    popid=data.popid
	    age=data.age
	    ebv=data.exbv_schlegel
	    ebv_inf=data.exbv_schlegel_inf
 
	    if (i eq 0) then begin
		    jt=mapp_j
		    ht=mapp_h
		    kt=mapp_k
		    gtm=mapp_g
		    rt=mapp_r
		    it=mapp_i
		    bt=mapp_b
		    vt=mapp_v
		    propmotiont=propmotion
		    loggt=logg
		    tefft=teff
		    feht=feh
		    radt=rad
		    part=par
		    rat=ra
		    det=dec
		    popidt=popid
		    aget=age
		    ebvt=ebv
		    ebvit=ebv_inf
                    
                    if strmatch(files[i],'*osamp*') then sid=replicate(1,n_elements(mapp_j))
                    if strmatch(files[i],'*usamp*') then sid=replicate(2,n_elements(mapp_j))

	    endif else begin
		    jt=[jt,mapp_j]
		    ht=[ht,mapp_h]
		    kt=[kt,mapp_k]
		    gtm=[gtm,mapp_g]
		    rt=[rt,mapp_r]
		    it=[it,mapp_i]
		    bt=[bt,mapp_b]
		    vt=[vt,mapp_v]
		    propmotiont=[propmotiont,propmotion]
		    loggt=[loggt,logg]
		    tefft=[tefft,teff]
		    feht=[feht,feh]
		    radt=[radt,rad]
		    part=[part,par]
		    rat=[rat,ra]
		    det=[det,dec]
		    popidt=[popidt,popid]
		    aget=[aget,age]
		    ebvt=[ebvt,ebv]
		    ebvit=[ebvit,ebv_inf]

                    if strmatch(files[i],'*osamp*') then sid=[sid,replicate(1,n_elements(mapp_j))]
                    if strmatch(files[i],'*usamp*') then sid=[sid,replicate(2,n_elements(mapp_j))]
		
	    endelse

    endfor

    save,file='galaxia/'+name+'_add.sav',jt,ht,kt,gtm,rt,it,bt,vt,propmotiont,loggt,tefft,$
    feht,radt,rat,det,part,popidt,aget,ebvt,ebvit

    ; now add in the extension
    restore,'galaxia/'+name+'_add.sav'
    restore,'galaxia/'+name+'.sav'

    ; put in structure
    n=n_elements(mapp_j)+n_elements(jt)
    model=replicate({$
    ;mapp_j_nored:0D,$
    mabs_j:0D,$
    mabs_v:0D,$
    mapp_j:0D,$
    mapp_h:0D,$
    mapp_k:0D,$
    mapp_g:0D,$
    mapp_r:0D,$
    mapp_i:0D,$
    mapp_b:0D,$
    mapp_v:0D,$
    propmotion:0D,$
    ;redpm:0D,$
    ebv:0D,$
    ebv_inf:0D,$
    teff:0D,$
    logg:0D,$
    feh:0D,$
    rad:0D,$
    mass:0D,$
    rho:0D,$
    age:0D,$
    dis:0D,$
    ra:0D,$
    dec:0D,$
    popid:0,$
    sid:0},n)

    ;model.mapp_j_nored = [mapp_j_nored,jtnr]
    model.mapp_j=[mapp_j,jt]
    model.mapp_h=[mapp_h,ht]
    model.mapp_k=[mapp_k,kt]
    model.mapp_g=[mapp_g,gtm]
    model.mapp_r=[mapp_r,rt]
    model.mapp_i=[mapp_i,it]
    model.mapp_b=[mapp_b,bt]
    model.mapp_v=[mapp_v,vt]
    model.propmotion=[propmotion,propmotiont]
    model.logg=[logg,loggt]
    model.teff=[teff,tefft]
    model.feh=[feh,feht]
    model.rad=[rad,radt]
    model.ra=[ra,rat]
    model.dec=[dec,det]
    model.dis=1000D/([par,part])
    model.popid=[popid,popidt]
    model.age=[age,aget]

    ; sid indicates which sampling population the model belongs to
    ; 0 = base model, no oversampling
    ; 1 = oversampled bright stars
    ; 2 = undersampled faint stars
    model.sid=[replicate(0,n_elements(mapp_j)),sid]
   
    model.mass=((model.rad*r_sun)^2*10^(model.logg)/gconst)/m_sun
    model.rho = alog10(model.mass/model.rad^3D)
    model.ebv=[ebv,ebvt]
    model.ebv_inf=[ebv_inf,ebvit]

endif else begin

    restore,'galaxia/'+name+'.sav'

    ; put in structure
    n=n_elements(mapp_j)+n_elements(jt)
    model=replicate({$
    ;mapp_j_nored:0D,$
    mabs_j:0D,$
    mabs_v:0D,$
    mapp_j:0D,$
    mapp_h:0D,$
    mapp_k:0D,$
    mapp_g:0D,$
    mapp_r:0D,$
    mapp_i:0D,$
    mapp_b:0D,$
    mapp_v:0D,$
    propmotion:0D,$
    ;redpm:0D,$
    ebv:0D,$
    ebv_inf:0D,$
    teff:0D,$
    logg:0D,$
    feh:0D,$
    rad:0D,$
    mass:0D,$
    rho:0D,$
    age:0D,$
    dis:0D,$
    ra:0D,$
    dec:0D,$
    popid:0,$
    sid:0},n)

    ;model.mapp_j_nored = [mapp_j_nored,jtnr]
    model.mapp_j=mapp_j
    model.mapp_h=mapp_h
    model.mapp_k=mapp_k
    model.mapp_g=mapp_g
    model.mapp_r=mapp_r
    model.mapp_i=mapp_i
    model.mapp_b=[mapp_b,bt]
    model.mapp_v=[mapp_v,vt]
    model.propmotion=propmotion
    model.logg=logg
    model.teff=teff
    model.feh=feh
    model.rad=rad
    model.ra=ra
    model.dec=dec
    model.dis=1000D/(par)
    model.popid=popid
    model.age=age
    model.sid=replicate(0,n_elements(mapp_j))
    model.mass=((model.rad*r_sun)^2*10^(model.logg)/gconst)/m_sun
    model.rho = alog10(model.mass/model.rad^3D)
    model.ebv=ebv
    model.ebv_inf=ebv_inf

endelse

;save,file=name+'_ext.sav',mapp_j_nored,mapp_j,mapp_h,mapp_k,propmotion,logg,$
;teff,feh,rad,ra,dec,par,popid,age,sid,magcut
;save,file='galaxia/'+name+'_ext.sav',model

; write final model with bright end extension
ebf_write,'galaxia/'+name+'_model.ebf','/model',model

; delete temporary files
spawn,'rm galaxia/'+name+'*sav'
spawn,'rm galaxia/'+name+'.ebf'
spawn,'rm -rf galaxia/'+dir
print,'done.'

end

; degrees to radians
function torad,in
return,!DPI*in/180D
end

; given ra, dec, distance and space velocities, calculate proper motions
; based on inverting the equations from Johnson & Soderblom (1987, ApJ)
function convtopm,ra,dec,dis,u,v,w 
cosd = cos(dec/!RADEG)
sind = sin(dec/!RADEG)
cosa = cos(ra/!RADEG)
sina = sin(ra/!RADEG)
input = [u,v,w]
k = 4.74047     ;Equivalent of 1 A.U/yr in km/s   
 A_G = [ [ 0.0548755604, +0.4941094279, -0.8676661490], $ 
         [ 0.8734370902, -0.4448296300, -0.1980763734], $
         [ 0.4838350155,  0.7469822445, +0.4559837762] ]
      
b = dblarr(3,3)
b[0,0] = cosa*cosd
b[1,0] = sina*cosd
b[2,0] = sind
b[0,1] = -sina
b[1,1] = cosa
b[2,1] = 0.
b[0,2] = -cosa*sind
b[1,2] = -sina*sind
b[2,2] = cosd

t=a_g#b
res=invert(t)#input
plx = 1e3/dis
pmra = res[1]*plx/k
pmde = res[2]*plx/k
totpm = sqrt(pmde^2D + (pmra^2D*cos(torad(dec))*cos(torad(dec)))) 
return,totpm
end


; read the Galaxia output and calculate quantities needed later
pro readgalaxia,file,name

gconst = double(6.6726e-8)
m_sun = double(1.9891e33)		; g
r_sun = double(6.9599e10)		; cm
r_jup = 71.398e8/r_sun
ebf_read,file,'/',data

; convert to apparent magnitudes
abs2app,data,/NOEXT

mapp_j=data.DCMC_J;+5*alog10(data.rad*100.0) 
mapp_h=data.DCMC_H;+5*alog10(data.rad*100.0) 
mapp_k=data.DCMC_Ks;+5*alog10(data.rad*100.0) 
mapp_g=data.SDSS_g;+5*alog10(data.rad*100.0) 
mapp_r=data.SDSS_r;+5*alog10(data.rad*100.0) 
mapp_i=data.SDSS_i;+5*alog10(data.rad*100.0) 
mapp_b=data.UBV_B;+5*alog10(data.rad*100.0) 
mapp_v=data.UBV_V;+5*alog10(data.rad*100.0)  

; convert to ra and dec
GLACTC, ra, dec, 2000D, data.glon, data.glat, 2, /DEGREE

; calculate proper motions from space velocities
propmotion=dblarr(n_elements(ra))
for q=0.,n_elements(ra)-1 do propmotion[q] = convtopm(ra[q],dec[q],(data.rad[q])*1000.,$
data.vx[q],data.vy[q],data.vz[q]) 

par = 1D/(data.rad)
radius = sqrt(gconst*data.smass*m_sun/10^data.grav)/r_sun
;numax = 3100.*(10^data.grav/10^4.44 * (10^data.teff/5777.)^(-0.5))
logg=data.grav
teff=10D^data.teff
feh=data.feh
rad=radius
propmotion=propmotion
ra=ra
dec=dec
popid=data.popid
age=data.age
ebv=data.exbv_schlegel
ebv_inf=data.exbv_schlegel_inf

save,file='galaxia/'+name+'.sav',mapp_j,mapp_h,mapp_k,mapp_g,mapp_r,mapp_i,mapp_b,mapp_v,$
propmotion,logg,teff,feh,rad,ra,dec,par,popid,age,ebv,ebv_inf

end

