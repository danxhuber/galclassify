pro makegalaxiasim,name,lowlim,uplim,glat,glon,area,nm,model

; make Galaxia simulation (no oversampling)
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
printf,1,'surveyArea                          '+strtrim(string(area,format='(f12.3)'),2)
printf,1,'fSample                             1'
printf,1,'popID                               -1'
printf,1,'warpFlareOn                         1'
printf,1,'seed                                12'
printf,1,'r_max                               1000'
printf,1,'starType                            0'
printf,1,'photoError                          0'
close,1
spawn,'galaxia -r galaxia/input'

ebf_read,'galaxia/'+name+'.ebf','/',data
nmodel=n_elements(data.rad)
os=nm/nmodel

if (os gt 1.) then begin

    print,'oversampling by',os

    ; make Galaxia simulation (no oversampling)
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
    printf,1,'surveyArea                          '+strtrim(string(area,format='(f12.3)'),2)
    printf,1,'fSample                             '+strtrim(string(os,format='(f12.3)'),2)
    printf,1,'popID                               -1'
    printf,1,'warpFlareOn                         1'
    printf,1,'seed                                12'
    printf,1,'r_max                               1000'
    printf,1,'starType                            0'
    printf,1,'photoError                          0'
    close,1
    spawn,'galaxia -r galaxia/input'

endif

spawn,'galaxia -a --psys=UBV '+'galaxia/'+name+'.ebf'
spawn,'galaxia -a --psys=SDSS '+'galaxia/'+name+'.ebf'
readgalaxia,'galaxia/'+name+'.ebf',name,model

ebf_write,'galaxia/'+name+'_model.ebf','/model',model

end

; read the Galaxia output and calculate quantities needed later
pro readgalaxia,file,name,model

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
mapp_u=data.UBV_U;+5*alog10(data.rad*100.0) 
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

; put in structure
n=n_elements(ra)
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
mapp_u:0D,$
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
model.mapp_u=mapp_u
model.mapp_b=mapp_b
model.mapp_v=mapp_v
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
model.mass=((model.rad*r_sun)^2*10^(model.logg)/gconst)/m_sun
model.rho = alog10(model.mass/model.rad^3D)
model.ebv=ebv
model.ebv_inf=ebv_inf

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
