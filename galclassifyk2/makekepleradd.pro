function torad,in
return,!DPI*in/180D
end

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


pro makekepleradd

gconst = double(6.6726e-8)
m_sun = double(1.9891e33)		; g
r_sun = double(6.9599e10)		; cm
r_jup = 71.398e8/r_sun
;nocor=1

readcol,'files',files,format='A'

for i=0.,n_elements(files)-1 do begin

	ebf_read,files[i],'/',data

	abs2app,data,/CORR
	mapp_j=data.DCMC_J;+5*alog10(data.rad*100.0) 
	mapp_k=data.DCMC_Ks;+5*alog10(data.rad*100.0) 
	mapp_h=data.DCMC_H;+5*alog10(data.rad*100.0) 
	
	if (min(mapp_j) gt 8.) then continue

	; convert to ra and dec
	GLACTC, ra, dec, 2000D, data.glon, data.glat, 2, /DEGREE

	propmotion=dblarr(n_elements(ra))
	for q=0.,n_elements(ra)-1 do propmotion[q] = convtopm(ra[q],dec[q],(data.rad[q])*1000.,$
	data.vx[q],data.vy[q],data.vz[q]) 
	dis = (data.rad)*1000.
	par = 1D/(data.rad)
	radius = sqrt(gconst*data.smass*m_sun/10^data.grav)/r_sun
	numax = 3100.*(10^data.grav/10^4.44 * (10^data.teff/5777.)^(-0.5))
	redpm = mapp_j+5.*alog10(propmotion/1000.)

	logg=data.grav
	teff=10^data.teff
	feh=data.feh
	
	u=where(data.popid lt 8)
	mapp_j=mapp_j[u]
	mapp_h=mapp_h[u]
	mapp_k=mapp_k[u]
	redpm=redpm[u]
	logg=logg[u]
	teff=teff[u]
	feh=feh[u]
	rad=radius[u]
	par=par[u]
	propmotion=propmotion[u]
	ra=ra[u]
	dec=dec[u]
	
	if (i eq 0) then begin
		jt=mapp_j
		ht=mapp_h
		kt=mapp_k
		propmotiont=propmotion
		loggt=logg
		tefft=teff
		feht=feh
		radt=rad
		part=par
		rat=ra
		dect=dec
	endif else begin
		jt=[jt,mapp_j]
		ht=[ht,mapp_h]
		kt=[kt,mapp_k]
		propmotiont=[propmotiont,propmotion]
		loggt=[loggt,logg]
		tefft=[tefft,teff]
		feht=[feht,feh]
		radt=[radt,rad]
		part=[part,par]
		rat=[rat,ra]
		dect=[dect,dec]
	endelse

endfor

save,file='keplergalaxia_add.sav',jt,ht,kt,propmotiont,loggt,tefft,$
feht,radt,rat,dect,part

stop
end