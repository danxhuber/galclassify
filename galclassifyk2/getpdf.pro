;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; various subroutines, mostly related to posteriors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; bin up PDF
pro binpdf,x,y,step,xax,yax

xax = min(x) + findgen(1.+((max(x)-min(x)))/step)*step
yax = dblarr(n_elements(xax))
for r=0.,n_elements(xax)-2 do begin
	u = where(x gt xax[r] and x lt xax[r+1])
	if (u[0] eq -1) then continue
	yax[r] = total(y[u])
endfor
xax=xax+step/2.
yax=yax/total(yax)
yax=gauss_smooth(yax,3,/silent)
yax=yax/total(yax)
end


; get mode, medium and confidence interval 
pro getstat,xax,yax,res,err1,err2,cprob
cprob = total(yax,/cumul)
yax2=yax

; mode
;dum=max(yax2,pos)
;res=xax[pos]

; median
dum = min(abs(cprob-0.5),pos)
res=xax[pos]

dum = min(abs(xax-res),pos)
temp = cprob[pos]
ll = temp-temp*0.683D
ul = temp+(1D - temp)*0.683D	
dum = min(abs(cprob-ll),pos5)
err1 = xax[pos]-xax[pos5]
dum = min(abs(cprob-ul),pos6)
err2 = xax[pos6]-xax[pos]
yax=yax2
end


; take probabilities and calculate the PDF and CDF
pro getpdf,xval,pval,res,err1,err2,xax,yax,pl=pl,log=log,no=no,fixed=fixed,leg,truth=truth,samp=samp

x=xval
y=pval

if strmatch(leg,'*dis*') or strmatch(leg,'*Rad*') then x=alog10(x)

; first, take fractional or absolute amount specified
if keyword_set(no) then begin
	dum=max(pval,pos)
	step = xval[pos]*no
endif else step=fixed

ostep=step
binpdf,x,y,step,xax,yax
getstat,xax,yax,res,err1,err2,cprob

; set stepsize to x bins within 1-sigma
sbin=5.0
if (err1 ne 0. and err2 ne 0.) then begin
   if (err1 lt err2) then step=err1/sbin else step=err2/sbin           
   binpdf,x,y,step,xax,yax
   getstat,xax,yax,res,err1,err2,cprob
endif

; if the mode is the first or last bin, assume something is wrong and
; use smaller subsets
try=4
ct=0
while (res eq min(xax) or res eq max(xax)) do begin
   step=step/2.
   binpdf,x,y,step,xax,yax
   getstat,xax,yax,res,err1,err2,cprob
   ct++
   if (ct eq try) then break
endwhile

binpdf,x,y,step,xax,yax
if strmatch(leg,'*dis*') or strmatch(leg,'*rho*') or strmatch(leg,'*Rad*') then begin
    x=10.^(x)
    xax=10.^(xax)
endif    
getstat,xax,yax,res,err1,err2,cprob

; sample the posterior
;y2=(y/max(y))*100.
;numsamp=samp
;samp=[0.]
;n=n_elements(y2)
;while (n_elements(samp) lt numsamp+1.) do begin	
;	ran=randomu(seed)*n
;	if (y2[ran] lt 1.) then continue
;	stop
;	samp=[samp,replicate(x[ran],y2[ran])]
;endwhile
;samp=samp[1:n_elements(samp)-1]
;n=n_elements(xax)
;if strmatch(leg,'*dis*') or strmatch(leg,'*rho*') or strmatch(leg,'*Rad*') then samp=alog10(samp)
;plothist,samp,xsamp,ysamp,bin=step,/noplot    
;ysamp=double(ysamp)
;if strmatch(leg,'*dis*') or strmatch(leg,'*rho*') or strmatch(leg,'*Rad*') then xsamp=10.^xsamp
;ysamp=ysamp/max(ysamp)
;yax=yax/max(yax)

if keyword_set(pl) then begin
;	!p.multi=[0,1,2]
	if keyword_set(log) then plot,xax,cprob,/xlog,title=leg else plot,xax,cprob,title=leg
	oplot,[res,res],[0,1e10],color=1
	oplot,[res-err1,res-err1],[0,1e10],color=1,linestyle=2
	oplot,[res+err2,res+err2],[0,1e10],color=1,linestyle=2
        if keyword_set(truth) then oplot,[truth,truth],[0,1e10],color=2,linestyle=3

	if keyword_set(log) then plot,xax,yax,/xlog,title=leg else $
           plot,xax,yax,title=leg
	;oplot,xsamp,ysamp,color=2
	oplot,[res,res],[0,1e10],color=1
	oplot,[res-err1,res-err1],[0,1e10],color=1,linestyle=2
	oplot,[res+err2,res+err2],[0,1e10],color=1,linestyle=2
        if keyword_set(truth) then oplot,[truth,truth],[0,1e10],color=2,linestyle=3      
	;!P.multi=0
endif

tmp1=err1
tmp2=err2
err1=tmp2
err2=tmp1

print,leg,res,err1,err2,ostep,step

end	


; integrate PDFs for all parameters
pro getsum,data,p,pars,q,steps,name,sample=sample,pl=pl,outpath=outpath

dum=max(p,bfpos)

if keyword_set(pl) then begin
    wset,0
    !p.multi=[0,2,8,0,0]
    !p.charsize=1.75
endif

getpdf,data.teff,p,res,err1,err2,xax,yax,pl=pl,no=steps[0],'Teff'
pars[q].teffm=res
pars[q].teffme1=err1
pars[q].teffme2=err2
pars[q].teffb=data[bfpos].teff
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_teff.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.logg,p,res,err1,err2,xax,yax,pl=pl,fixed=steps[1],'log(g)'
pars[q].loggm=res
pars[q].loggme1=err1
pars[q].loggme2=err2
pars[q].loggb=data[bfpos].logg
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_logg.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.feh,p,res,err1,err2,xax,yax,pl=pl,fixed=steps[2],'[Fe/H]'
pars[q].fehm=res
pars[q].fehme1=err1
pars[q].fehme2=err2
pars[q].fehb=data[bfpos].feh
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_feh.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.rad,p,res,err1,err2,xax,yax,pl=pl,fixed=steps[3],'Rad',/log
pars[q].radm=res
pars[q].radme1=err1
pars[q].radme2=err2
pars[q].radb=data[bfpos].rad
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_rad.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.mass,p,res,err1,err2,xax,yax,pl=pl,no=steps[4],'mass'
pars[q].massm=res
pars[q].massme1=err1
pars[q].massme2=err2
pars[q].massb=data[bfpos].mass
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_mass.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.rho,p,res,err1,err2,xax,yax,pl=pl,fixed=steps[5],'rho',/log
pars[q].rhom=res
pars[q].rhome1=err1
pars[q].rhome2=err2
pars[q].rhob=10.^data[bfpos].rho
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_rho.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.dis,p,res,err1,err2,xax,yax,pl=pl,fixed=steps[6],'distance',/log
pars[q].dism=res
pars[q].disme1=err1
pars[q].disme2=err2
pars[q].disb=data[bfpos].dis
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_dis.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

getpdf,data.ebv,p,res,err1,err2,xax,yax,pl=pl,fixed=steps[7],'E(B-V)',/log
pars[q].ebvm=res
pars[q].ebvme1=err1
pars[q].ebvme2=err2
pars[q].ebvb=data[bfpos].ebv
if keyword_set(sample) then begin
    openw,1,outpath+strtrim(string(name),2)+'_ebv.txt'
    for i=0.,n_elements(xax)-1 do printf,1,xax[i],yax[i]
    close,1
endif

end


; check for bogus gri colors
;pro checkcol,lh1,lh2,lh3,lh4
pro checkcol,model,col1,col2,col1s,col2s,lh3,lh4

;goto,skip2
u1=where(model.mapp_g-model.mapp_r lt col1+5.*col1s and model.mapp_g-model.mapp_r gt col1-5.*col1s)
u2=where(model.mapp_r-model.mapp_i lt col2+5.*col2s and model.mapp_r-model.mapp_i gt col2-5.*col2s)
nmodel=100.
if n_elements(u1) lt nmodel then begin
   print,'dropping g-r'
   lh3[*]=1D
endif

if n_elements(u2) lt nmodel then begin
   print,'dropping r-i'
   lh4[*]=1D
endif
;skip2:

;if (float(n_elements(where(lh1 gt 0.01)))/float(n_elements(lh1)) lt 0.1) then begin
;if (float(n_elements(where(lh2 gt 0.01))) lt nmodel) then begin
;	lh1[*]=1D
;	print,'dropping J-H'
;endif
	
;if (float(n_elements(where(lh2 gt 0.01)))/float(n_elements(lh2)) lt 0.1) then begin
;if (float(n_elements(where(lh2 gt 0.01))) lt nmodel) then begin
;	lh2[*]=1D
;	print,'dropping H-K'
;endif

goto,skip
if (float(n_elements(where(lh3 gt 0.01)))/float(n_elements(lh3)) lt 0.1) then begin
;(float(n_elements(where(lh4 gt 0.01)))/float(n_elements(lh4)) lt 0.1) then begin
;if (float(n_elements(where(lh3 gt 0.01))) lt nmodel) then begin
	lh3[*]=1D
;        lh4[*]=1D
	print,'dropping g-r'
endif

if (float(n_elements(where(lh4 gt 0.01)))/float(n_elements(lh4)) lt 0.1) then begin
;if (float(n_elements(where(lh4 gt 0.01))) lt nmodel) then begin
	lh4[*]=1D
	print,'dropping r-i'
endif
skip:
end


pro psamp,x,y,numsamp,post

post = replicate($
	{teff:0D, logg:0D, feh:0D, rad:0D, $
	mass:0D, rho:0D, dis:0D},numsamp)

; sample the posterior
y2=(y/max(y))*100.
tsamp=[0.]
gsamp=[0.]
rsamp=[0.]
msamp=[0.]
rhosamp=[0.]
dsamp=[0.]
fsamp=[0.]
n=n_elements(y2)
while (n_elements(tsamp) lt numsamp+1.) do begin	
       ran=randomu(seed)*n
       if (y2[ran] lt 1.) then continue
       tsamp=[tsamp,replicate(x[ran].teff,y2[ran])]
	gsamp=[gsamp,replicate(x[ran].logg,y2[ran])]
	rsamp=[rsamp,replicate(x[ran].rad,y2[ran])]
	msamp=[msamp,replicate(x[ran].mass,y2[ran])]
	fsamp=[fsamp,replicate(x[ran].feh,y2[ran])]
	rhosamp=[rhosamp,replicate(x[ran].rho,y2[ran])]
	dsamp=[dsamp,replicate(x[ran].dis,y2[ran])]
endwhile
	
if n_elements(tsamp) gt numsamp then begin
	post.teff=tsamp[1:numsamp]
	post.logg=gsamp[1:numsamp]
	post.rad=rsamp[1:numsamp]
	post.mass=msamp[1:numsamp]
	post.feh=fsamp[1:numsamp]
	post.rho=rhosamp[1:numsamp]
	post.dis=dsamp[1:numsamp]
endif else begin
	post.teff=tsamp[1:numsamp-1]
	post.logg=gsamp[1:numsamp-1]
	post.rad=rsamp[1:numsamp-1]
	post.mass=msamp[1:numsamp-1]
	post.feh=fsamp[1:numsamp-1]
	post.rho=rhosamp[1:numsamp-1]
	post.dis=dsamp[1:numsamp-1]
endelse

end
