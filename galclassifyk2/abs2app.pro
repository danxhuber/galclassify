; The program converts absolute mag to apparent mag 
; and also applies extinction to Galaxia data 
;-------------------------------------------------------------------
; for Galaxia data read as 
; ebf_read,filename,'/',data 
; abs2app,data,/CORR
;-------------------------------------------------------------------
; CORR- Recommended. will appply a correction factor to Schlegel masp 
; for low latitude stars where Schlegel maps overestimate extinction
; DERED: will simulated dereddeing as people generally do , i.e.,
; apply 3d extinction and then correct it using extinction at infinity  
; NOEXT: if ypou do not want to apply extinction
pro abs2app,data,NOEXT=NOEXT,DERED=DERED,CORR=CORR
ON_ERROR,2
tnames=tag_names(data)
surveys=['SDSS','DCMC','UBV','STROEMGREN','KEPLER','2MASS','DENIS']
;surveys=['SDSS','DCMC','UBV','STROEMGREN']
;surveys=['SDSS_','DCMC_','UBV_','STROEMGREN_']
dmod=5*alog10(100.0*data.rad)
j=0L
if keyword_set(DERED) then ebv=data.exbv_schlegel-data.exbv_schlegel_inf else ebv=data.exbv_schlegel
if keyword_set(CORR) then ebv=ebv*extcorr(data.exbv_schlegel_inf)

for i=0,n_elements(tnames)-1 do begin
ind=where((strsplit(tnames[i],'_',/EXTRACT))[0] eq surveys,count)
;; ind=where(strcmp(tnames[i],strupcase(surveys[0]),strlen(surveys[0])) eq 1,count1)
;; ind=where(strcmp(tnames[i],strupcase(surveys[1]),strlen(surveys[1])) eq 1,count2)
;; ind=where(strcmp(tnames[i],strupcase(surveys[2]),strlen(surveys[2])) eq 1,count3)
;; ind=where(strcmp(tnames[i],strupcase(surveys[3]),strlen(surveys[3])) eq 1,count4)
;if (count1+count2+count3+count4) gt 0 then begin
if count gt 0 then begin
if keyword_set(NOEXT) then ebvs=0.0 else ebvs=aebv_factor(tnames[i])*ebv
;print,aebv_factor_old(tnames[i]),aebv_factor(tnames[i])
data.(i)=data.(i)+dmod+ebvs
j++
endif
endfor
if j eq 0 then message,'No quantity to add extinction'
end
