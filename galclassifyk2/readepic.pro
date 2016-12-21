;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in an EPIC catalog
;
; required input:
;   cam=cam 	    	    ... K2 campaign number
;   pathtoepic=pathtoepic   ... path to EPIC files from MAST
;   
; optional input:
;   magcut=magcut   	    ... upper magnitude limit 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro readepic

file='k2_epic_Sep_2016.csv'

n=long(0)
s=''
openr,1,file
readf,1,s
while not eof(1) do begin
	readf,1,s
	;t=strsplit(s,'|',/extract,/preserve_null)
	;kepmag=double(t[45])
	;ra=double(t[9])
	;stop
	;if (abs(ra-centerra) gt 10.) then continue
	;if (kepmag le magcut) then 
	n++ 
endwhile
close,1

data=replicate({epic:0L,$
hip:bytarr(6),$
tyco:bytarr(11),$
ucac:bytarr(10),$
mass:bytarr(17),$
nomad:bytarr(12),$
sdss:bytarr(20),$
obj:bytarr(8),$
kepflag:bytarr(3),$
bmag:0D,vmag:0D,bmage:0D,vmage:0D,jmag:0D,hmag:0D,kmag:0D,jmage:0D,hmage:0D,kmage:0D,$
gmag:0D,gmage:0D,rmag:0D,rmage:0D,imag:0D,image:0D,$
kepmag:0D,pmra:0D,pmrae:0D,pmde:0D,pmdee:0D,pmt:0D,pmte:0D,ra:0D,dec:0D,plx:0D,plxe:0D},n)

i=long(0)
openr,1,file
readf,1,s
while not eof(1) do begin
	readf,1,s
	t=strsplit(s,'|',/extract,/preserve_null)
	; parameter indices are listed in the readme file

	kepmag=double(t[45])
	ra=double(t[9])
	;if (kepmag gt magcut) then continue
	;if (abs(ra-centerra) gt 10.) then continue

	data[i].hip = byte(t[1])
	data[i].tyco = byte(t[2])
	data[i].ucac = byte(t[3])
	data[i].mass = byte(t[4])
	data[i].sdss = byte(t[5])
	;if (n_elements(t) gt 63) then data[i].nomad = byte(t[63])

	data[i].obj = byte(t[6])
	data[i].kepflag = byte(t[7])
	data[i].epic=long(t[0])
	data[i].ra=double(t[9])
	data[i].dec=double(t[10])
	data[i].kepmag=double(t[45])	
	data[i].gmag=double(t[23])
	data[i].gmage=double(t[24])
	data[i].rmag=double(t[25])
	data[i].rmage=double(t[26])
	data[i].imag=double(t[27])
	data[i].image=double(t[28])
	data[i].bmag=double(t[17])
	data[i].bmage=double(t[18])
	data[i].vmag=double(t[19])
	data[i].vmage=double(t[20])
	data[i].jmag=double(t[31])
	data[i].jmage=double(t[32])
	data[i].hmag=double(t[33])
	data[i].hmage=double(t[34])
	data[i].kmag=double(t[35])
	data[i].kmage=double(t[36])	
	data[i].pmra=double(t[11])
	data[i].pmrae=double(t[12])
	data[i].pmde=double(t[13])
	data[i].pmdee=double(t[14])
	data[i].plx=double(t[15])
	data[i].plxe=double(t[16])
	;stop
	i++
endwhile
close,1

data.pmt = sqrt(data.pmde^2D + data.pmra^2D)
x = data.pmde^2D + data.pmra^2D
data.pmte = sqrt( (x^(-0.5D)*data.pmde*data.pmdee)^2D + (x^(-0.5D)*data.pmra*data.pmrae)^2D )

u=where(data.pmt eq 0.)
data[u].pmte=0.

stop
ebf_write,'epic.ebf','/data',data

;save,file=pathtoepic+name+'epic.sav',data

print,'done.'

end
