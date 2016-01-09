; this is the same as Sanjib's routine, adapted to deal with IDL sav files that contain
; the bright end mag extensions

@extcorr.pro
;@abs2app.pro
;@aebv_factor.pro

pro addextinction,model

print,'adding extinction ...'

;MAPP_J MAPP_H MAPP_K MAPP_G MAPP_R MAPP_I MAPP_B MAPP_V
;;; GDL can't read this EBF table, so use restore instead
;ebf_read,'aebv_factor.ebf','/',schlegel_table
restore,'aebv_factor.idl'
;stop

ind_j=where(strmatch(schlegel_table.filter,'DCMC_J'))
ind_h=where(strmatch(schlegel_table.filter,'DCMC_H'))
ind_k=where(strmatch(schlegel_table.filter,'DCMC_Ks'))
ind_g=where(strmatch(schlegel_table.filter,'SDSS_g'))
ind_r=where(strmatch(schlegel_table.filter,'SDSS_r'))
ind_i=where(strmatch(schlegel_table.filter,'SDSS_i'))
ind_b=where(strmatch(schlegel_table.filter,'UBV_B'))
ind_v=where(strmatch(schlegel_table.filter,'UBV_V'))

ebv=model.ebv*extcorr(model.ebv_inf)	
	
model.mapp_j=model.mapp_j+(ebv*schlegel_table.aebv_factor[ind_j[0]])
model.mapp_h=model.mapp_h+(ebv*schlegel_table.aebv_factor[ind_h[0]])
model.mapp_k=model.mapp_k+(ebv*schlegel_table.aebv_factor[ind_k[0]])
model.mapp_g=model.mapp_g+(ebv*schlegel_table.aebv_factor[ind_g[0]])
model.mapp_r=model.mapp_r+(ebv*schlegel_table.aebv_factor[ind_r[0]])
model.mapp_i=model.mapp_i+(ebv*schlegel_table.aebv_factor[ind_i[0]])
model.mapp_b=model.mapp_b+(ebv*schlegel_table.aebv_factor[ind_b[0]])
model.mapp_v=model.mapp_v+(ebv*schlegel_table.aebv_factor[ind_v[0]])

; compare this to Sanjib's result
;ebf_read,'galaxia/kepler.ebf','/',test
;abs2app,test,/CORR
;u2=where(test.popid lt 8)
;u=where(model.sid eq 0)
;test1=model[u].mapp_j
;test2=test.DCMC_J[u2]
;stop

end
