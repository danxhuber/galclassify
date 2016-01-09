function extcorr,ebv
;return,1.0
;return, (1-tanh((ebv-0.15)/0.03))*0.5*(1-0.6)+0.6
;return, (1-tanh((ebv-0.15)/0.3))*0.5*(1-0.6)+0.6
;return, (1-tanh((ebv-0.15)/0.075))*0.5*(1-0.6)+0.6
;return, (1-tanh((ebv-0.15)/0.1))*0.5*(1-0.6)+0.6

;return, (1-tanh((ebv-0.15)/0.1))*0.5*(1-0.5)+0.5

; default recommended by Sanjib
return, (1-tanh((ebv-0.15)/0.075))*0.5*(1.0-0.6)+0.6


end
