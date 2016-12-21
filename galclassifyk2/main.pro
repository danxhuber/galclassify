;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main wrapper for K2 stellar classification code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@readepic.pro
@makegalaxiasim.pro
@classifyK2.pro
@consolidate.pro

pro main,cam=cam,maglim=maglim,faintcut=faintcut

if not keyword_set(cam) then begin
    print,'no campaign number specified, stopping.'
    stop
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read in the EPIC from ASCII, and store for later use
;readepic,cam=cam,magcut=magcut
;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; make Galaxia simulation of the chosen campaign
;makegalaxiasim,cam=cam,maglim=maglim,faintcut=faintcut
;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; classify targets
classifyK2,cam=cam,/pl,/verbose;,epicin=212763021,/pl;,/verbose,/sample

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; consolidate the output of classifyK2
consolidate,cam=cam,/pl;,outpath=outpath,/pl

stop

end
