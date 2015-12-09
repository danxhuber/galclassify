; define some colors in IDL/GDL
pro loadcolors
device,decomposed=0
tvlct,0,0,0,0
tvlct,255,255,255,255	
tvlct,255,0,0,1     	; red
tvlct,0,255,0,2     	; green
tvlct,0,0,255,3     	; blue
tvlct,255,255,0,4   	; yellow
tvlct,0,255,255,5   	; cyan
tvlct,255,0,255,6   	; magenta
tvlct,0,200,0,7     	; purple

tvlct,140,140,140,11	; greytones
tvlct,180,180,180,12
tvlct,200,200,200,13
tvlct,255,255,255,20	; black
end
