; FFT gauss smooth by Dennis
function gauss_smooth, arr, fwhm, FFT=fft, SILENT=silent
;PURPOSE: Smooth an array of equidistant elements with a moving gauss
;         (instead of a boxcar). Note this program's FFT approach avoids
;         edge effects by gluing a "mirror" of arr on both ends of arr
;         before taking fft. This is however, not the Nurep apporach to
;         deal with edge effects in fft convolution. 
;INPUT  : 1: name of array to be smoothed.
;         2: fwhm of gaussian filter (in units of #elements in array).
;         3: set keyword FFT if smoothing should be done in fouier space
;            (which is much faster).
;OUTPUT : 1: smoothed array.
;CALL BY: e.g. gauss_smooth, mag, 10, magsmooth, /FFT
;
;
;14/11-07 by Dennis Stello 
;****************************************************************************************
ON_ERROR,2

n    = 2*N_ELEMENTS(arr)
sig  = fwhm/SQRT(8*alog(2))

IF keyword_set( FFT ) THEN BEGIN
 IF keyword_set( SILENT ) EQ 0 THEN print,'     ->FFT approach'
 w   = n*GAUSSIAN( FINDGEN(n), [1,n/2.,sig])/TOTAL(GAUSSIAN( FINDGEN(n), [1,n/2.,sig]))
 rev = REVERSE(arr)
 out = SHIFT(FFT(FFT([rev(n/4:n/2-1),arr,rev(0:n/4-1)],-1)*FFT(w,-1),1),LONG(n/2.))
 out = out(n/4:n*3/4)
 if (n/2. mod 2. eq 0.) then out = out[0.:n_elements(out)-2.] else out = out[1.:n_elements(out)-1.]
ENDIF ELSE BEGIN
 IF keyword_set( SILENT ) EQ 0 THEN print,'     ->Weighted mean approach'
 n   = N_ELEMENTS(arr)
 out = arr
 FOR i=0L,n-1 DO BEGIN
  w      = GAUSSIAN( FINDGEN(n), [1,i,sig])
  out(i) = TOTAL(arr*w)/TOTAL(w)
 ENDFOR
ENDELSE

return,out

END
