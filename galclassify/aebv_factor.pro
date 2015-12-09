;+
; NAME:
;       AEBV_FACTOR
;
; PURPOSE:
;       This computes the value of factor A/E(B-V) that needs to be 
;       multiplied to redenning E(B-V) to get extinction in a given band.
;
; CALLING SEQUENCE:
;       aebv_factor('Landolt_V')
;       aebv_factor('Sloan_r')
;       Look at band names in the code, space under filter name is
;       replaced by '_'
; OUTPUT:
;       A/E(B-V) factor 
;
; REFERENCE
;       David Schlegel, Princeton, 16 April 1999
;        Doug Finkbeiner, Berkeley
;; Extinction in Different Bandpasses
;; ----------------------------------
;;    Assuming an R_V=3.1 extinction curve, the dust maps should be multiplied
;;    by the value in the final column to determine the extinction in a given
;;    passband.  The standard optical-IR bandpasses are represented by the CTIO
;;    and UKIRT transmission curves.  For further details, see Appendix B of
;;    the text.

;;    Filter name       LamEff  A/A(V)   A/E(B-V
;;    ----------------  ------  -------  -------
;;    Landolt U           3372    1.664    5.434
;;    Landolt B           4404    1.321    4.315
;;    Landolt V           5428    1.015    3.315
;;    Landolt R           6509    0.819    2.673
;;    Landolt I           8090    0.594    1.940
;;    CTIO U              3683    1.521    4.968
;;    CTIO B              4393    1.324    4.325
;;    CTIO V              5519    0.992    3.240
;;    CTIO R              6602    0.807    2.634
;;    CTIO I              8046    0.601    1.962
;;    UKIRT J            12660    0.276    0.902
;;    UKIRT H            16732    0.176    0.576
;;    UKIRT K            22152    0.112    0.367
;;    UKIRT L'           38079    0.047    0.153
;;    Gunn g              5244    1.065    3.476
;;    Gunn r              6707    0.793    2.590
;;    Gunn i              7985    0.610    1.991
;;    Gunn z              9055    0.472    1.540
;;    Spinrad R           6993    0.755    2.467
;;    APM b_J             4690    1.236    4.035
;;    Stromgren u         3502    1.602    5.231
;;    Stromgren b         4676    1.240    4.049
;;    Stromgren v         4127    1.394    4.552
;;    Stromgren beta      4861    1.182    3.858
;;    Stromgren y         5479    1.004    3.277
;;    Sloan u'            3546    1.579    5.155
;;    Sloan g'            4925    1.161    3.793
;;    Sloan r'            6335    0.843    2.751
;;    Sloan i'            7799    0.639    2.086
;;    Sloan z'            9294    0.453    1.479
;;    WFPC2 F300W         3047    1.791    5.849
;;    WFPC2 F450W         4711    1.229    4.015
;;    WFPC2 F555W         5498    0.996    3.252
;;    WFPC2 F606W         6042    0.885    2.889
;;    WFPC2 F702W         7068    0.746    2.435
;;    WFPC2 F814W         8066    0.597    1.948
;;    DSS-II g            4814    1.197    3.907
;;    DSS-II r            6571    0.811    2.649
;;    DSS-II i            8183    0.580    1.893
;-------------------------------------------------------------------------------
; MODIFICATION HISTORY:
;       Written, 2011, Sanjib Sharma
;-

function aebv_factor,name
ebf_read,'aebv_factor.ebf','/',schlegel_table
ind=where(strupcase(schlegel_table.filter) eq strupcase(name),count)
if count ne 1 then message,'Filter not found in Schlegel Table: '+name
return, schlegel_table.aebv_factor[ind[0]]
end

;; function aebv_factor,name
;; ebf_read,'/work1/sharma/GsynthData/Isochrones/filter_lambda_eff.ebf','/',schlegel_table
;; ind=where(strupcase(schlegel_table.filter) eq strupcase(name),count)
;; if count ne 1 then message,'Filter not found in Schlegel Table: '+name
;; ;print,name,' ',schlegel_table.filter[ind[0]]
;; return, alambda_ebv(schlegel_table.lambda_eff[ind[0]])
;; end


