begin


  linlog = 1
; pressure (mb): bottom [near surface] to top
  p      = (/ 1000.,925.,850.,700.,600.,500., \
               400.,300.,250.,200.,150.,100., \
                70.,50.,30.,20.,10. /) 
  nlvl   = dimsizes(p) 
  psfc   = 1013.
  pbot   = 1100.
  ptop   = 10
; temperatures (K) at each p level
  x      = (/  292.,285.,283.,277.,270.,260., \
               250.,235.,225.,215.,207.,207., \
               213.,220.,225.,228.,230. /) 
 
  ; vint = vibeta (p,x,linlog,psfc,pbot,ptop)  ; units= K-mb
  vint = vibeta (p,x,linlog,psfc,700,0)  ; units= K-mb
  print("vint="+vint)                         ; vint ===> vint=255368  K-mb

  ; vint = vint/(psfc-ptop)                    ; normalize (units=K)
  vint = vint/700                    ; normalize (units=K)
  print("vint="+vint)                        ; vint ===> 254.605 K

; This matches an alternative approach
                                             ; dp = layer thickness
  dp   = dpres_plevel(p(3:), psfc, ptop, 0)    
                                             ; weighted sum/sum_of_layer_thickness
; print(dp)
; exit()
  ; fint = dim_sum_n(x*dp, 0)/dim_sum_n(dp,0)
  ; fint = dim_sum_n(x*dp, 0)
  fint = dim_sum_n(x(3:)*dp, 0)/700

  print(dp)  ; (/50.5, 75, 112.5, 125, 100, 100, 100, 75, 50, 50, 40, 25, 20, 15, 10, 5/)
  ; print("dim_sum_n(dp,0)="+dim_sum_n(dp,0))  ; 1003
  print("fint="+fint)                        ; fint ===> 254.605 K
end 