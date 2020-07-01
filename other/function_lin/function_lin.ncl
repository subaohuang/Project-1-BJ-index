;;;;********************************************************
;;This script is for  scientific research, it contain the 
;; basic functions for calculating some statistic or math questions   
;;                             CopyRight created by Shuheng Lin  
;; 
;;                                           2020.06.30 
;**********************************************************

;;funtion 1
;this funtion is for testing Linear trend derived from 
;ncl function regCoef_n(x, y, dims_x, dims_y)
undef("rc_test")
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end
;;;;;;;

;;funtion 2
;this funtion is for calculating season trend 
undef("cal_season_trend")
function cal_season_trend(data,season)
  begin
    p_djf_r=month_to_season(data, season)
    dim = dimsizes(p_djf_r)
    p_rc=regCoef_n(ispan(1, dim(0), 1), p_djf_r, 0, 0)
    copy_VarCoords(p_djf_r(1,:,:), p_rc)
    p_rc=(/p_rc*dim(0)/)
  return p_rc
end 
;;;;;;;;;

;;funtion 3
;this funtion is for removing season trend 
;;;nino3_4 is standardize index and remove trend 
undef("remove_enso_signal")
function remove_enso_signal(data,nino3_4)
  begin  
    data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
    data_rc = regCoef_n(nino3_4,data_de, 0, 0)
    copy_VarCoords(data(1,:,:), data_rc)
    data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
    m = data_rc_3d 
    m = (/data_rc@yintercept(0)/)

    data_relate_enso  = data_de*data_rc_3d + m 
    data_new = data 
    data_new = (/data - data_relate_enso/)
  return data_new
end 

