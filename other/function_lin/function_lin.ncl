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
  ; ;;;读入HadiSST海温数据
  ;     diri="/home/ys17-19/lsh/data/sst/"
  ;     fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  ;     f1=addfile(fils1, "r")
  ;     date:=cd_calendar(f1->time, 1)

  ;     timeind:=ind(date.le.201512.and.date.ge.197901)

  ;     sst =  lonFlip(f1->sst(timeind,:,:))
   
  ;     ssta=(/rmMonAnnCycTLL(sst)/)
  ;     copy_VarCoords(sst, ssta)

  ;     ssta_de = ssta
  ;     ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
  ;     ssta_index = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
  ;     printVarSummary(ssta_index)

  ;     nino3_4 = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
  ;     printVarSummary(nino3_4)
      
  ;     ; nino3_4_new = nino3_4                ;;;;前期ENSO
  ;     ; nino3_4_new(0:33) = (/nino3_4(1:)/)
  ;     ; nino3_4_new(34) =ssta_index(35*12-1) 
  ;     nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
  ;;;;;
      data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_rc = regCoef_n(nino3_4,data_de, 0, 0)
      copy_VarCoords(data(1,:,:), data_rc)
      data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
      m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
      nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))
  
      data_relate_enso  = nino3_4_3d*data_rc_3d + m 
      data_new = data 
      data_new = (/data - data_relate_enso/)
  return data_new
end 

