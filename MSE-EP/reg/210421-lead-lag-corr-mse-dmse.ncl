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

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

begin
  ; yearstart = 200001
  ; yearend = 201912
    yearstart = 195001
  yearend = 201912

  lat1 = -30
  lat2 = 30
  lon1 = 0
  lon2 = 360


 f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")
 date=cd_calendar(f1->time, 1)
 timeind=ind(date.ge.yearstart.and.date.le.yearend)
 MSE = f1->MSE(timeind,:,:,:)

  ;;;MSE 
  plev = MSE&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

  MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

  MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

  MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

  ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)

  MSE_col_nep = dim_avg_n_Wrap(MSE_col(:,{0:20},{240:280}), (/1,2/))
  MSE_col_nep = (/MSE_col_nep/10^6/)


  ; MSE_col_MAM = month_to_season(MSE_col, "MAM")
  ; MSE_col_JJA = month_to_season(MSE_col, "JJA")
  ; MSE_col_SON = month_to_season(MSE_col, "SON")

 
  MSE_col_tendency = MSE_col_nep
  MSE_col_tendency = center_finite_diff_n(MSE_col_nep, 1, False, 1, 0)

   mxlag    = 12
   
  x = MSE_col_nep
  y = MSE_col_tendency

  x_Lead_y = esccr(x,y,mxlag)
  y_Lead_x = esccr(y,x,mxlag)    ; switch the order of the series

  ccr = new ( 2*mxlag+1, float)    
  ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
  ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"


  wks=gsn_open_wks("pdf", "MSE和DMSE超前滞后相关")
  resp=True
  resp@gsnXRefLine = 0
  plot = gsn_csm_xy(wks,fspan(-12,12,25),ccr,resp)
  end 

