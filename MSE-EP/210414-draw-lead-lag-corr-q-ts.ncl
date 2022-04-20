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
  yearstart = 195001
  yearend = 201912

  lat1 = -30
  lat2 = 30
  lon1 = 0
  lon2 = 360

;;;;read data
    var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(ind(date.eq.yearstart)))
    delete(date)     
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    q= q/1000.0

    hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;

    Cp = 1004
    Lv = 2.5*10^6
    g = 9.8
;;;;

;;; 
 

  ;;;MSE 
  plev = q&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)


  q_col =  dim_sum_wgt_n_Wrap(q, dz, 1, 1)

  q_col = (/rmMonAnnCycTLL(q_col)/)

  t_col =  dim_sum_wgt_n_Wrap(t, dz, 1, 1)

  t_col = (/rmMonAnnCycTLL(t_col)/)
  

  q_col_nep = dim_avg_n_Wrap(q_col(:,{0:20},{210:280}), (/1,2/))
  t_col_nep = dim_avg_n_Wrap(t_col(:,{0:20},{210:280}), (/1,2/))


  mxlag    = 12
   
  x =q_col_nep
  y = t_col_nep

  x_Lead_y = esccr(x,y,mxlag)
  y_Lead_x = esccr(y,x,mxlag)    ; switch the order of the series

  ccr = new ( 2*mxlag+1, float)    
  ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
  ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"


  wks=gsn_open_wks("pdf", "超前滞后相关")
  resp=True
  resp@gsnXRefLine = 0
  plot = gsn_csm_xy(wks,fspan(-12,12,25),ccr,resp)
  end 





