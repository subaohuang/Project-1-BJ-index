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

function mse_diff(data,year1,year2,year3,year4)
begin
  ; rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  ; copy_VarCoords(data(1,:,:), rc1)
  ; rc1 = rc1*(year2-year1+1)
  rc2 = dim_avg_n_Wrap(data({year1:year2},:,:),0)
  rc2 = dim_avg_n_Wrap(data({year3:year4},:,:),0) - dim_avg_n_Wrap(data({year1:year2},:,:),0)
  return rc2
end 


begin
  yearstart = 195001
  yearend = 201912

  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360

  f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")

  time = f1->time 
  date = cd_calendar(time, 1)
  delete(time)
  timeind :=ind(date.le.yearend.and.date.ge.yearstart)

  MSE = f1->MSE(timeind,:,:,:)


;;;MSE 
  plev = MSE&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)


  MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

  MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

  MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))
  ; MSE_nep = MSE_col
  MSE_nep = (/MSE_nep/10^6/)

  MSE_nep_annual = month_to_annual(MSE_nep, 1)
  MSE_nep_annual!0 = "time"

  MSE_nep_MAM = month_to_season(MSE_nep, "MAM")
  MSE_nep_JJA = month_to_season(MSE_nep, "JJA")
  MSE_nep_SON = month_to_season(MSE_nep, "SON")
  MSE_nep_DJF = month_to_season(MSE_nep, "DJF")


  MSE_nep_DJF(0:dimsizes(MSE_nep_DJF&time)-2) = MSE_nep_DJF(1:)
  MSE_nep_DJF(dimsizes(MSE_nep_DJF&time)-1) = MSE_nep(dimsizes(MSE_nep&time)-1)

 
 ;;;计算趋势
  ; MSE_nep_annual&time = ispan(1950,2019,1)
  ; MSE_nep_MAM&time = ispan(1950,2019,1)
  ; MSE_nep_JJA&time = ispan(1950,2019,1)
  ; MSE_nep_SON&time = ispan(1950,2019,1)
  ; MSE_nep_DJF&time = ispan(1950,2019,1)

  MSE_nep_annual&time = ispan(1950,2019,1)
  MSE_nep_MAM&time = ispan(1950,2019,1)
  MSE_nep_JJA&time = ispan(1950,2019,1)
  MSE_nep_SON&time = ispan(1950,2019,1)
  MSE_nep_DJF&time = ispan(1950,2019,1)
  
  ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1
  
  totalyear = ispan(1950, 2019, 1)
  ely_out_ind = get1Dindex_Exclude(MSE_nep_MAM&time, ely)
  noninoyear = totalyear(ely_out_ind)
  x = ispan(1,dimsizes(ely_out_ind), 1)

  MSE_nep_MAM_enino = MSE_nep_MAM(ely_out_ind)
  MSE_nep_JJA_enino = MSE_nep_JJA(ely_out_ind)
  MSE_nep_SON_enino = MSE_nep_SON(ely_out_ind)
  MSE_nep_DJF_enino = MSE_nep_DJF(ely_out_ind)

    
  wks   = gsn_open_wks("pdf", "非ENnioMSE的时间序列")

  res               = True                   ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False
  res@tiMainString  = " MSE anomalies (80W-150W,0-20N) during no Enino year"   ; title
  res@tiMainFontHeightF = 0.02
  ; res@tiYAxisFontHeightF = 0.02
  ; res@tiXAxisFontHeightF = 0.02
  res@tiYAxisString = "MSE anomaly"           ; y axis title
  res@tiXAxisString = "Year"

  res@xyLineColors      = (/"black","red","blue","green"/)  ; line colors
   res@xyMonoLineThickness = True
  res@xyLineThicknessF = (/4.0/)        ; line thicknesses
  res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
  res@gsnYRefLine = 0.0

  ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
  ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
  ; res@xyMarkers         = (/2,0/)          ; style

  res@tmXBMode="Explicit"

  res@trXMaxF = dimsizes(ely_out_ind)+1
  ; res@trXMinF = 1950
  res@trXMinF = 0
  ; res@trYMinF = 0.2
  ; res@trYMaxF = 1.8

  res@tmXTOn = False
  res@tmYROn  = False

  res@tmXBValues= x
  res@tmXBLabels= noninoyear
  res@tmXBLabelAngleF = -60.0
  res@tmXBLabelJust = "CenterCenter"
  res@tmXBLabelFontHeightF = 0.01
  res@tmXBMajorLengthF = 0.01
  res@tmXBMajorOutwardLengthF = 0.01
  res@tmYRMode="Automatic"


  plot = gsn_csm_xy(wks, x, (/MSE_nep_MAM_enino,MSE_nep_JJA_enino,MSE_nep_SON_enino,MSE_nep_DJF_enino/), res)


    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"MAM","JJA","SON","DJF"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot = simple_legend(wks,plot,gres,lineres,textres)   

draw(plot)
frame(wks)
end
