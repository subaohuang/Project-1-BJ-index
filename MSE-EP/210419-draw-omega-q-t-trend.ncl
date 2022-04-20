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

function mse_trend(data)
begin
  rc1 =regCoef_n(ispan(1,dimsizes(data&time),1), data,0,0)
  copy_VarCoords(data(1,:,:), rc1)
  rc1 = rc1*dimsizes(data&time)
  return rc1
end 

begin
  yearstart = 197901
  yearend = 201912

  lat1 = -30
  lat2 = 30
  lon1 = 0
  lon2 = 360

; ;;;;read data
    var= (/"air.mon.mean.nc","shum.mon.mean.nc","omega.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{-50:50},{lon1:lon2})

    ; q= q/1000.0


    time = f3->time 
    date := cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    omega = f3->omega(timeind,{1000:300},{-50:50},{lon1:lon2})
    omega = (/omega*-10^2/)

;;;
  plev = q&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

    
;;;choose season
    plot  = new(4,graphic)
    iplot = new(4,graphic)

    plot1  = new(4,graphic)
    iplot1 = new(4,graphic)


    plot2  = new(4,graphic)
    iplot2 = new(4,graphic)


    str = "" + ispan(1,4,1) + ""
    wks1 = gsn_open_wks("pdf", "q-trend-79-19.ncl")
    wks2 = gsn_open_wks("pdf", "w-trend-79-19.ncl")

   level = (/925,750,925,750/)


do i=0,3 
  if i.gt.1 then 
   omega_level = omega(:,{level(i)},:,:)
   omega_level  = (/rmMonAnnCycTLL(omega_level)/)
   omega_level = (/dtrend_msg_n(ispan(1,dimsizes(omega_level&time),1),omega_level,False,True,0)/)
   omega_mam = month_to_season(omega_level, "MAM")
   rc_mam := mse_trend(omega_mam)
   prob_mam := rc_test(rc_mam)

 else
    q_col = q(:,{level(i)},:,:)
    q_col = (/rmMonAnnCycTLL(q_col)/)
    q_col_MAM  = month_to_season(q_col, "MAM")
    rc_mam := mse_trend(q_col_MAM)
    prob_mam := rc_test(rc_mam)
end if

;;;;

;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn=True

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  res@mpFillDrawOrder = "PostDraw"


  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "NCV_blu_red"
  if i.lt.2 then 
  res@gsnLeftString="MAM"+ level(i)+"hPa q trend"
  plot(i) = gsn_csm_contour_map(wks1, rc_mam ,res)
  else
  res@gsnLeftString="MAM"+ level(i)+"hPa omega trend"
  plot(i) = gsn_csm_contour_map(wks2, rc_mam ,res)
  end if

  
  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 


   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 210
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

 
   if i.lt.2 then 
   plot@$str(i)$ = gsn_add_polyline(wks1, plot(i), pgx_pool, pgy_pool, pgres)
   else
   plot@$str(i)$ = gsn_add_polyline(wks2, plot(i), pgx_pool, pgy_pool, pgres)
   end if


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

   if i.lt.2 then
    iplot(i) = gsn_csm_contour(wks1,prob_mam,res_pro)    
  else 
      iplot(i) = gsn_csm_contour(wks2,prob_mam,res_pro)    
    end if
    ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)    
    ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)    
    ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)    

    ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)
    ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)
    ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)

    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數
     ; do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
     overlay(plot(i), iplot(i))
end do

    ;;; 
    pres = True
    pres@gsnPanelLabelBar  = False
    pres@pmLabelBarOrthogonalPosF = -0.01
    gsn_panel(wks1, (/plot(0),plot(1)/), (/2,1/), pres)
    gsn_panel(wks2, (/plot(2),plot(3)/), (/2,1/), pres)

end


