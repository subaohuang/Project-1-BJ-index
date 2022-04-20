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
;;;;
begin
 ;;;time span
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   ; lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160
   lat1 = -5
   lat2 = 30
   lon1 = 100
   lon2 = 170
;;;;
;;;read NCEP-NCAR heat flux data
    
    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/annual/ERAI-heatflux-1degree.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
  
    lhtfl = f3->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f3->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f3->nlwrs(timeind,{lat1:lat2},{lon1:lon2})

    lhtfl = (/lhtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    shtfl = (/shtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nswrs = (/nswrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nlwrs = (/nlwrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2


    nswrs_MAM = month_to_season(nswrs, "MAM")
    nswrs_JJA = month_to_season(nswrs, "JJA")
    nswrs_SON = month_to_season(nswrs, "SON")
    nswrs_DJF = month_to_season(nswrs, "DJF")

  
   

  
print("*******开始计算趋势系数*****************")
; ;;;calculate correlation coef
;       lhtfl_annual = lhtfl_annual/(365*86400)   ;;; w/m^2
;       shtfl_annual = shtfl_annual/(365*86400)   ;;; w/m^2
;       nswrs_annual = nswrs_annual/(365*86400)   ;;; w/m^2
;       nlwrs_annual = nlwrs_annual/(365*86400)   ;;; w/m^2
; ;;;

;;;
    ;;;大气

    r_ns_MAM =  regCoef_n(ispan(1, yearnum, 1), nswrs_MAM, 0, 0)

    r_ns_JJA =  regCoef_n(ispan(1, yearnum, 1), nswrs_JJA, 0, 0)
    r_ns_SON =  regCoef_n(ispan(1, yearnum, 1), nswrs_SON, 0, 0)
    r_ns_DJF =  regCoef_n(ispan(1, yearnum, 1), nswrs_DJF, 0, 0)

   
    copy_VarCoords(nswrs(0,:,:), r_ns_MAM)
    copy_VarCoords(nswrs(0,:,:), r_ns_JJA)
    copy_VarCoords(nswrs(0,:,:), r_ns_SON)
    copy_VarCoords(nswrs(0,:,:), r_ns_DJF)

    prob_mam = rc_test(r_ns_MAM)
    prob_jja = rc_test(r_ns_JJA)
    prob_son = rc_test(r_ns_SON)
    prob_djf = rc_test(r_ns_DJF)

    

    print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "短波辐射通量的线性趋势-ERAI-fourseason-pilipian")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="MPL_seismic"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillDrawOrder = "PostDraw"


  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             =lon1
  res@mpMaxLonF             = lon2

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.02
  res@tmYLLabelFontHeightF=0.02
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.1
 

  plot = new(4,graphic)
  res@gsnLeftString = "(a)MAM ShortWave radiation trend"
  plot(0) = gsn_csm_contour_map(wks, r_ns_MAM({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "(b)JJA ShortWave radiation trend"
  plot(1) = gsn_csm_contour_map(wks, r_ns_JJA({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "(c)SON ShortWave radiation trend"
  plot(2) = gsn_csm_contour_map(wks, r_ns_SON({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "(d)DJF ShortWave radiation trend"
  plot(3) = gsn_csm_contour_map(wks, r_ns_DJF({lat1:lat2},{lon1:lon2}), res)

  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,4,1)

   do i = 0,3
     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   end do


  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   str2 = "0" + ispan(1,4,1)
   do i = 0,3
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do

    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1      ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(4,graphic)

    iplot(0) = gsn_csm_contour(wks,prob_mam({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_jja({lat1:lat2},{lon1:lon2}) ,res_pro)
    iplot(2) = gsn_csm_contour(wks,prob_son({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(3) = gsn_csm_contour(wks,prob_djf({lat1:lat2},{lon1:lon2}),res_pro)
 

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


    do i = 0,3
     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    end do

    ;;;  
    do i = 0,3
      overlay(plot(i), iplot(i))
    end do 


  pres = True
  pres@gsnPanelLabelBar = True
  ; pres@lbOrientation = "Vertical"
  pres@lbLabelFontHeightF = 0.01

  pres@lbOrientation = "Vertical"
  pres@pmLabelBarHeightF = 0.3
  pres@lbLabelStride = 2

  gsn_panel(wks, plot, (/2,2/), pres)
 


end