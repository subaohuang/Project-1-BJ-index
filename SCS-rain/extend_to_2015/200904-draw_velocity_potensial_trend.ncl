;;********************************************************
;;画HadIsst的sst趋势
;;
;;lnshuheng 2020.02.17
;**********************************************************
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
;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/era/"
    fils1=systemfunc("ls "+diri+"era_vp_sf.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201612.and.date.ge.197901)
    

    vp_850 = f1->vp(:,{925},:,:)
    vp_200 = f1->vp(:,{200},:,:)

    
    vp_850 = (/vp_850/10^6/) 
    vp_200 = (/vp_200/10^6/) 

    vp_200_annual = month_to_annual(vp_200, 1)
    vp_850_annual = month_to_annual(vp_850, 1)
    

    printMinMax(vp_850_annual, 1)

    dim = dimsizes(vp_850_annual)
    vp_850_rc=regCoef_n(ispan(1, dim(0), 1), vp_850_annual, 0, 0)
    copy_VarCoords(vp_850_annual(1,:,:), vp_850_rc)

    vp_200_rc=regCoef_n(ispan(1, dim(0), 1), vp_200_annual, 0, 0)
    copy_VarCoords(vp_200_annual(1,:,:), vp_200_rc)

    prob_p_850 = rc_test(vp_850_rc)
    prob_p_200 = rc_test(vp_200_rc)


    vp_850_rc = (/vp_850_rc*dim(0)/)
    vp_200_rc = (/vp_200_rc*dim(0)/)
    
    ; ssta_rc = where(prob_p.le.0.05, ssta_rc, ssta_rc@_FillValue)
  
 
  
    wks   = gsn_open_wks("pdf", "年平均velocity_potensial趋势-SCS-1979-2016")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(2,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  
  ; res@cnFillPalette = "NCV_jaisnd"
  res@cnFillPalette = "precip4_diff_19lev"
  res@cnFillOn=True
  res@cnLinesOn=True

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -30
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             =0
  ; res@mpMaxLonF             = 360

  res@mpMinLatF             = -10
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 70
  res@mpMaxLonF             = 170

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  ; res@tmXBValues=(/30,90,150,210,270,330,387.5/)
  ; res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
  ;             ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

  ; res@tmYLValues=(/30,20,10,0,-10,-20,-30/)
  ; res@tmYLLabels=(/"30~S~o~N~N","20~S~o~N~N","10~S~o~N~N","0"\
  ;       ,"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S"/)
  ; res@tmYLLabels=(/"30~S~o~C~~Z130~~H5~N","20~S~o","10~S~o","0"\
  ;               ,"10~S~o","20~S~o~","~Z70~~H5~~V2~S~C~~Z100~~V20~30~S~o"/)

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   1.0
  res@cnMaxLevelValF       =   1.6
  res@cnLevelSpacingF      =  0.1
    
  res@gsnLeftString="850hpa velocity potensial"
  res@cnFillColors = (/11,11,12,14,16,17,19,19/)

  plot(0) = gsn_csm_contour_map(wks,vp_850_rc, res)

  ; res@cnMinLevelValF       =   -4.0
  ; res@cnMaxLevelValF       =   0.0
  ; res@cnLevelSpacingF      =  0.5
   ; res@cnFillColors = (/2,3,4,6,8,11,12/)
  plot(1) = gsn_csm_contour_map(wks,vp_200_rc, res)

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    


  ; ;;;添加赤道直线
  ;  pgx=ispan(30,390,10)
  ;  pgy=new(dimsizes(pgx),typeof(pgx))
  ;  pgy=0
  ;  str1=unique_string("dum1")
  ;  str2=unique_string("dum2")
  ;   str3=unique_string("dum3")
  ;  pgres=True
  ;  pgres@gsLineColor="black"
  ;  pgres@gsLineThicknessF=1.0
  ;  pgres@gsLineDashPattern=1
  ;  plot@$str1$ = gsn_add_polyline(wks, plot(0), pgx, pgy, pgres)
  ; ;;;添加indo-pacific warm pool
  ;   ; ;;添加区域
  ;  ;;warm pool
  ;  pgres@gsLineThicknessF=2.0
  ;  pgres@gsLineDashPattern=0
  ;  str3=unique_string("dum3")
  ;  str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1
    str3 = "" + ispan(1, 2, 1)
    plot@$str3(0)$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)
    plot@$str3(1)$ = gsn_add_polyline(wks, plot(1), pgx_pool, pgy_pool, pgres)
    
   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.003
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  
  
     ; iplot = gsn_csm_contour(wks,prob_p  ,res_pro)

    ;;;  
    ; overlay(plot, iplot)
 
 pres = True
 gsn_panel(wks, plot, (/2,1/), pres)


end