;;********************************************************
;;画GPCP的降水趋势
;;
;;lnshuheng 2019.12.24
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
   begin

    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi = p

   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    ; p_djf_r=month_to_season(p_pi(time|:,lat|:,lon|:), "MAM")
    p_djf_r=month_to_annual(p_pi(time|:,lat|:,lon|:), 1)
    ;;注意点
    p_djf=p_djf_r(:,:,:)

    delete(date)
    delete(timeind)
    delete(p)
    delete(p_pi)
    ;;;
    
    dim = dimsizes(p_djf_r)
    p_rc=regCoef_n(ispan(1, dim(0), 1), p_djf, 0, 0)
    copy_VarCoords(p_djf_r(1,:,:), p_rc)
    prob_p=rc_test(p_rc)
    p_rc=(/p_rc*dim(0)/)
    
    wks   = gsn_open_wks("pdf", "年平均GPCP降水的趋势-scs-1979-2016")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  ; res@cnFillPalette="NCV_blu_red"
  res@cnFillPalette="MPL_BrBG"
  res@cnFillPalette =  "precip_diff_12lev"
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -5
  res@mpMaxLatF             = 45
  res@mpMinLonF             =90
  res@mpMaxLonF             = 165
  ; res@mpMinLatF             = -60
  ; res@mpMaxLatF             = 60
  ; res@mpMinLonF             =0
  ; res@mpMaxLonF             = 360



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

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.7

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.5
  res@cnMaxLevelValF       =   1.5
  res@cnLevelSpacingF      =  0.25


  res@gsnLeftString="GPCP Precipitation trend"
  plot = gsn_csm_contour_map(wks, p_rc, res)
  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
  ; ;;添加区域
  ;  ;;indian
  ;  pgx_indian=(/70,100,100,70,70/)
  ;  pgy_indian=(/-5,-5,5,5,-5/)
  ;  ;;;
  ;  ;;;pacific
  ;  pgx_pa=(/150,210,210,150,150/)
  ;  pgy_pa=(/-5,-5,5,5,-5/)
  ;  ;;;
  ;  ;;;altantic
  ;  pgx_alt=(/320,340,340,320,320/)
  ;  pgy_alt=(/-5,-5,5,5,-5/)
  ;  ;;;
   
  ;  str1=unique_string("dum1")
  ;  str2=unique_string("dum2")
  ;  str3=unique_string("dum3")

  ;  pgres=True
  ;  pgres@gsLineColor="black"
  ;  pgres@gsLineThicknessF=2.0
   
   
  ; plot@$str1$ = gsn_add_polyline(wks, plot(0), pgx_indian, pgy_indian, pgres)
  ; plot@$str2$ = gsn_add_polyline(wks, plot(0), pgx_pa, pgy_pa, pgres)
  ; plot@$str3$ = gsn_add_polyline(wks, plot(0), pgx_alt, pgy_alt, pgres)
  ; ;;

  ;;;添加赤道直线

    
   pgx=ispan(30,390,10)
   pgy=new(dimsizes(pgx),typeof(pgx))
   pgy=0
   str1=unique_string("dum1")
   str2=unique_string("dum2")
    str3=unique_string("dum3")
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1
   plot@$str1$ = gsn_add_polyline(wks, plot(0), pgx, pgy, pgres)
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
     plot@$str3$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)

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

    iplot=new(1,graphic)



    iplot = gsn_csm_contour(wks,prob_p,res_pro)

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeMid     = 6         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.1         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.2      ; make dots larger
  iplot = gsn_contour_shade(iplot, 0,0.1, opt)
  overlay(plot, iplot)

  ; res_pro=True
  ; res_pro@gsnDraw=False
  ; res_pro@gsnFrame=False
  ; res_pro@cnFillOn = True

  ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
  ; res_pro@cnLevels = (/0,0.05/)
  ; res_pro@lbLabelBarOn = False
  ; res_pro@cnLinesOn=False
  ; res_pro@cnMonoLineColor =False
  ; res_pro@cnLineColors = ((/"Transparent","black","Transparent"/))

  ; res_pro@cnMonoFillPattern = True
  ; res_pro@cnFillPattern = 17
  ; res_pro@cnFillDotSizeF = 0.0045
  ; res_pro@cnMonoFillColor = False
  ; res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
  ; ; res_pro@cnLevels = ispan(-6,6,1)
  ; res_pro@cnLineLabelsOn=False
  ; res_pro@tmYRMode = "Automatic"
  ; res_pro@cnInfoLabelOn=False
  ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; ; res_pro@cnLineThicknessF = 3.
  ; ; res_pro@cnFillOn        = True
  ; iplot(0)=gsn_csm_contour(wks,prob_p,res_pro)
  ; overlay(plot, iplot)
  ;;;  

draw(plot)
frame(wks)


end