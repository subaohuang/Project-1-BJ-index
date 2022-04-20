   
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
   diri="/mnt/f/share/NOAA/"
   fils1=systemfunc("ls "+diri+"olr.mon.mean.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    year1=1979
    year2=2012
    yearend=201212
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    olr_r=f1->olr(timeind,:,:)
    olr=(/olr_r@add_offset+olr_r@scale_factor*olr_r/)
    copy_VarCoords(olr_r, olr)
    copy_VarAtts(olr_r, olr)
    delete(olr&time)
    olr&time=date(timeind)
    printVarSummary(olr)
    olr_pi=lonPivot(olr, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_djf_r=month_to_season(olr_pi(time|:,lat|:,lon|:), "MAM")

    ; olr_djf_r=(/month_to_annual(olr_pi(time|:,lat|:,lon|:), 1)/)
    olr_djf=olr_djf_r(1:,:,:)
    ;;;
    delete(date)
    delete(timeind)

    f2=addfile("/mnt/f/share/GPCP/precip.mon.mean.nc", "r")
    ; f2=addfile("/mnt/f/share/GPCC/precip.mon.total.v7.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi=lonPivot(p, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    p_djf_r=month_to_season(p_pi(time|:,lat|:,lon|:), "MAM")
    ; p_djf_r=(/month_to_annual(p_pi(time|:,lat|:,lon|:), 1)/)
    ;;注意点
    p_djf=p_djf_r(1:,:,:)

    delete(date)
    delete(timeind)
    delete(p)
    delete(p_pi)
    ;;;
    
    ;;;GPCC 陆地降水
    f3=addfile("/mnt/f/share/GPCC/precip.mon.total.v7.nc", "r")
    date=cd_calendar(f3->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    p=f3->precip(timeind,:,:)
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi=lonPivot(p, 30)
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    p_djf_r_land=month_to_season(p_pi(time|:,lat|:,lon|:), "MAM")
    ; p_djf_r_land=(/month_to_annual(p_pi(time|:,lat|:,lon|:), 1)/)

    ;;注意点
    p_djf_land=p_djf_r_land(1:,:,:)
    printVarSummary(p_djf_land)
    ;;;

    ;;;计算趋势
    olr_rc=regCoef_n(ispan(1, 34, 1), olr_djf, 0, 0)
    copy_VarCoords(olr_djf(1,:,:), olr_rc)
    prob=rc_test(olr_rc)
    olr_rc=(/olr_rc*10.0/)
    
    p_rc=regCoef_n(ispan(1, 34, 1), p_djf, 0, 0)
    copy_VarCoords(p_djf(1,:,:), p_rc)
    prob_p=rc_test(p_rc)
    p_rc=(/p_rc*10.0/)
    
    p_land_rc=regCoef_n(ispan(1, 34, 1), p_djf_land, 0, 0)
    copy_VarCoords(p_djf_land(1,:,:), p_land_rc)
    prob_p_land=rc_test(p_land_rc)
    p_land_rc=(/p_land_rc*10.0/)
    

    wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/MAM的OLR与GPCP降水的趋势")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(3,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnLeftString="(a) OLR"
  res@cnFillPalette="NCV_blu_red"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -8.0
  res@cnMaxLevelValF       =   8.0
  res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             =30
  res@mpMaxLonF             = 390
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBValues=(/30,90,150,210,270,330,387.5/)
  res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

  res@tmYLValues=(/30,20,10,0,-10,-20,-30/)
  res@tmYLLabels=(/"30~S~o~N~N","20~S~o~N~N","10~S~o~N~N","0"\
          ,"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S"/)
  ; res@tmYLLabels=(/"30~S~o~C~~Z130~~H5~N","20~S~o","10~S~o","0"\
  ;               ,"10~S~o","20~S~o~","~Z70~~H5~~V2~S~C~~Z100~~V20~30~S~o"/)

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  res@mpShapeMode="FreeAspect"
  res@vpHeightF=0.3
  res@vpWidthF=0.7
  plot(0)=gsn_csm_contour_map(wks, olr_rc, res)

   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.2
  res@cnMaxLevelValF       =   1.2
  res@cnLevelSpacingF      =  0.2

   res@cnFillPalette="MPL_BrBG"
   res@gsnLeftString="(b)GPCP Precipitation"
  plot(1)=gsn_csm_contour_map(wks, p_rc, res)
  delete(res@cnLevelSelectionMode)
  delete(res@cnMaxLevelValF)
  delete(res@cnMinLevelValF)
  delete(res@cnLevelSpacingF)
  delete(res@tmYLValues)
  delete(res@tmXBValues)
  delete(res@tmXBLabels)
  delete(res@tmYLLabels)

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -10
  res@mpMaxLatF             = 10
  res@mpMinLonF             =80
  res@mpMaxLonF             = 150
  ; res@mpMinLonF             =30
  ; res@mpMaxLonF             = 60
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -40
  res@cnMaxLevelValF       =   40
  res@cnLevelSpacingF      =  10

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@tmXBValues=ispan(80, 150, 10)
  res@tmXBLabels=(/"80~S~o~N~E","90~S~o~N~E","100~S~o~N~E","110~S~o~N~E"\
            ,"120~S~o~N~E","130~S~o~N~E","140~S~o~N~E","150~S~o~N~E"/)
  res@tmYLValues=(/10,5,0,-5,-10/)
  res@tmYLLabels=(/"10~S~o~N~N","5~S~o~N~N","0"\
            ,"5~S~o~N~S","10~S~o~N~S"/)
  
  ; res@tmXBLabels=(/"80E","90E","100E","110E"\
  ;         ,"120E","130E","140E","150E"/)
  res@tmXBMinorOn=True
  res@tmXBMinorValues=ispan(80,150,5)


  res@gsnLeftString="(c)GPCC Precipitation"
  plot(2)=gsn_csm_contour_map(wks, p_land_rc, res)

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
   plot@$str2$ = gsn_add_polyline(wks, plot(1), pgx, pgy, pgres)
   plot@$str3$ = gsn_add_polyline(wks, plot(2), pgx, pgy, pgres)
  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   pgx_pool=(/80,150,150,80,80/)
   pgy_pool=(/-10,-10,10,10,-10/)
   ;;;
     plot@$str3$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)
     plot@$str4$ = gsn_add_polyline(wks, plot(1), pgx_pool, pgy_pool, pgres)




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

    iplot=new(3,graphic)

    iplot(0)=gsn_csm_contour(wks,prob,res_pro)
 
    iplot(1)=gsn_csm_contour(wks,prob_p,res_pro)
    
    iplot(2)=gsn_csm_contour(wks,prob_p_land,res_pro)

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


   do i=0,2
    iplot(i)= gsn_contour_shade(iplot(i), 0,0.05,opt)
end do
     

    ;;;  
do i=0,2
overlay(plot(i), iplot(i))
end do

pres=True
gsn_panel(wks, plot, (/3,1/), pres)



end