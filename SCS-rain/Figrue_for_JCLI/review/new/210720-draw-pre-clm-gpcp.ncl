;;********************************************************
;;画GPCP的降水趋势
;;
;;lnshuheng 2019.12.24
;**********************************************************

begin
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    ; f2=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201512)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi = p
    
    
    rc_mam = dim_avg_n_Wrap(month_to_season(p_pi, "MAM"), 0)
    rc_jja = dim_avg_n_Wrap(month_to_season(p_pi, "JJA"), 0)
    rc_son = dim_avg_n_Wrap(month_to_season(p_pi, "SON"), 0)
    rc_djf = dim_avg_n_Wrap(month_to_season(p_pi, "DJF"), 0)


wks   = gsn_open_wks("pdf", "不同季节GPCP-clm-1979-2015")
;;**************************************************************
;画图
;;***************************************************************
      
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnLeftStringFontHeightF=20.0
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
    res@lbLabelStride=2
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 100
    res@mpMaxLonF             = 170
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025
    res@gsnLeftStringFontHeightF = 0.03
    res@tmYLMinorOn = True

    res@lbBoxLinesOn              = True
    res@lbBoxLineThicknessF       = 0.


    res@tmXBValues = (/100,120,140,160/) 
    res@tmXBLabels = (/100,120,140,160/) +""+"E" 


    ; res@tmXBLabelFontHeightF=0.015
    ; res@tmYLLabelFontHeightF=0.015
    ; res@tmXBMajorThicknessF=1.0
    ; res@tmYLMajorThicknessF=1.0
    ; res@tmBorderThicknessF=1.0

    ; res@mpShapeMode="FreeAspect"
    ; res@vpHeightF=0.3
    ; res@vpWidthF=0.7


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   0
  res@cnMaxLevelValF       =   10
  res@cnLevelSpacingF      =  1

  res@cnFillPalette="CBR_wet"


  res@gsnLeftString="(a) MAM climatological Precip  " 
  plot = new(4,graphic)
  plot(0) = gsn_csm_contour_map(wks, rc_mam, res)
  res@gsnLeftString= "(b) JJA climatological Precip  "
  plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  res@gsnLeftString="(c) SON climatological Precip  "
  plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  res@gsnLeftString="(d) DJF climatological Precip  "
  plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 

  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=1

   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

     str2 = ispan(21, 24, 1) +""
   do i =0,3 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 


  
   
pres = True
pres@gsnPanelLabelBar  = True
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,2/), pres)

end