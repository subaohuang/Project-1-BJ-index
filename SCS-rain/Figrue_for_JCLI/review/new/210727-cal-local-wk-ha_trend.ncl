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
    diri="/home/yangsong3/data-observation/linshh/data/era/"
    fils1=systemfunc("ls "+diri+"era_m_zonal_meri_ERAI_500hpa.nc")

    startmon = 197901
    endmon = 201512

    f1=addfile(fils1, "r")
    time = f1->time 
    date = cd_calendar(time, 1)

    timeind = ind(date.ge.startmon.and.date.le.endmon)

    m_zonal=f1->m_zonal(timeind,:,:)
    m_meri=f1->m_meri(timeind,:,:)

    m_zonal = (/m_zonal*10^4/)
    m_meri = (/m_meri*10^4/)
    
    m_zonal_clm = dim_avg_n_Wrap(m_zonal, 0)
    m_meri_clm = dim_avg_n_Wrap(m_meri, 0)

    m_zonal_500=rmMonAnnCycTLL(m_zonal)
    m_meri_500=rmMonAnnCycTLL(m_meri)
    
    ; m_zonal_annual=month_to_annual(m_zonal_500, 1)
    ; m_meri_annual=month_to_annual(m_meri_500, 1)


    season = (/"MAM","JJA","SON","DJF"/)

    title = (/"(a)","(b)","(c)","(d)"/)

    wks=gsn_open_wks("pdf", "local-wk-trend")
    wks1=gsn_open_wks("pdf", "local-hc-trend")

    plot=new(4,graphic)
    plot1=new(4,graphic)
    iplot = plot
    iplot1 = plot
    str = ispan(1,4,1)+""

do i=0,3
;;;;
    m_zonal_annual=month_to_season(m_zonal_500, season(i))
    m_meri_annual=month_to_season(m_meri_500, season(i))
    
    dims=dimsizes(m_zonal_annual)
;;;;;

 ;;;trend 
    m_zonal_rc = regCoef_n(ispan(1,dims(0),1),m_zonal_annual,0,0)  
    copy_VarCoords(m_zonal_annual(0,:,:), m_zonal_rc)
    m_zonal_rc = (/m_zonal_rc*dims(0)/)

    m_meri_rc = regCoef_n(ispan(1,dims(0),1),m_meri_annual,0,0)  
    copy_VarCoords(m_meri_annual(0,:,:), m_meri_rc)
    m_meri_rc = (/m_meri_rc*dims(0)/)

    print("趋势计算结束")


    ;;;显著性检验
    prob_m_zonal = rc_test(m_zonal_rc)
    prob_m_meri = rc_test(m_meri_rc)
    ;;;


    
    print("检验结束")


;;;画图 
 
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False

    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -30.0
    res@cnMaxLevelValF       =   30.0
    res@cnLevelSpacingF      =  5.0

    res@cnLineLabelsOn = False

    res@lbLabelBarOn = False

    res@cnFillPalette = "BlueWhiteOrangeRed"
  
    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = -5
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 100
    res@mpMaxLonF             = 170

    ;  res@mpMinLatF             = -30
    ; res@mpMaxLatF             = 30
    ; res@mpMinLonF             = 0
    ; res@mpMaxLonF             = 360
    res@mpCenterLonF = (/res@mpMinLonF + res@mpMaxLonF /)*0.5


    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF=0.025
    res@tmYLLabelFontHeightF=0.025
    res@tmXBMajorThicknessF=2.0
    res@tmYLMajorThicknessF=2.0
    ; res@tmBorderThicknessF=1.0
    res@tmXBValues = (/100,120,140,160/) 
    res@tmXBLabels = (/100,120,140,160/) +""+"E" 
    res@tmYLValues = (/0,10,20,30/) 
    res@tmYLLabels = (/"0","10N","20N","30N"/) 

    res@tmYLMinorValues = ispan(-5,30,5)

    res@tmXTOn=False
    res@tmYROn=False

    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.9
    ; res@vpHeightF = 0.2 
    ;;; *************************************************
    res@gsnLeftString = title(i)+" "+season(i)  +" Local WC trend"


    plot(i)= gsn_csm_contour_map(wks,m_zonal_rc,res)

    res@gsnLeftString = title(i)+" "+season(i)  +" Local HC trend"
    plot1(i)= gsn_csm_contour_map(wks1,m_meri_rc,res)
  ;*********************************************************
   ;plot clm
  ; wks1=gsn_open_wks("pdf", "ERA-in的SCS-local沃克及哈得来环流的气候态")
  ; plot_clm = new(2,graphic)
  ; res@cnMinLevelValF       =   -0.005
  ; res@cnMaxLevelValF       =   0.005
  ; res@cnLevelSpacingF      =  0.0005
  ; res@cnFillPalette = "NCV_blue_red"
  ; res@gsnLeftString = "local WC climatology"
  ; plot_clm(0) = gsn_csm_contour_map(wks1,m_zonal_clm, res)
  ; res@gsnLeftString = "local HC climatology"
  ; plot_clm(1) = gsn_csm_contour_map(wks1,m_meri_clm ,res)
  ; pres=True
  ; gsn_panel(wks1, plot_clm, (/2,1/), pres)

;;;打点
    res_pro4=True
    res_pro4@gsnAddCyclic = False
    res_pro4@gsnDraw=False
    res_pro4@gsnFrame=False
    res_pro4@cnFillOn = True
    res_pro4@cnLinesOn  = False
    res_pro4@cnLineLabelsOn = False

    res_pro4@cnLevelSelectionMode = "ExplicitLevels"
    res_pro4@cnLevels = (/0.0,0.1/)
    res_pro4@lbLabelBarOn = False

    res_pro4@cnLinesOn=False
    ; res_pro4@cnMonoLineColor =False
    ; res_pro4@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro4@cnMonoFillPattern = True
    res_pro4@cnFillPattern = 17
    res_pro4@cnFillDotSizeF = 0.001
    res_pro4@cnMonoFillColor = False
    res_pro4@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro4@cnLevels = ispan(-6,6,1)
    res_pro4@cnLineLabelsOn=False
    res_pro4@tmYRMode = "Automatic"
    res_pro4@cnInfoLabelOn=False

    iplot(i)=gsn_csm_contour(wks,prob_m_zonal,res_pro4)

      ; res_pro4@cnLevels = (/0.0,0.4/)
      ;  res_pro4@cnFillColors = (/"Transparent","Black","Transparent"/)
    iplot1(i)=gsn_csm_contour(wks1,prob_m_meri,res_pro4)

;;;  
 
  ; ;;添加区域
   pgres = True
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
  
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

   
   plot1@$str(i)$ = gsn_add_polyline(wks1, plot1(i), pgx_pool, pgy_pool, pgres)

  overlay(plot(i), iplot(i))
  overlay(plot1(i), iplot1(i))

  delete(m_zonal_rc)
  delete(m_meri_rc)
  
end do
  pres=True
  pres@gsnPanelLabelBar = True
  pres@lbLabelFontHeightF = 0.02
  gsn_panel(wks, plot, (/2,2/), pres)
  gsn_panel(wks1, plot1, (/2,2/), pres)

end 

