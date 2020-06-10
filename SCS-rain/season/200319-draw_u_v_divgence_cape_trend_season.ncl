;;********************************************************
;;计算ERA的散度趋势以及稳定度趋势
;;
;;lnshuheng 2020.03.17
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

;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"uwnd","vwnd","air","shum"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
    
   ;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160


    ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    u = short2flt(f1->u(timeind,{925},{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{925},{lat1:lat2},{lon1:lon2}))

    
    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    t = short2flt(f3->t(timeind,::-1,{lat1:lat2},{lon1:lon2}))        
    printVarSummary(t)


    f4 = addfile(filein(3),"r")
    date :=cd_calendar(f4->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    q = short2flt(f4->q(timeind,::-1,{lat1:lat2},{lon1:lon2}))        
    printVarSummary(q)

    f5 = addfile("/home/ys17-19/lsh/Project/SCS-rain/CAPE_PI_PIL-197901-201412.nc","r")
    date :=cd_calendar(f5->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    cape = short2flt(f5->cape(timeind,{lat1:lat2},{lon1:lon2}))        
    printVarSummary(cape)
    print("=====finish read data==========")
   ;;;



   ;;；计算散度

    d =  uv2dv_cfd(u,v,u&latitude,u&longitude,2)
    
    copy_VarCoords(u, d)
    printVarSummary(d)

   ;;;计算大气不稳定度 ：假相当位温：difference between 1000 and 700 hPa
      p = t&level
      pre = p*100.0
      printVarSummary(t)
      printVarSummary(q)

      sigmae = pot_temp_equiv(pre,t,q,1,"q")
      printVarSummary(sigmae)
      
     sta_s = sigmae(:,{1000},:,:) - sigmae(:,{950},:,:) 
   
     copy_VarCoords(u(:,:,:), sta_s)
     printVarSummary(sta_s)
    ;;;
  
;;;choose season
   season = (/"MAM","JJA","SON","DJF"/)
   plot  = new(4,graphic)
   iplot = new(4,graphic)

    plot1  = new(4,graphic)
    iplot1 = new(4,graphic)


   str = "" + ispan(1,4,1) + ""
   wks = gsn_open_wks("pdf", "散度trend-925hpa-scs")
   wks1   = gsn_open_wks("pdf", "cape-trend-925hpa-scs")


do i = 0,3
   ;;;处理成年数据 
     u_annual = month_to_season(u, season(i))
     v_annual = month_to_season(v, season(i))
     d_annual = month_to_season(d, season(i))
     sta_s_annual = month_to_season(sta_s, season(i))
     cape_annual = month_to_season(cape, season(i))
   ;;;

   ;;;区域平均 
      sigmae_scs = dim_avg_n_Wrap(sigmae(:,{1000:200},{20:35},{110:125}), (/2,3/))
      sigmae_scs_clm  = dim_avg_n_Wrap(sigmae_scs, 0)
      rc = regCoef_n(ispan(1, yearnum*12, 1),sigmae_scs, 0, 0)
      sigmae_scs_change = rc*yearnum*12 + sigmae_scs_clm
    ;;;

        ; wks = gsn_open_wks("pdf","ERA-In华南地区相当位温的气候态和趋势的profile")                  ; send graphics to PNG file
        ; resxy = True
        ; resxy@gsnDraw = False
        ; resxy@gsnFrame  =False
        ; resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
        ; resxy@xyMonoLineThickness = True
        ; resxy@xyLineColors       =  (/"black","blue"/)          ; change line color
        ; resxy@xyMarkLineModes    =  (/"lines","lines"/)
        ; resxy@xyExplicitLegendLabels  =  (/"main","dyn"/)     ;加线的标注
        ; resxy@trYReverse        = True               
        ; ; resxy@trYMinF = 100        ; reverse Y-axis

        ; plot1 = gsn_csm_xy(wks, (/sigmae_scs_clm,sigmae_scs_change/),sigmae_scs&level,resxy)
        ; draw(plot1)
        ; frame(wks)

 

   ;;回归
   u_rc := regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
   v_rc := regCoef_n(ispan(1, yearnum, 1), v_annual, 0, 0)
   d_rc := regCoef_n(ispan(1, yearnum, 1), d_annual, 0, 0)
   sta_s_rc := regCoef_n(ispan(1, yearnum, 1), sta_s_annual, 0, 0)
   cape_rc := regCoef_n(ispan(1, yearnum, 1), cape_annual, 0, 0)

   copy_VarCoords(u_annual(1,:,:), u_rc)
   copy_VarCoords(v_annual(1,:,:), v_rc)
   copy_VarCoords(d_annual(1,:,:), d_rc)
   copy_VarCoords(sta_s_annual(1,:,:), sta_s_rc)
   copy_VarCoords(cape_annual(1,:,:),cape_rc)

   d_rc = (/d_rc*yearnum*10^8/)
   u_rc = (/u_rc*yearnum/)
   v_rc = (/v_rc*yearnum/)
   sta_s_rc = (/sta_s_rc*yearnum/)
   cape_rc = (/cape_rc*yearnum/)

   ; sta_s_rc = dim_avg_n_Wrap(sta_s_annual, 0)


   ;;;

   ;;;显著性检验
    
    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)
    prob_d = rc_test(d_rc)
    prob_sta_s = rc_test(sta_s_rc)
    prob_cape_rc = rc_test(cape_rc)

   ;;;
 
   ; d_rc = where(prob_d.le.0.1, d_rc, d_rc@_FillValue)
   ; cape_rc = where(prob_cape_rc.le.0.1, cape_rc, cape_rc@_FillValue)
   ; sta_s_rc = where(prob_sta_s.le.0.1, sta_s_rc, sta_s_rc@_FillValue)
   ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

;==============================================================
; Mission 4. draw plot
; ==============================================================


  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette ="CBR_coldhot"

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -80.0
  res@cnMaxLevelValF       =   80.0
  res@cnLevelSpacingF      =  10.0
  ;res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1
  ; res@cnMaxLevelValF       =   1
  ; res@cnLevelSpacingF      =  .1
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@lbLabelStride = 2 
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = lat1
  ; res@mpMaxLatF             = lat2
  ; res@mpMinLonF             = lon1
  ; res@mpMaxLonF             = lon2
  res@mpMinLatF             = 0
  res@mpMaxLatF             = 15
  res@mpMinLonF             = 110
  res@mpMaxLonF             = 155

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
  
  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.02           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  res@vcRefMagnitudeF          = 3.0
  ; res@vcRefMagnitudeF          = 0.1
  res@vcRefLengthF             = 0.1     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.0
  res@vcLineArrowHeadMaxSizeF=0.005
  res@vcLineArrowHeadMinSizeF=0.005
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8

  res@gsnLeftString  = season(i)+"-Divergence trend"
  plot(i) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,d_rc ,res)
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -2.0
  ; res@cnMaxLevelValF       =   8.0
  ; res@cnLevelSpacingF      =  2.0
    res2 = res 
    res2@cnFillPalette  = "ViBlGrWhYeOrRe"

    res2@gsnSpreadColors      = True               ; use full color map
    res2@gsnSpreadColorStart  = 40              ; start at color 17
    res2@gsnSpreadColorEnd    = 100                ; end at color 200

    res2@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res2@cnMinLevelValF        = -100.0        ; set min contour level
    res2@cnMaxLevelValF        = 200.0          ; set max contour level
    res2@cnLevelSpacingF       = 25.0          ; set contour interval

  res@gsnLeftString  = season(i)+"-Convective available potential energy trend"
  plot1(i) = gsn_csm_contour_map(wks1,cape_rc,res2)
  ; delete(res@cnFillPalette)
  ; delete(res@gsnSpreadColorStart)
  ; delete(res@gsnSpreadColorEnd)
  ;  res@gsnSpreadColors      = False               ; use full color map
 

  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

 
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   plot1@$str(i)$ = gsn_add_polyline(wks1, plot1(i), pgx_pool, pgy_pool, pgres)


   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
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
    iplot(i)=gsn_csm_contour(wks,prob_d,res_pro)
    iplot1(i)=gsn_csm_contour(wks1,prob_cape_rc,res_pro)
   
    overlay(plot(i), iplot(i))
    overlay(plot1(i), iplot1(i))
end do 
    pres = True
    pres@gsnPanelLabelBar  = False
    gsn_panel(wks, plot, (/2,2/), pres)
    gsn_panel(wks1, plot1, (/2,2/), pres)


end
    
  







    
   