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
    var_name = (/"uwnd","vwnd","slp"/)
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
   lon2 = 170


    ;;;time span
    startmon =197901
    endmon = 201612

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    
    f3 = addfile(filein(2),"r")
    ; f3 = addfile("/home/ys17-19/lsh/data/Hadley-slp/slp.mnmean.real.nc","r")
    date :=cd_calendar(f3->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    slp = short2flt(f3->msl(timeind,{lat1:lat2},{lon1:lon2}))      
    ; slp = short2flt(f3->slp(timeind,{lat1:lat2},{lon1:lon2}))        
    printVarSummary(slp)
   
    ; slp = (/slp*100/)
  
;;;choose season
    season = (/"MAM","JJA","SON","DJF"/)
    plot  = new(4,graphic)
    iplot = new(4,graphic)

    plot1  = new(4,graphic)
    iplot1 = new(4,graphic)


    plot2  = new(4,graphic)
    iplot2 = new(4,graphic)


    str = "" + ispan(1,4,1) + ""
    wks = gsn_open_wks("pdf", "u,v,slp趋势-scs-16")



do i = 0,3
   ;;;处理成年数据 
     u_annual = month_to_season(u, season(i))
     v_annual = month_to_season(v, season(i))
     slp_annual = month_to_season(slp, season(i))
   ;;;


   ;;趋势
   u_rc := regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
   v_rc := regCoef_n(ispan(1, yearnum, 1), v_annual, 0, 0)
   slp_rc := regCoef_n(ispan(1, yearnum, 1), slp_annual, 0, 0)
 
   copy_VarCoords(u_annual(1,:,:), u_rc)
   copy_VarCoords(v_annual(1,:,:), v_rc)
   copy_VarCoords(slp_annual(1,:,:), slp_rc)
 
   slp_rc = (/slp_rc*yearnum/)
   u_rc = (/u_rc*yearnum/)
   v_rc = (/v_rc*yearnum/)
  


   ;;;显著性检验
    
    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)
    prob_slp = rc_test(slp_rc)
   
   ;;;
 
   ; d_rc = where(prob_d.le.0.1, d_rc, d_rc@_FillValue)
   ; cape_rc = where(prob_cape_rc.le.0.1, cape_rc, cape_rc@_FillValue)
   ; sta_s_rc = where(prob_sta_s.le.0.1, sta_s_rc, sta_s_rc@_FillValue)

   u_rc_test = u_rc
   v_rc_test = v_rc
   u_rc_test = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

;==============================================================
; Mission 4. draw plot
; ==============================================================

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
 
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "BlueRed"

  res@cnFillOn=True
  res@cnLinesOn= True
  
  res@cnLineThicknessF = 2.0
  res@gsnContourNegLineDashPattern = 14
  res@gsnContourPosLineDashPattern = 0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -60.0
  res@cnMaxLevelValF       =   60.0
  res@cnLevelSpacingF      =  10.0

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@lbLabelStride = 2 
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             = 125
  ; res@mpMaxLonF             = 160
  res@mpMinLatF             = -5
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 100
  res@mpMaxLonF             = 170
  ;  ;   lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160


  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.5
  
  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.022           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  res@vcRefMagnitudeF          = 4.0
  res@vcRefLengthF             = 0.1     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =True
  res@mpFillColor = "White"
  res@mpFillDrawOrder = "PostDraw"

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8

  res@gsnLeftString  = season(i)+" SLP trend"
  plot(i) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,slp_rc ,res)



  ;;;;;;;LineArrow
  resv = True
  resv@gsnDraw = False
  resv@gsnFrame = False
  resv@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  resv@vcMinDistanceF           = 0.022           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  resv@vcRefMagnitudeF          = 4.0
  resv@vcRefLengthF             = 0.1     ; define length of v
  resv@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  resv@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  resv@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  resv@vcLineArrowHeadMinSizeF = 0.006
  resv@vcVectorDrawOrder="PostDraw"

  resv@vcMonoLineArrowColor = True
  resv@vcLineArrowColor = "green"
  resv@vcLineArrowThicknessF=1.5

  resv@vcRefAnnoOn    = False
  resv@gsnLeftString         = " "
  resv@gsnRightString        = " "
  

  plot_vect := gsn_csm_vector(wks,u_rc_test({-5:30},{100:170}),v_rc_test({-5:30},{100:170}),resv)
  overlay (plot(i),plot_vect)

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -2.0
  ; res@cnMaxLevelValF       =   8.0
  ; res@cnLevelSpacingF      =  2.0
 
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
 

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
    iplot(i)=gsn_csm_contour(wks,prob_slp,res_pro)
  
    overlay(plot(i), iplot(i))

end do 
    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar  = True
    pres@lbOrientation = "Vertical"
  
    gsn_panel(wks, plot, (/2,2/), pres)



end
    
  







    
   