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

    u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    
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
    print("=====finish read data==========")
   ;;;



   ;;；计算散度

    d =  uv2dv_cfd(u,v,u&latitude,u&longitude,2)
    
    copy_VarCoords(u, d)
    printVarSummary(d)

   ;;;计算大气不稳定度 ：假相当位温：difference between 1000 and 700 hPa
      p = t&level
      pre = p*100.0
      sigmae = pot_temp_equiv(pre,t,q,1,"q")
      printVarSummary(sigmae)
      
     sta_s = sigmae(:,{1000},:,:) - sigmae(:,{700},:,:) 
     copy_VarCoords(u(:,:,:), sta_s)
     printVarSummary(sta_s)
    ;;;

    
   ;;;处理成年异常数据
     u = (/rmMonAnnCycTLL(u)/)
     v = (/rmMonAnnCycTLL(v)/)
     d = (/rmMonAnnCycTLL(d)/)
     sta_s = (/rmMonAnnCycTLL(sta_s)/)



     u_annual = month_to_annual(u, 1)
     v_annual = month_to_annual(v, 1)
     d_annual = month_to_annual(d, 1)
     sta_s_annual = month_to_annual(sta_s, 1)
   ;;;
  
   ;;回归
   u_rc = regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
   v_rc = regCoef_n(ispan(1, yearnum, 1), v_annual, 0, 0)
   d_rc = regCoef_n(ispan(1, yearnum, 1), d_annual, 0, 0)
   sta_s_rc = regCoef_n(ispan(1, yearnum, 1), sta_s_annual, 0, 0)

   copy_VarCoords(u_annual(1,:,:), u_rc)
   copy_VarCoords(v_annual(1,:,:), v_rc)
   copy_VarCoords(d_annual(1,:,:), d_rc)
   copy_VarCoords(sta_s_annual(1,:,:), sta_s_rc)

   d_rc = (/d_rc*yearnum*10^8/)
   u_rc = (/u_rc*yearnum/)
   v_rc = (/v_rc*yearnum/)
   sta_s_rc = (/sta_s_rc*yearnum/)


   ;;;

   ;;;显著性检验
    
    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)
    prob_d = rc_test(d_rc)
    prob_sta_s = rc_test(sta_s_rc)
   ;;;
 
   ; slp_rc = where(prob_slp.le.0.1, slp_rc, slp_rc@_FillValue)
   ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "散度,uv以及大气净力稳定度_trend")
  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -5.0
  res@cnMaxLevelValF       =   5.0
  res@cnLevelSpacingF      =  1
  ;   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1
  ; res@cnMaxLevelValF       =   1
  ; res@cnLevelSpacingF      =  .1
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2
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
  ; res@vcRefMagnitudeF          = 1.0
  ; res@vcRefMagnitudeF          = 0.1
  res@vcRefLengthF             = 0.08     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.0
  res@vcLineArrowHeadMaxSizeF=0.005
  res@vcLineArrowHeadMinSizeF=0.005
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
  plot = new(2, graphic)
  plot(0) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,d_rc ,res)


  plot(1) = gsn_csm_contour_map(wks,sta_s_rc,res)


  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

 
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   i = 0
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   

   ; ;;;打点
   ;  res_pro=True
   ;  res_pro@gsnAddCyclic = False
   ;  res_pro@gsnDraw=False
   ;  res_pro@gsnFrame=False
   ;  res_pro@cnFillOn = True

   ;  res_pro@cnLevelSelectionMode = "ExplicitLevels"
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  res_pro@lbLabelBarOn = False
   ;  res_pro@cnLinesOn=False
   ;  ; res_pro@cnMonoLineColor =False
   ;  ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

   ;  res_pro@cnMonoFillPattern = True
   ;  res_pro@cnFillPattern = 17
   ;  res_pro@cnFillDotSizeF = 0.0045
   ;  res_pro@cnMonoFillColor = False
   ;  res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
   ;  ; res_pro@cnLevels = ispan(-6,6,1)
   ;  res_pro@cnLineLabelsOn=False
   ;  res_pro@tmYRMode = "Automatic"
   ;  res_pro@cnInfoLabelOn=False
   ;  ; res_pro@cnLineLabelBackgroundColor="Transparent"
   ;  ; res_pro@cnLineThicknessF = 3.
   ;  ; res_pro@cnFillOn        = True
   ;  iplot = new(1,graphic)

   ;  iplot(0)=gsn_csm_contour(wks,prob_slp,res_pro)
   

   ;  i = 0
   ;  overlay(plot(i), iplot(i))
    
    draw(plot)
    frame(wks)

end
    
  







    
   