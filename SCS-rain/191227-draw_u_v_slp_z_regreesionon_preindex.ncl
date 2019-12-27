;;********************************************************
;;回归ERA的u,v,slp,z到降水的变化上
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
    timeind=ind(date.ge.197901.and.date.le.201312)
    p=f2->precip(timeind,:,:)
    
    p = (/rmMonAnnCycTLL(p)/)
    p_annual = month_to_annual(p, 1)
    dim = dimsizes(p_annual)
    p_rc = regCoef_n(ispan(1, dim(0), 1), p_annual, 0, 0)
    copy_VarCoords(p_annual(1,:,:), p_rc)
    prob_p=rc_test(p_rc)
    p_rc=(/p_rc*dim(0)/)

    ;;;计算区域的降水平均，得到一个index,先看降水由哪里因素决定
     pre_index =dim_avg_n_Wrap(p_annual(:,{0:15},{110:155}), (/1,2/))
    ;;;
    delete(f2)
    delete(timeind)
    delete(date)
;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"uwnd","vwnd","slp","hgt"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ;;;region 
      lat1 = -10
      lat2 = 90
      lon1 = 0
      lon2 = 360
    ;;;

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
    slp = short2flt(f3->msl(timeind,{lat1:lat2},{lon1:lon2}))        

    f4 = addfile(filein(3),"r")
    date :=cd_calendar(f4->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    z = short2flt(f4->z(timeind,{500},{lat1:lat2},{lon1:lon2}))   
    print("=====finish read data==========")
   ;;;

    
   ;;;处理成年异常数据
     u = (/rmMonAnnCycTLL(u)/)
     v = (/rmMonAnnCycTLL(v)/)
     z = (/rmMonAnnCycTLL(z)/)
     slp = (/rmMonAnnCycTLL(slp)/)

     u_annual = month_to_annual(u, 1)
     v_annual = month_to_annual(v, 1)
     z_annual = month_to_annual(z, 1)
     slp_annual = month_to_annual(slp, 1)

   ;;;
  
   ;;;回归
   ; u_rc = regCoef_n(pre_index, u_annual, 0, 0)
   ; v_rc = regCoef_n(pre_index, v_annual, 0, 0)
   ; slp_rc = regCoef_n(pre_index, slp_annual, 0, 0)
   ; z_rc = regCoef_n(pre_index, z_annual, 0, 0)

   u_rc = regCoef_n(u_annual, pre_index, 0, 0)
   v_rc = regCoef_n(v_annual, pre_index, 0, 0)
   slp_rc = regCoef_n(slp_annual, pre_index, 0, 0)
   z_rc = regCoef_n(z_annual, pre_index, 0, 0)
   
   copy_VarCoords(u_annual(1,:,:), u_rc)
   copy_VarCoords(v_annual(1,:,:), v_rc)
   copy_VarCoords(slp_annual(1,:,:), slp_rc)
   copy_VarCoords(z_annual(1,:,:), z_rc)

   ;;;

   ;  ;;;回归
   ; u_rc = escorc_n(pre_index, u_annual, 0, 0)
   ; v_rc = escorc_n(pre_index, v_annual, 0, 0)
   ; slp_rc = escorc_n(pre_index, slp_annual, 0, 0)
   ; z_rc = escorc_n(pre_index, z_annual, 0, 0)
   
   ; copy_VarCoords(u_annual(1,:,:), u_rc)
   ; copy_VarCoords(v_annual(1,:,:), v_rc)
   ; copy_VarCoords(slp_annual(1,:,:), slp_rc)
   ; copy_VarCoords(z_annual(1,:,:), z_rc)
   ; ;;;  
   ;;;趋势一致性
   ;;;


   ; ;;;显著性检验
   ;  prob_slp = rc_test(slp_rc)
   ;  prob_u = rc_test(u_rc)
   ;  prob_v = rc_test(v_rc)
   ;  prob_z = rc_test(z_rc)
   ; ;;;
    prob_slp = rtest(slp_rc,35,0)
    prob_u = rtest(u_rc,35,0)
    prob_v = rtest(v_rc,35,0)
    prob_z = rtest(z_rc,35,0)

  

   slp_rc = where(prob_slp.le.0.05, slp_rc, slp_rc@_FillValue)

   u_rc = where(prob_u.le.0.05.or.prob_v.le.0.05, u_rc, u_rc@_FillValue)





;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "uv_slp_z_regression_on_pre-index-回归到各个变量场")
  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueWhiteOrangeRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -100
  ; res@cnMaxLevelValF       =   100
  ; res@cnLevelSpacingF      =  10
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
  res@vcMinMagnitudeF          =0.04
  res@vcRefMagnitudeF          = 1.0
  res@vcRefLengthF             = 0.04     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.0
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
 
  plot = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,slp_rc ,res)


  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

    ;;;画z的等直线
    res_z=True
    res_z@gsnDraw=False
    res_z@gsnFrame=False

    res_z@cnLinesOn=True
    ; res_z@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_z@cnMinLevelValF      = -140.0        ; set min contour level
    ; res_z@cnMaxLevelValF      = 140.0      ; set max contour level
    ; res_z@cnLevelSpacingF     = 40      ; set contour spacin
    res_z@cnLineLabelsOn     =True
    ; res_z@cnLineLabelBackgroundColor="Transparent"

    res_z@gsnContourNegLineDashPattern=2
    res_z@gsnContourPosLineDashPattern = 0

    res_z@cnLineLabelInterval        = 2             ; default = 2
    res_z@cnInfoLabelOn=False
    res_z@cnLineLabelFont=3.0
    res_z@cnLineLabelFontHeightF=0.008
    res_z@cnInfoLabelOn=False
    res_z@cnLineThicknessF=1.0
    res_z@gsnContourZeroLineThicknessF=0.001

    plot_z = gsn_csm_contour(wks , z_rc ,res_z)

    overlay(plot, plot_z)


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
    
  







    
   