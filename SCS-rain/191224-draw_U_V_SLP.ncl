;;********************************************************
;;画U,V,SLP的的趋势
;; ERA-interim
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
;-----------------------------------------------------------------
begin
    ;;;read 
    path = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"uwnd","vwnd","slp"/)
    filein = path + var_name +".mon.mean.nc"

    
    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)

    u = short2flt(f1->u(timeind,{850},:,:))
    
    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{850},:,:))
    
    f3 = addfile(filein(2),"r")
    slp = short2flt(f3->msl(timeind,:,:))
    ;;;
     
    ;;;处理成年平均的时间序列

    u_annual = month_to_annual(u, 1)
    v_annual = month_to_annual(v, 1)
    slp_annual = month_to_annual(slp, 1)
   
    ;;;
    
    ;;;计算趋势

    dim = dimsizes(u_annual)
    u_rc = regCoef_n(ispan(1, dim(0), 1), u_annual, 0, 0)
    copy_VarCoords(u_annual(1,:,:), u_rc)

    dim = dimsizes(v_annual)
    v_rc = regCoef_n(ispan(1, dim(0), 1), v_annual, 0, 0)
    copy_VarCoords(v_annual(1,:,:), v_rc)
    
    dim = dimsizes(slp_annual)
    slp_rc = regCoef_n(ispan(1, dim(0), 1), slp_annual, 0, 0)
    copy_VarCoords(slp_annual(1,:,:), slp_rc)
    
    ;;;

    ;;;显著性检验
    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)
    prob_slp = rc_test(slp_rc)

    ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)
    ; v_rc = where(prob_v.le.0.1.or.prob_v.le.0.1, v_rc, v_rc@_FillValue)
  
 
    ;;;转化成35year-1
    u_rc=(/u_rc*dim(0)/)
    v_rc=(/v_rc*dim(0)/)
    slp_rc=(/slp_rc*dim(0)/)
  


    wks   = gsn_open_wks("pdf", "年平均-era的U-V-SLP的趋势")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(1,graphic)

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnLeftString="(a) OLR"
  res@cnFillPalette="NCV_blu_red"
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
  res@mpMaxLonF             = 160
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -100
  res@cnMaxLevelValF       =   100
  res@cnLevelSpacingF      =  10

  res@cnFillPalette="BlueWhiteOrangeRed"

  res@gsnScalarContour        = True

  ;;;LineArrow
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
 
  plot = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,slp_rc ,res)
  
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

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

    
 
    iplot = gsn_csm_contour(wks,prob_slp,res_pro)
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


  
    iplot = gsn_contour_shade(iplot , 0,0.05,opt)

    ;;;  
    overlay(plot, iplot)

    draw(plot)
    frame(wks)



end