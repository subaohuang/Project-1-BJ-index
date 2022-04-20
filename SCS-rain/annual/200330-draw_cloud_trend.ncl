;;********************************************************
;;计算ICLOUD的云量趋势
;;
;;lnshuheng 2020.03.17
;**********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/crop.ncl"
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
    ;;;Icoads 云量
    f4=addfile("/home/ys17-19/lsh/data/ICOADS/cldc.mean.nc", "r")
    date=cd_calendar(f4->time, 1)
    timeind=ind(date.ge.198001.and.date.le.201312)
    cldc_r=f4->cldc(timeind,:,:)
    cldc=cldc_r*cldc_r@scale_factor+cldc_r@add_offset
    copy_VarCoords(cldc_r, cldc)
    delete(cldc&time)
    cldc&time=date(timeind)
    printVarSummary(cldc)
    
    ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
    opt=(/0,5/)
    ;;先插值lon方向的
    cldc=(/linmsg_n(cldc, opt, 2)/)
    ;;插值lat方向的
    opt=(/0,1/)
    cldc=(/linmsg_n(cldc, opt, 1)/)
    delete(opt)
    cldc_annual  = month_to_annual(cldc, 1)

    cldc_annual&year = ispan(1, 34, 1)
    printVarSummary(cldc_annual)
  

    ;;注意点,去掉热带整体的趋势
 
    cldc_annual_tropical=dim_avg_n_Wrap(cldc_annual(:,{30:-30},:), (/1,2/))
    ; do i=0,dimsizes(cldc_annual&year)-1
    ;     cldc_annual(i,:,:)=(/cldc_annual(i,:,:)-cldc_annual_tropical(i)/)
    ; end do
    ;;;

     ;;;region
     lat1 = 0
     lat2 = 15
     lon1 = 110
     lon2 = 155
     ;;;
   print("=====finish read data==========")

    
   ;;回归
   cldc_rc = regCoef_n(ispan(1, 34, 1), cldc_annual, 0, 0)
   copy_VarCoords(cldc_annual(1,:,:), cldc_rc)
   cldc_rc = (/cldc_rc*34/)

   ;;;

   ;;;显著性检验
    
    prob_u = rc_test(cldc_rc)
 

   ;;;
 
   ; d_rc = where(prob_d.le.0.1, d_rc, d_rc@_FillValue)
   ; sta_s_rc = where(prob_sta_s.le.0.1, sta_s_rc, sta_s_rc@_FillValue)
   ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "云量在南海区域的趋势")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlRe"

  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -80.0
  ; res@cnMaxLevelValF       =   80.0
  ; res@cnLevelSpacingF      =  10.0
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1
  res@cnMaxLevelValF       =   1
  res@cnLevelSpacingF      =  0.1
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
  
  res@mpFillOn =False

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
  plot = new(1, graphic)
  plot = gsn_csm_contour_map(wks,cldc_rc,res)

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
    i =0
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
    iplot = new(1,graphic)
  
    iplot = gsn_csm_contour(wks,prob_u,res_pro)
    overlay(plot, iplot)
    
    draw(plot)
    frame(wks)


end
    
  







    
   