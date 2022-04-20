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

function cal_season_trend(data,season)
  begin
    p_djf_r=month_to_annual(data, 1)
    dim = dimsizes(p_djf_r)
    p_rc=regCoef_n(ispan(1, dim(0), 1), p_djf_r, 0, 0)  
    copy_VarCoords(p_djf_r(1,:,:), p_rc)
    p_rc=(/p_rc*dim(0)/)
  return p_rc
end 

; function cal_season_trend(data,season)
;   begin
;     p_djf_r=month_to_season(data, season)
;     dim = dimsizes(p_djf_r)
;     p_mk=trend_manken(p_djf_r, False, 0)
;     p_rc = p_mk(1,:,:)
;     copy_VarCoords(p_djf_r(1,:,:), p_rc)
;     p_rc=(/p_rc*dim(0)/)
;     p_rc@prob = p_mk(0,:,:)
;   return p_rc
; end 
begin
    f2=addfile("/home/yangsong3/data-observation/ERA-interim/monthly/precipitation/ERAI-pre-1degree.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    lsp=f2->lsp(timeind,:,:)  ;;;units are m
    cp=f2->cp(timeind,:,:)   
    
    ;;convert units to mm/day
    p_pi = cp
    p_pi = (/lsp+cp/)
    p_pi = p_pi*10^3

    rc_mam = cal_season_trend(p_pi,"MAM")
  

    prob_mam=rc_test(rc_mam)

    
    ; prob_mam = 1- rc_mam@prob
    ; prob_jja = 1- rc_jja@prob
    ; prob_son = 1- rc_son@prob
    ; prob_djf = 1- rc_djf@prob


    copy_VarCoords(p_pi(1,:,:), prob_mam)
   

    
wks   = gsn_open_wks("pdf", "ERA-降水的趋势-scs-1979-2013")
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
  res@cnFillColors = (/0,1,2,3,4,5,6,6,7,8,9,10,11,12/)

  res@gsnLeftString="ERA Precipitation trend"
  
  plot = new(1,graphic)
  plot(0) = gsn_csm_contour_map(wks, rc_mam, res)
  
  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 
  ;;;添加赤道直线

   pgx=ispan(30,390,10)
   pgy=new(dimsizes(pgx),typeof(pgx))
   pgy=0

   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1
   
   str = ispan(1, 4, 1) + ""
   
   plot@$str(0)$ = gsn_add_polyline(wks, plot(0), pgx, pgy, pgres)
   



  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
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

     str2 = ispan(21, 24, 1) +""
     i = 0 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  


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

    
 
    iplot(0) = gsn_csm_contour(wks,prob_mam,res_pro)    
    
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

   
     i = 0
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
     overlay(plot(i), iplot(i))
  

    ;;;  
 draw(plot)
 frame(wks)
end