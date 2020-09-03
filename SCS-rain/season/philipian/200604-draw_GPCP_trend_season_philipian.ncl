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
  p_djf_r=month_to_season(data, season)
   if season.eq."DJF"
      p_djf_r(0:dimsizes(p_djf_r&time)-2,:,:) = p_djf_r(1:,:,:)
      p_djf_r(dimsizes(p_djf_r&time)-1,:,:) = data(dimsizes(data&time)-1,:,:)
    end if
    p_djf_r=month_to_season(data, season)
    dim = dimsizes(p_djf_r)
    p_rc=regCoef_n(ispan(1, dim(0), 1), p_djf_r, 0, 0)
    copy_VarCoords(p_djf_r(1,:,:), p_rc)
    p_rc=(/p_rc*dim(0)/)
  return p_rc
end 
; function cal_season_trend(data,season)
;   begin
;     p_djf_r=month_to_season(data, season)

;     if season.eq."DJF"
;       p_djf_r(0:dimsizes(p_djf_r&time)-2,:,:) = p_djf_r(1:,:,:)
;       p_djf_r(dimsizes(p_djf_r&time)-1,:,:) = data(dimsizes(data&time)-1,:,:)
;     end if
;     dim = dimsizes(p_djf_r)
;     p_mk=trend_manken(p_djf_r, False, 0)
;     p_rc = p_mk(1,:,:)
;     copy_VarCoords(p_djf_r(1,:,:), p_rc)
;     p_rc=(/p_rc*dim(0)/)
;     p_rc@prob = p_mk(0,:,:)
;   return p_rc
; end 

begin
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201512)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    p_pi = p
    
    rc_mam = cal_season_trend(p_pi,"MAM")
    rc_jja = cal_season_trend(p_pi,"JJA")
    rc_son = cal_season_trend(p_pi,"SON")
    rc_djf = cal_season_trend(p_pi,"DJF")

    prob_mam=rc_test(rc_mam)
    prob_jja=rc_test(rc_jja)
    prob_son=rc_test(rc_son)
    prob_djf=rc_test(rc_djf)
    
    ; prob_mam = 1- rc_mam@prob
    ; prob_jja = 1- rc_jja@prob
    ; prob_son = 1- rc_son@prob
    ; prob_djf = 1- rc_djf@prob


    copy_VarCoords(p_pi(1,:,:), prob_mam)
    copy_VarCoords(p_pi(1,:,:), prob_jja)
    copy_VarCoords(p_pi(1,:,:), prob_son)
    copy_VarCoords(p_pi(1,:,:), prob_djf)

    
wks   = gsn_open_wks("pdf", "不同季节GPCP降水的趋势-phipipian-1979-2015")
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="NCV_blu_red"
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  res@mpMaxLatF             = 30
  res@mpMinLonF             =100
  res@mpMaxLonF             = 170



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

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.7

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -2.5
  res@cnMaxLevelValF       =   2.5
  res@cnLevelSpacingF      =  0.25

  res@cnFillPalette="MPL_BrBG"
  res@gsnLeftString="Spring"

  plot = new(4,graphic)
  plot(0) = gsn_csm_contour_map(wks, rc_mam, res)
  res@gsnLeftString="Summer"
  plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  res@gsnLeftString="Autumn"
  plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  res@gsnLeftString="Winter"
  plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


  

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
   do i =0,3 
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   end do 



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
   do i =0,3 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 


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

    iplot=new(4,graphic)

    
 
    iplot(0) = gsn_csm_contour(wks,prob_mam,res_pro)    
    iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)
    iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)
    iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)

    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

   
     do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
     overlay(plot(i), iplot(i))
  end do

    ;;;  
   
pres = True
pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,2/), pres)

end