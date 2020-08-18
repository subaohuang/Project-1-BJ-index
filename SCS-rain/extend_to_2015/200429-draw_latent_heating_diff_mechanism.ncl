;;;画ERA 的0-12小时forecast月平均的潜热通量数据的牛顿冷却，风蒸发，
;;稳定度对潜热趋势的贡献的图
;;
;;
;;
;;                                         lnshuheng 2020.04.29
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
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   ; lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160
   lat1 = -10
   lat2 = 20
   lon1 = 125
   lon2 = 160
;;;;

;;;读入风场数据
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/ERA_interim_10m_wind_speed.nc"
    f1 = addfile(path_wind,"r")

    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)

    ; u = short2flt(f1->si10(timeind,{lat1:lat2},{lon1:lon2}))
    u = short2flt(f1->si10(timeind,{lat1:lat2},{lon1:lon2}))
    u_annual = month_to_annual(u, 1)
    printVarSummary(u)
    u_clm = dim_avg_n_Wrap(u, 0)

;;;;

;;;read ERA heat flux data

    fils1 = "/home/ys17-19/lsh/Project/SCS-rain/annual/ERAI-heatflux-1degree.nc"
    f1 = addfile(fils1, "r")

    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.198001)
  
    
    lhtfl_r = f1->lhtfl(timeind,{lat1-1:lat2+1},{lon1-1:lon2+1})
    lhtfl_r = (/lhtfl_r/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    printVarSummary(lhtfl_r)
    
    print("*****************开始水平插值************************")  
    xo  =  u&longitude
    yo  =  u&latitude
    lhtfl = linint2_Wrap(lhtfl_r&lon, lhtfl_r&lat, lhtfl_r, \
                  True, xo, yo, 0)
    printVarSummary(lhtfl_r)
    delete(lhtfl_r)

    lhtfl_annual = month_to_annual(lhtfl, 1)
    lhtfl_clm = dim_avg_n_Wrap(lhtfl, 0)
;;;;

 

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils2=systemfunc("ls "+diri+"HadISST_sst.nc")
    f2=addfile(fils2, "r")
    date:=cd_calendar(f2->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.198001)

    sst_r1 = lonFlip(f2->sst)
    sst_r = sst_r1(timeind,{lat1-1:lat2+1},{lon1-1:lon2+1})
    delete(sst_r1)
    printVarSummary(sst_r)

    print("*****************开始水平插值************************")  
    xo  =  u&longitude
    yo  =  u&latitude
    sst = linint2_Wrap(sst_r&longitude, sst_r&latitude, sst_r, \
                  True, xo, yo, 0)
    printVarSummary(sst)
    delete(sst_r)

    ; ssta=(/rmMonAnnCycTLL(sst)/)
    ; copy_VarCoords(sst, ssta)
    ssta = sst
    ssta = (/ssta + 273.6/)
    ssta_annual = month_to_annual(ssta, 1)
   
    ;;;
    delete(date)
    delete(timeind)

;;;



;;;;read ERA-dew temp data 
    yearin  = ispan(1980, 2015, 1)
    fils3 = "/home/yangsong3/data-observation/ERA-interim/monthly/2m_dew_tem/2m_dew_temp_monthly_" + yearin + ".nc"
    f3 = addfiles(fils3, "r")

    date:=cd_calendar(f3[:]->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude({lat1:lat2})
    lon = f3[1]->longitude({lon1:lon2})
  
    d2m = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
   
    do i = 0,dimsizes(timeind)/12-1
      d2m(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->d2m(:,{lat1:lat2},{lon1:lon2}))
    end do

   ;;;计算相对湿度
    rh =  relhum_ttd(sst+273.6,d2m,1)
    copy_VarMeta(d2m, rh)
    printVarSummary(rh)
    rh_annual = month_to_annual(rh, 1)
    rh_clm = dim_avg_n_Wrap(rh_annual, 0)

    
   ;;;
;;;


;;;;读入2-m温度数据
    path_t = "/home/yangsong3/data-observation/ERA-interim/monthly/2_metre_temperature.mon.mean.nc"
    f1 = addfile(path_t,"r")

    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)

    t = short2flt(f1->t2m(timeind,{lat1:lat2},{lon1:lon2}))

    t_diff  = t 
    t_diff = t - sst -273.6

    t_diff_annual  = month_to_annual(t_diff, 1)
    t_diff_clm = dim_avg_n_Wrap(t_diff_annual, 0)

print("*******开始计算趋势系数*****************")
;;;calculate correlation coef
      
;;;
    ;;;潜热
    r_lh =  regCoef_n(ispan(1, yearnum, 1), lhtfl_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_lh)
    r_lh = (/r_lh*yearnum/)
    ;;;

    ;;;sst
    r_sst =  regCoef_n(ispan(1, yearnum, 1), ssta_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst)
    r_sst = (/r_sst*yearnum/)


    ;;;风速
    r_w =  regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_w)
    r_w = (/r_w*yearnum/)

    ;;大气-海洋温差
    r_t =  regCoef_n(ispan(1, yearnum, 1), t_diff_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_t)
    r_t = (/r_t*yearnum/)
    ;;;

 
    ;;;相对湿度的变化
    r_rh =  regCoef_n(ispan(1, yearnum, 1), rh_annual, 0, 0)
    printVarSummary(r_rh)

    copy_VarCoords(ssta_annual(0,:,:), r_rh)
    r_rh = (/r_rh*yearnum/)
    ;;;
;;;

;;;计算不同物理机制对潜热变化的贡献
  alfa = 0.06 ;;;k-1

  ;;;牛顿冷却
   Qo = r_sst
   Qo = (/alfa*lhtfl_clm*r_sst/)

  ;;;

  ;;;风速的效应
   Qw = u_clm
   Qw = lhtfl_clm*r_w/(u_clm)
   

  ;;;

  ;;;相对湿度效应
    alfa2 = exp(-1*alfa*t_diff_clm)  ;;; alfa2 = exp(beta(T)/beta(T+S)) ;;beta(T) = alfa*T
    Qrh = rh_clm
    Qrh  = lhtfl_clm*r_rh/(rh_clm - alfa2)
   ;;;

   ;;;稳定度(海表温差的作用)
    Qs = rh_clm
    Qs = alfa*lhtfl_clm*rh_clm*r_t/(rh_clm - alfa2)
   ;;;
;;;;

; Qo = (/Qo + Qw + Qrh + Qs/)

; ;;;;

print("*******显著性检验*****************")

;;;t test 
   prob_lh= rc_test(r_lh)
   prob_sst = rc_test(r_sst)
   prob_w = rc_test(r_w)
   prob_rh = rc_test(r_rh)
   prob_s = rc_test(r_t)


    copy_VarCoords(r_lh, prob_lh)
    copy_VarCoords(r_sst, prob_sst)
    copy_VarCoords(r_w, prob_w)
    copy_VarCoords(r_rh, prob_rh)
    copy_VarCoords(r_t, prob_s)




print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "ERA-潜热通量及各机制分解的趋势-2015")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  ; res@cnFillPalette="ViBlGrWhYeOrRe"
  res@cnFillPalette="BlueRed"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillDrawOrder = "PostDraw"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             =lon1
  res@mpMaxLonF             = lon2

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
  res@cnMinLevelValF       =   -20.0
  res@cnMaxLevelValF       =   20.0
  res@cnLevelSpacingF      =  2.5
 

  plot = new(5,graphic)

  res@gsnRightString = ""
  res@gsnLeftString = "Latent Heat Flux"
  plot(0) = gsn_csm_contour_map(wks, r_lh, res)
 ; res@cnFillPalette="BlRe"
  ; res@cnMinLevelValF       =   -6.0
  ; res@cnMaxLevelValF       =   -6.0
  ; res@cnLevelSpacingF      =  2.5
  ; res@gsnSpreadColorStart = 4
  ; res@gsnSpreadColorEnd = 41

  ; res@gsnSpreadColorStart = 3
  ; res@gsnSpreadColorEnd = 53


  res@gsnLeftString = "Newtonian cooling effect"
  plot(1) = gsn_csm_contour_map(wks,Qo, res)
  ;   res@cnFillPalette="ViBlGrWhYeOrRe"
  ; res@cnMinLevelValF       =   -20.0
  ; res@cnMaxLevelValF       =   20.0
  ; res@cnLevelSpacingF      =  2.5
  ; delete(res@gsnSpreadColorStart)
  ; delete(res@gsnSpreadColorEnd)
  res@gsnLeftString = "WES effect"
  plot(2) = gsn_csm_contour_map(wks, Qw, res)
  
  res@gsnLeftString = "Relative humidity effect"
  plot(3) = gsn_csm_contour_map(wks, Qrh, res)

  res@gsnLeftString = "Stability effect"
  plot(4) = gsn_csm_contour_map(wks, Qs, res)
  
  

  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,dimsizes(plot),1)

   do i = 0,4
     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   end do


  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ; ;;;画box
   ; lat1 = 0
   ; lat2 = 15
   ; lon1 = 110
   ; lon2 = 155

   ; pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   ; pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ; ;;;
   ; str2 = "0" + ispan(1,dimsizes(plot),1)
   ; do i = 0,dimsizes(plot)-1
   ;   plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   ; end do

    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1      ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(5,graphic)

    iplot(0) = gsn_csm_contour(wks,prob_lh({0:15},{125:160}),res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_sst({0:15},{125:160}),res_pro)
    iplot(2) = gsn_csm_contour(wks,prob_w({0:15},{125:160}),res_pro)
    iplot(3) = gsn_csm_contour(wks,prob_rh({0:15},{125:160}),res_pro)
    iplot(4) = gsn_csm_contour(wks,prob_s({0:15},{125:160}),res_pro)

 
    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    opt@gsnShadeFillDotSizeF =  0.005
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


    do i = 0,4
     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    end do

    ;;;  
    do i = 0,4
      overlay(plot(i), iplot(i))
    end do 

  pres = True
  pres@gsnPanelLabelBar  = False
  pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h"/)

  gsn_panel(wks, plot, (/3,2/), pres)
 


end