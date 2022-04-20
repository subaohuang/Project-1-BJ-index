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

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

begin
  yearstart = 195001
  yearend = 201912

  lat1 = -30
  lat2 = 30
  lon1 = 0
  lon2 = 360

; ;;;;read data
    var= (/"air.mon.mean.nc","shum.mon.mean.nc","omega.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    q= q/1000.0


    time = f3->time 
    date := cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    omega = f3->omega(timeind,{1000:300},{-50:50},{lon1:lon2})
    
    omega = (/omega*-1/)

;;;

      ;;;MSE 
  plev = q&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

 
  
  ;;;读入HadiSST海温数据
  diri="/home/yangsong3/data-observation/linshh/data/sst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date:=cd_calendar(f1->time, 1)

  timeind:=ind(date.le.201912.and.date.ge.195001)

  sst  = lonFlip(f1->sst(timeind,{-50:50},:))
  
  ssta=(/rmMonAnnCycTLL(sst)/)
  copy_VarCoords(sst, ssta)
  
  ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

  ; ssta = (/runave_n_Wrap(ssta, 9*12, 0, 0)/)

 
    
;;;choose season
    season = (/"JJA","SON","DJF","MAM"/)
    plot  = new(4,graphic)
    iplot = new(4,graphic)

    plot1  = new(4,graphic)
    iplot1 = new(4,graphic)


    plot2  = new(4,graphic)
    iplot2 = new(4,graphic)


    str = "" + ispan(1,4,1) + ""
    wks = gsn_open_wks("pdf", "q-w-sst-corr.ncl")

   level = (/925,750,925,750/)

   ssta_mam = month_to_season(ssta, "MAM")

do i=0,3 

    q_col = q(:,{level(i)},:,:)
    q_col = (/rmMonAnnCycTLL(q_col)/)
    q_col = (/dtrend_msg_n(ispan(1,dimsizes(q_col&time),1),q_col,False,True,0)/)
    q_col_nep = dim_avg_n_Wrap(q_col(:,{0:20},{210:280}), (/1,2/))
    q_col_MAM  = month_to_season(q_col_nep, "MAM")
    q_col_MAM  = month_to_season(q_col_nep, "MAM")

  if i.gt.1 then 
      omega_level = omega(:,{level(i)},:,:)
      omega_level  = (/rmMonAnnCycTLL(omega_level)/)
      omega_level = (/dtrend_msg_n(ispan(1,dimsizes(omega_level&time),1),omega_level,False,True,0)/)
      omega_mam = month_to_season(omega_level, "MAM")
      rc_mam := escorc_n(q_col_MAM, omega_mam, 0, 0)
      copy_VarCoords(omega_mam(0,:,:), rc_mam)
      ;;;计算检验系数
      prob_mam := rtest(rc_mam, dimsizes(ssta_mam&time)-1, 0)
      copy_VarCoords(omega_mam(0,:,:), prob_mam)
 else 
    rc_mam := escorc_n(q_col_MAM, ssta_mam, 0, 0)
    copy_VarCoords(ssta_mam(0,:,:), rc_mam)

  ;;;计算检验系数
  prob_mam := rtest(rc_mam, dimsizes(ssta_mam&time)-1, 0)
  copy_VarCoords(ssta_mam(0,:,:), prob_mam)
 
  end if

;;;;

;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn=False

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  res@mpFillDrawOrder = "PostDraw"


  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "NCV_blu_red"
  if i.lt.2 then 
  res@gsnLeftString="Corr MAM SST &"+"MAM"+ level(i)+"hPa q"
  else
  res@gsnLeftString="Corr MAM omega &"+"MAM"+ level(i)+"hPa q"
  end if
  plot(i) = gsn_csm_contour_map(wks, rc_mam ,res)
  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 


   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 210
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

 

     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)



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

    iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)    
    ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)    
    ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)    
    ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)    

    ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)
    ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)
    ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)

    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數
     ; do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
     overlay(plot(i), iplot(i))
end do

    ;;; 
    pres = True
    pres@gsnPanelLabelBar  = True
    pres@pmLabelBarOrthogonalPosF = -0.01
    gsn_panel(wks, plot, (/2,2/), pres)

end


