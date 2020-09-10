;;********************************************************
;;画HadIsst的sst趋势
;;
;;lnshuheng 2020.02.17
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
;;********************************************************
;;计算速度势函数以及辐散风场
;; ERA-interim:u,v，
;; used ncl function: uv2dv_cfd,
;；；grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2020.09.04
;**********************************************************
begin
    
;;;time and region
 timestart = 198001
 timeend = 201512
;;;

lat1 = -30 
lat2 = 30
lon1 = 70
lon2 = 170

;;;;read data
    diri="/home/yangsong3/data-observation/ERA-interim/monthly/"
    fils1=systemfunc("ls "+diri+"uwnd.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"vwnd.mon.mean.nc")

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)

    timeind=ind(date.le.timeend.and.date.ge.timestart)
    timeind2=ind(date2.le.timeend.and.date2.ge.timestart)        


    u = short2flt(f1->u(timeind,{1000:100},::-1,:))
    v = short2flt(f2->v(timeind,{1000:100},::-1,:))

    u_850 = u(:,{850},{lat1:lat2},{lon1:lon2})
    v_850 = v(:,{850},{lat1:lat2},{lon1:lon2})

    u_200 = u(:,{200},{lat1:lat2},{lon1:lon2})
    v_200 = v(:,{200},{lat1:lat2},{lon1:lon2})
    


    ; u_850_annual = month_to_annual(u_850, 1)
    ; u_200_annual = month_to_annual(u_200, 1)
    
    ; v_850_annual = month_to_annual(v_850, 1)
    ; v_200_annual = month_to_annual(v_200, 1)

    u_850_annual = month_to_season(u_850, "MAM")
    u_200_annual = month_to_season(u_200, "MAM")
    
    v_850_annual = month_to_season(v_850, "MAM")
    v_200_annual = month_to_season(v_200, "MAM")

    ;;;风切变情况
    stear = sqrt((u_200_annual- u_850_annual)^2 + (v_200_annual- v_850_annual)^2)
    copy_VarCoords(u_200_annual, stear)
    ;;;
    

   dim = dimsizes(stear)

    ;;;计算趋势
    stear_rc=regCoef_n(ispan(1, dim(0), 1), stear, 0, 0)
    copy_VarCoords(stear(1,:,:), stear_rc)


    prob_p=rc_test(stear_rc)

;;;;

;;;; 
   wks   = gsn_open_wks("pdf", "stear趋势-SCS-1980-2015-MAM")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(1,graphic)
  res=True
  res@gsnAddCyclic = False
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  
  ; res@cnFillPalette = "NCV_jaisnd"
 res@cnFillPalette = "BlueWhiteOrangeRed"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2

  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             =1
  ; res@mpMaxLonF             = 155



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

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.06
  res@cnMaxLevelValF       =   0.06
  res@cnLevelSpacingF      =  0.01
  
  ; ; res@cnMinLevelValF       =   0.2
  ; ; res@cnMaxLevelValF       =   0.6
  ; ; res@cnLevelSpacingF      =  0.05

 
  ; res@cnMinLevelValF       =   -0.2
  ; res@cnMaxLevelValF       =   0.2
  ; res@cnLevelSpacingF      =  0.025

  ; res@gsnLeftString="HadISST sst trend"
  res@gsnLeftString="SST trend (remove regional average)"
  plot = gsn_csm_contour_map(wks, stear_rc, res)

  ; res2 = res
  ; res2@cnLevelSelectionMode = "ExplicitLevels"
  ; res2@cnFillOn = False
  ; res2@cnLinesOn = True
  ; res2@cnLevels = (/0.55/)
  ; res2@cnLineColors = "black"
  ; res2@cnMonoLineColor  = True
  ; res2@cnLineThicknessF = 5
  ; res2@cnMonoLineThickness  = True
  ; res2@gsnLeftString = ""
  ; res2@cnInfoLabelOn = False

  ; plot1 = gsn_csm_contour(wks, stear_rc, res2)
  ; overlay(plot, plot1)  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    


  ;;;添加赤道直线
   pgx=ispan(30,390,10)
   pgy=new(dimsizes(pgx),typeof(pgx))
   pgy=0
   str1=unique_string("dum1")
   str2=unique_string("dum2")
    str3=unique_string("dum3")
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1
   plot@$str1$ = gsn_add_polyline(wks, plot(0), pgx, pgy, pgres)
  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
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
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
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
  
  
     ; iplot = gsn_csm_contour(wks,prob_p  ,res_pro)

    ;;;  
    ; overlay(plot, iplot)
draw(plot)
frame(wks)

end 


