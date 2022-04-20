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
  ; yearstart = 200001
  ; yearend = 201912
    yearstart = 195001
  yearend = 201912

  lat1 = -30
  lat2 = 30
  lon1 = 0
  lon2 = 360


 f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")
 date=cd_calendar(f1->time, 1)
 timeind=ind(date.ge.yearstart.and.date.le.yearend)
 MSE = f1->MSE(timeind,:,:,:)

  ;;;MSE 
  plev = MSE&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

  MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

  MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

  MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

  ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)

  ; MSE_col = dim_avg_n_Wrap(MSE_col, (/1,2/))
  MSE_col = (/MSE_col/10^6/)

  MSE_col_annual = month_to_annual(MSE_col, 1)
  MSE_col_annual!0 = "time"

  ; MSE_col_MAM = month_to_season(MSE_col, "MAM")
  ; MSE_col_JJA = month_to_season(MSE_col, "JJA")
  ; MSE_col_SON = month_to_season(MSE_col, "SON")

 
  MSE_col_tendency = MSE_col
  MSE_col_tendency = center_finite_diff_n(MSE_col, 1, False, 1, 0)

  MSE_col_tendency_JFM = month_to_season(MSE_col_tendency, "DJF")

  MSE_nep_JFM = dim_avg_n_Wrap(MSE_col_tendency_JFM(:,{0:20},{240:280}), (/1,2/))
  MSE_nep_JFM = (/dim_standardize_n(MSE_nep_JFM, 1, 0)/)


  season = (/"MAM","JJA","SON","DJF"/)
  plot = new(4,graphic)
  iplot=new(4,graphic)

wks   = gsn_open_wks("pdf", "MSE异常回归到JFM dMSE-"+yearstart+"-"+yearend)
  ; wks   = gsn_open_wks("pdf", "DMSE异常回归到JFM MSE")

do i=0,3
   MSE_col_MAM = month_to_season(MSE_col, season(i))

   if season(i).eq."DJF" then
    MSE_col_MAM(0:dimsizes(MSE_col_MAM&time)-2,:,:)  = (/MSE_col_MAM(1:dimsizes(MSE_col_MAM&time)-1,:,:)/)
    MSE_col_MAM(dimsizes(MSE_col_MAM&time)-1,:,:)  = (/MSE_col(dimsizes(MSE_col&time)-1,:,:)/)
   end if
    
    ; if season(i).eq."MAM" then 
    ;   rc_mam := regCoef_n(MSE_nep_JFM(0:dimsizes(MSE_nep_JFM&time)-2), MSE_col_MAM(1:,:,:), 0, 0)
    ;   copy_VarCoords(MSE_col_MAM(1,:,:), rc_mam)

    ; else 
    rc_mam := regCoef_n(MSE_nep_JFM, MSE_col_MAM, 0, 0)
    copy_VarCoords(MSE_col_MAM(0,:,:), rc_mam)
    ; end if
    prob_mam := rc_test(rc_mam)
 
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
  res@lbLabelBarOn=False
  ; res@lbLabelBarOn=False

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -40
  res@mpMaxLatF             = 40
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
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.1

; 
  ;  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; res@cnLevels := (/-8,-6,-4,-2,2,4,6,8/)*0.05

  ;   res@cnMinLevelValF       =   -0.6
  ; res@cnMaxLevelValF       =   0.6
  ; res@cnLevelSpacingF      =  0.1
  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_blu_red"

  ; res@gsnLeftString=" DJF SST regression onto MAM MSE  " 
  res@gsnLeftString=  "regression of "+ season(i) + " MSE aganist JFM dMSE"

  plot(i) = gsn_csm_contour_map(wks, rc_mam ,res)
  ; res@gsnLeftString="(b)JJA MSE trend " 
  ; plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  ; res@gsnLeftString="(c)SON MSE trend " 
  ; plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  ; res@gsnLeftString="(d)DJF MSE trend " 
  ; plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 
   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

  
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 240
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

     str2 = ispan(21, 24, 1) +""
     str3 = ispan(25, 28, 1) +""

     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

   ; lat1 = 0
   ; lat2 = -10
   ; lon1 = 270
   ; lon2 = 280

   ; pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   ; pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;  pgres@gsLineColor="red"


   ; plot@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)


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

    ; i=0
     ; do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
     overlay(plot(i), iplot(i))
  ; end do

    ;;;  
end do 
pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01
gsn_panel(wks, plot, (/2,2/), pres)

end


