;;********************************************************
;;计算EKMAN输运和风应力旋度的趋势
;;数据来源soda3.3.1
;;
;;                                         lnshuheng 2020.02.18
;**********************************************************
function chazhi(data)
 begin
   xi = data&lon
   yi = data&lat
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
   data_chazhi&lat = yo
   data_chazhi&lon = xo
   copy_VarAtts(data, data_chazhi)
 
   return data_chazhi
end 
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
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160

;;;read soda data and then calculate EKMAN transport and curl of wind stress
    
       
    f7=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_taux_mn_1980-2015_chazhi.nc", "r")   
    f8=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_tauy_mn_1980-2015_chazhi.nc", "r")   


    indyear = (2013-1980+1)*12-1
   

    taux_r = f7->taux(0:indyear,{lat1:lat2},{lon1:lon2})
    tauy_r = f8->tauy(0:indyear,{lat1:lat2},{lon1:lon2})

    taux = tofloat(taux_r)
    tauy = tofloat(tauy_r)
    
   copy_VarCoords(taux_r, taux) 
   copy_VarCoords(tauy_r, tauy)
   
    delete(taux_r)
    delete(tauy_r)
    
    
    ; taux = tofloat(taux_r)
    ; tauy = tofloat(tauy_r)
    ; copy_VarCoords(taux_r, taux)
    ; copy_VarCoords(tauy_r, tauy)

    taux&lat@units = "degrees_north"
    taux&lon@units = "degrees_east"
    tauy&lat@units = "degrees_north"
    tauy&lon@units = "degrees_east"


    
    lat = taux&lat
    pi = get_pi("float")
    omega = 2*pi/86400
    f = 2*omega*sin(pi*lat/180)
    density = 1025
    f_3d = conform_dims(dimsizes(taux), f, 1)

    taux = (/taux/(density*f_3d)/)
    tauy = (/tauy/(density*f_3d)/)

    ;;parameter 
     dx = lat
     dx = 6378388.*cos(0.0174533*lat)
     dx!0 = "lat"
     copy_VarCoords(taux(1,:,1), dx)

     dy = (pi/180)*1*6378388.
     

    ;;;calculate Horizontal and Vertical advction

    
    ; gradLatLon1  = grad_latlon_cfd (taux, taux&lat, taux&lon, True, False)
                             
    ; dxdy = gradLatLon1[0]

    ; gradLatLon2  = grad_latlon_cfd (tauy, tauy&lat, tauy&lon, True, False)
                             
    ; dydx = gradLatLon2[1]

    ; vr = dydx - dxdy

    vr = uv2vr_cfd(taux,tauy,taux&lat,taux&lon,2)






    copy_VarCoords(taux, vr)
    printVarSummary(vr)
    
    wg = vr 
    copy_VarCoords(taux, wg)

     ue = tauy              ;;;Ekman volume transport
     ve = taux              ;;;Ekman volume transport
     ; ve = -1*ve
     
    ue = taux*density*f_3d
    ve = tauy*density*f_3d

    ;;;年平均
       ue_annual = month_to_annual(ue, 1)
       ve_annual = month_to_annual(ve, 1)
       wg_annnual = month_to_annual(wg, 1)
    ;;;


  
 
print("*******计算趋势*****************")
;;;calculate correlation coef

    ;;;hadv correlation with sst
    rc_ue =   regCoef_n(ispan(1,yearnum,1), ue_annual, 0, 0)  
    copy_VarCoords(ue(0,:,:), rc_ue)
    
    ;;;hadv correlation with sst
    rc_ve =   regCoef_n(ispan(1,yearnum,1), ve_annual, 0, 0)    
    copy_VarCoords(ve(0,:,:), rc_ve)
    
    rc_wg =   regCoef_n(ispan(1,yearnum,1), wg_annnual, 0, 0)    
    copy_VarCoords(wg(0,:,:), rc_wg)
  
    ;;;
   
    rc_ue = (/rc_ue*yearnum/)
    rc_ve = (/rc_ve*yearnum/)
    rc_wg = (/rc_wg*yearnum/)

    ;;;


    ; rc_ue  = (/rc_ue*10^5/)
    ; rc_ve  = (/rc_ve*10^5/)
    rc_wg  = (/rc_wg*10^5/)
   ;;;;


  
       


;;;;

print("*******显著性检验*****************")

;;;t test 
  prob_p1= rc_test(rc_ue)
  copy_VarCoords(rc_ue, prob_p1)

  prob_p2= rc_test(rc_ve)
  copy_VarCoords(rc_ve, prob_p2)
  prob_p3= rc_test(rc_wg)
  copy_VarCoords(rc_wg, prob_p3)


  ; prob_p1 = where(prob_p1.eq.0, 2, prob_p1)
  ; prob_p2 = where(prob_p2.eq.0, 2, prob_p2)

  ; prob_p3 = where(prob_p3.eq.0, 2, prob_p3)
  ; prob_p4 = where(prob_p4.eq.0, 2, prob_p4)

;;;;

 ; r_sst_hadv = where(prob_p1.le.0.5,r_sst_hadv, r_sst_hadv@_FillValue)
 ; r_dsst_hadv = where(prob_p2.le.0.5,r_dsst_hadv, r_dsst_hadv@_FillValue)
 ; r_sst_vadv = where(prob_p3.le.0.5,r_sst_vadv, r_sst_vadv@_FillValue)
 ; r_dsst_vadv = where(prob_p4.le.0.5,r_dsst_vadv, r_dsst_vadv@_FillValue)


print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "风应力旋度和风应力的趋势-SODA.3.4.2")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn= True
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
  res@cnMinLevelValF       =   -0.1
  res@cnMaxLevelValF       =   0.1
  res@cnLevelSpacingF      =  0.01
 
  res@gsnScalarContour        = True
   ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.01           ; thin vectors
  ; res@vcMinMagnitudeF          =0.01
  ; res@vcMaxMagnitudeF          =0.1
  res@vcRefLengthF             = 0.02 ; define length of v

  res@vcRefMagnitudeF = 0.008
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.5

  res@vcLineArrowHeadMaxSizeF=0.005
  res@vcLineArrowHeadMinSizeF=0.005
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  
  plot = gsn_csm_vector_scalar_map(wks,rc_ue ,rc_ve,rc_wg ,res)

  


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,3,1)

   i = 0
     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   


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
   str2 = "0" + ispan(1,3,1)
   
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
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

      iplot=new(3,graphic)

    
      iplot(0) = gsn_csm_contour(wks,prob_p3,res_pro)
    

      

  

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


    i = 0
     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    

    ;;;  
     i = 0
      overlay(plot(i), iplot(i))
  draw(plot)
  frame(wks)   

 


end