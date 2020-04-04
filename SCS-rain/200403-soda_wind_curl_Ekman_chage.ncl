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
    
    ue = tauy
    ve = taux  
    ue  = (/tauy/(density*f_3d)/)   ;;;EKMAN
    ve = (/-1*taux/(density*f_3d)/)   ;;;Ekman


    ;;;calculate Horizontal and Vertical advction

    
    ; gradLatLon1  = grad_latlon_cfd (taux, taux&lat, taux&lon, True, False)
                             
    ; dxdy = gradLatLon1[0]

    ; gradLatLon2  = grad_latlon_cfd (tauy, tauy&lat, tauy&lon, True, False)
                             
    ; dydx = gradLatLon2[1]

    ; vr = dydx - dxdy

    vr = uv2vr_cfd((/taux/(density*f_3d)/) ,(/tauy/(density*f_3d)/) ,taux&lat,taux&lon,2)
     ; vr = uv2vr_cfd(taux ,tauy ,taux&lat,taux&lon,2)






    copy_VarCoords(taux, vr)
    printVarSummary(vr)
    
    wg = vr 
    copy_VarCoords(taux, wg)


    ;;;年平均
       ue_annual = month_to_annual(taux, 1)
       ve_annual = month_to_annual(tauy, 1)
       wg_annual = month_to_annual(wg, 1)
    ;;;
    

    ;;;clm diff 
    
    ue_annual_1 = dim_avg_n_Wrap(ue_annual(0:12,:,:), 0)
    ue_annual_2 = dim_avg_n_Wrap(ue_annual(13:,:,:), 0)

    ue_annual_diff = ue_annual_1 
    ue_annual_diff = (/ue_annual_2 - ue_annual_1/)

    ve_annual_1 = dim_avg_n_Wrap(ve_annual(0:12,:,:), 0)
    ve_annual_2 = dim_avg_n_Wrap(ve_annual(13:,:,:), 0)

    ve_annual_diff = ve_annual_1 
    ve_annual_diff = (/ve_annual_2 - ve_annual_1/)


    wg_annual_1 = dim_avg_n_Wrap(wg_annual(0:12,:,:), 0)
    wg_annual_2 = dim_avg_n_Wrap(wg_annual(13:,:,:), 0)

    wg_annual_diff = wg_annual_1 
    wg_annual_diff = (/wg_annual_2 - wg_annual_1/)
    ;;;;


    ;;;;clm 
    ue_clm = dim_avg_n_Wrap(ue_annual, 0)
    ve_clm = dim_avg_n_Wrap(ve_annual, 0)
    wg_clm = dim_avg_n_Wrap(wg_annual, 0)
    wg_clm  = (/wg_clm*10^5/)
    ;;;


  
 
;;;calculate correlation coef

   rc_ue = ue_annual_diff
   rc_ve = ve_annual_diff
   rc_wg = wg_annual_diff


    ; rc_ue  = (/rc_ue*10^5/)
    ; rc_ve  = (/rc_ve*10^5/)
    rc_wg  = (/rc_wg*10^5/)
   ;;;;

print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "风应力旋度和风应力的变化-SODA.3.4.2-除f")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="GMT_polar"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn= True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
  res@mpMinLatF             = 0
  res@mpMaxLatF             = 15
  res@mpMinLonF             = 110
  res@mpMaxLonF             = 155



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
  res@cnMinLevelValF       =   -0.05
  res@cnMaxLevelValF       =   0.05
  res@cnLevelSpacingF      =  0.01



  res@gsnScalarContour        = True
   ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.01           ; thin vectors
  res@vcMinMagnitudeF          =0.002
  ; res@vcMaxMagnitudeF          =0.1
  res@vcRefLengthF             = 0.02 ; define length of v

  res@vcRefMagnitudeF = 0.005
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.5

  res@vcLineArrowHeadMaxSizeF=0.005
  res@vcLineArrowHeadMinSizeF=0.005
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  plot = new(2,graphic)
  plot(0) = gsn_csm_vector_scalar_map(wks,rc_ue ,rc_ve,rc_wg ,res)

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.5
  res@cnMaxLevelValF       =   0.5
  res@cnLevelSpacingF      =  0.1
  res@vcRefMagnitudeF = 0.04
  res@vcMinDistanceF           = 0.015         ; thin vectors
  ;   res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
  plot(1) = gsn_csm_vector_scalar_map(wks,ue_clm ,ve_clm,wg_clm ,res)



  


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,3,1)

    do i = 0,1
     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
    end do


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
   str2 = "0" + ispan(1,2,1)
   do i = 0,1 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do

 
  pres = True
  gsn_panel(wks, plot, (/2,1/), pres)

 


end