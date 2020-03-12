;;********************************************************
;;回归SODA的u,v趋势
;;
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
begin
  
;==========================================================
;;;计算海水U,V,T的趋势
;;==========================================================
    
  
    ;;region 
      lat1 = -10
      lat2 = 40
      lon1 = 90
      lon2 = 160
    ;;;

    ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;

    f2=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_mlt_mn_1980-2015_chazhi.nc", "r")     
    
      
    indyear = (2013-1980+1)*12-1
    u = f2->u(0:indyear,:,{lat1:lat2},{lon1:lon2})
    v = f3->v(0:indyear,:,{lat1:lat2},{lon1:lon2})
    wt = f4->wt(0:indyear,:,{lat1:lat2},{lon1:lon2})
    T = f5->temp(0:indyear,:,{lat1:lat2},{lon1:lon2})
    mlt = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})
    
    u&lat@units = "degree_north"
    v&lat@units = "degree_north"
    T&lat@units = "degree_north"

    u&lon@units = "degree_east"
    v&lon@units = "degree_east"
    T&lon@units = "degree_east"
    
  ;;;混合层内平均
      ;;;;垂直积分平均，在混合层内进行
    mlt_clm = dim_avg_n_Wrap(mlt, 0)

    z =  u&depth
    z&depth = z
    print(z)    
    dz = dpres_plevel(z, dim_avg_n_Wrap(mlt, 0), 0, 0)
    copy_VarCoords(u(1,:,:,:), dz)
     
    dz_4d = conform_dims(dimsizes(u), dz, (/1,2,3/))

    u_mean = dim_sum_n_Wrap(u*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    v_mean  = dim_sum_n_Wrap(v*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    T_mean = dim_sum_n_Wrap(T*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    copy_VarCoords(u(:,1,:,:), u_mean)
    copy_VarCoords(v(:,1,:,:), v_mean)
    copy_VarCoords(T(:,1,:,:), T_mean)
  ;;;

   u_annual  = month_to_annual(u_mean, 1)
   v_annual  = month_to_annual(v_mean, 1)
   T_annual  = month_to_annual(T_mean, 1)

  ;;计算趋势
    rc_u = regCoef_n(ispan(1,34, 1), u_annual, 0, 0)
    rc_v = regCoef_n(ispan(1,34, 1), v_annual, 0, 0)
    rc_T = regCoef_n(ispan(1,34, 1), T_annual, 0, 0)
  ;;;
   
   rc_u = yearnum*rc_u
   rc_v = yearnum*rc_v
   rc_T = yearnum*rc_T
   copy_VarCoords(u(1,1,:,:), rc_u)
   copy_VarCoords(v(1,1,:,:), rc_v)
   copy_VarCoords(T(1,1,:,:), rc_T)
   
  ;;;

  prob_T = rc_test(rc_T)
  prob_u = rc_test(rc_u)
  prob_v = rc_test(rc_v)
  rc_T = where(prob_T.le.0.05, rc_T, rc_T@_FillValue)
  ; rc_v = where(prob_v.le.0.1, rc_v, rc_T@_FillValue)
  ; rc_u = where(prob_u.le.0.1, rc_u, rc_T@_FillValue)




;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "海洋u,v,T的趋势-SCS")
  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueWhiteOrangeRed"

  res@cnFillOn=True
  res@cnLinesOn=False
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
  
  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.01           ; thin vectors
  res@vcMinMagnitudeF          =0.01
  res@vcMaxMagnitudeF          =0.1
  res@vcRefLengthF             = 0.015 ; define length of v
  res@vcRefMagnitudeF = 0.05
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.5

  res@vcLineArrowHeadMaxSizeF=0.005
  res@vcLineArrowHeadMinSizeF=0.005
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  
  plot = gsn_csm_vector_scalar_map(wks,rc_u,rc_v,rc_T ,res)


  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

 

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
   i = 0
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   

   ; ;;;打点
   ;  res_pro=True
   ;  res_pro@gsnAddCyclic = False
   ;  res_pro@gsnDraw=False
   ;  res_pro@gsnFrame=False
   ;  res_pro@cnFillOn = True

   ;  res_pro@cnLevelSelectionMode = "ExplicitLevels"
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  res_pro@lbLabelBarOn = False
   ;  res_pro@cnLinesOn=False
   ;  ; res_pro@cnMonoLineColor =False
   ;  ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

   ;  res_pro@cnMonoFillPattern = True
   ;  res_pro@cnFillPattern = 17
   ;  res_pro@cnFillDotSizeF = 0.0045
   ;  res_pro@cnMonoFillColor = False
   ;  res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
   ;  ; res_pro@cnLevels = ispan(-6,6,1)
   ;  res_pro@cnLineLabelsOn=False
   ;  res_pro@tmYRMode = "Automatic"
   ;  res_pro@cnInfoLabelOn=False
   ;  ; res_pro@cnLineLabelBackgroundColor="Transparent"
   ;  ; res_pro@cnLineThicknessF = 3.
   ;  ; res_pro@cnFillOn        = True
   ;  iplot = new(1,graphic)

   ;  iplot(0)=gsn_csm_contour(wks,prob_T,res_pro)
   

   ;  i = 0
   ;  overlay(plot(i), iplot(i))
    
    draw(plot)
    frame(wks)

end
    
  







    
   