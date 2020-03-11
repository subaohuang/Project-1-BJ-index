;;********************************************************
;;计算海洋动力过程项,包括平流项和垂直扩散项的趋势，将其分解成流场的变化和
;;温度的变化，区别分析
;;使用soda数据
;;
;;                                         lnshuheng 2020.03.01
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
   ;;;
 

;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
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
    
    u_annual = month_to_annual(u, 1)
    v_annual = month_to_annual(v, 1)
    wt_annual = month_to_annual(wt, 1)
    T_annual = month_to_annual(T, 1)
    
    u_clm = dim_avg_n_Wrap(u, 0)
    v_clm = dim_avg_n_Wrap(v, 0)
    wt_clm = dim_avg_n_Wrap(wt, 0)
    T_clm = dim_avg_n_Wrap(T, 0)

    u_delta =  dim_rmvmean_n_Wrap(u_annual,0)
    v_delta =  dim_rmvmean_n_Wrap(v_annual,0)
    wt_delta =  dim_rmvmean_n_Wrap(wt_annual,0)
    T_delta =  dim_rmvmean_n_Wrap(T_annual,0)


     dim = dimsizes(T_annual(1,:,:,:))
     T_mon = T
     u_mon = u
     v_mon = v
     do i = 0,yearnum-1,1
     T_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), T_annual(i,:,:,:), (/1,2,3/))
     u_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), u_annual(i,:,:,:), (/1,2,3/))
     v_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), v_annual(i,:,:,:), (/1,2,3/))
     end do 
     delete(u_annual)
     delete(v_annual)
     delete(T_annual)
     delete(wt_annual)
     
      T_eddy = T 
      T_eddy = (/T - T_mon/)     ;;;;eddy term
      u_eddy = u 
      u_eddy = (/u - u_mon/)     ;;;;eddy term
      v_eddy = v 
      v_eddy = (/v - v_mon/)     ;;;;eddy term
      delete(T_mon)
      delete(u_mon)
      delete(v_mon)
    
      u_delta&lat@units = "degrees-north"
      u_delta&lon@units = "degrees-east"
      u_eddy&lat@units = "degrees-north"
      u_eddy&lon@units = "degrees-east"

  ;;;calculate Horizontal and Vertical advction

     ;;;total 
      dT  = grad_latlon_cfd (T, T&lat, T&lon, True, False)
      hadv_total_month = u
      hadv_total_month = u*dT[1] + v*dT[0]
      hadv_total_annual = month_to_annual(hadv_total_month, 1)
      hadv_total = hadv_total_annual
      hadv_total =  dim_rmvmean_n_Wrap(hadv_total_annual,0)
      hadv_total = (/-1*hadv_total/)


     ;;;;σu*DT/dx + σv*DT/dy

      dT_clm  = grad_latlon_cfd (T_clm, T_clm&lat, T_clm&lon, True, False)
                           
      hadv_dynamic = u_delta

      hadv_dynamic = u_delta*conform_dims(dimsizes(u_delta),dT_clm[1], (/1,2,3/)) + v_delta*conform_dims(dimsizes(v_delta),dT_clm[0], (/1,2,3/))
      hadv_dynamic = (/-1*hadv_dynamic/)
     ;;;;

     ;;;;u*DσT/dx + v*DσT/dy

      dT_delta  = grad_latlon_cfd (T_delta, T_delta&lat, T_delta&lon, True, False)
      hadv_thermaldy = T_delta
      hadv_thermaldy = dT_delta[1]*conform_dims(dimsizes(dT_delta[0]),u_clm, (/1,2,3/)) + dT_delta[0]*conform_dims(dimsizes(dT_delta[1]),v_clm, (/1,2,3/))  
      hadv_thermaldy = (/-1*hadv_thermaldy/)

     ;;;;σu*DσT/dx + σv*DσT/dy
      hadv_quad  = u_delta 
      hadv_quad  = dT_delta[1]*u_delta + dT_delta[0]*v_delta
      hadv_quad = (/-1*hadv_quad/)
     ;;;;

     ;;;eddy 
      dT_eddy  = grad_latlon_cfd (T_eddy, T_eddy&lat, T_eddy&lon, True, False)
      hadv_e = u_eddy
      hadv_e = u_eddy*dT_eddy[1] + v_eddy*dT_eddy[0]
      hadv_eddy = month_to_annual(hadv_e, 1)
      hadv_eddy = (/dim_rmvmean_n(hadv_eddy, 0) /)
      hadv_eddy = (/-1*hadv_eddy/)
     ;;;
  ;;;
  
  ;;;垂直积分平均，在混合层内进行
    mlt_clm = dim_avg_n_Wrap(mlt, 0)

    z = u&depth
    z&depth = z
    print(z)    
    dz = dpres_plevel(z, dim_avg_n_Wrap(mlt, 0), 0, 0)
    copy_VarCoords(u(1,:,:,:), dz)
   
     
    dz_4d = conform_dims(dimsizes(u_delta), dz, (/1,2,3/))

    hadv_dynamic_mean = dim_sum_n_Wrap(hadv_dynamic*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    hadv_thermaldy_mean  = dim_sum_n_Wrap(hadv_thermaldy*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    hadv_quad_mean  = dim_sum_n_Wrap(hadv_quad*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    hadv_eddy_mean  = dim_sum_n_Wrap(hadv_eddy*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
    hadv_total_mean = dim_sum_n_Wrap(hadv_total*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
   
   
    copy_VarCoords(u_delta(:,1,:,:), hadv_dynamic_mean)
    copy_VarCoords(u_delta(:,1,:,:), hadv_thermaldy_mean)
    copy_VarCoords(u_delta(:,1,:,:), hadv_quad_mean)
    copy_VarCoords(u_delta(:,1,:,:), hadv_eddy_mean)
    copy_VarCoords(u_delta(:,1,:,:), hadv_total_mean)

      ;;;convert to units degree ℃/year
      hadv_dynamic_mean =365*86400*hadv_dynamic_mean
      hadv_thermaldy_mean =365*86400*hadv_thermaldy_mean
      hadv_quad_mean = 365*86400*hadv_quad_mean
      hadv_eddy_mean  = 365*86400*hadv_eddy_mean
      hadv_total_mean  = 365*86400*hadv_total_mean
    ;;;

;;;


;;;;calculate trend 
print("*******计算趋势*****************")
;;;calculate correlation coef
    
     ;;;total 
    
    r_hadv_total =   regCoef_n(ispan(1,yearnum,1), hadv_total_mean, 0, 0)  
    copy_VarCoords(hadv_dynamic_mean(0,:,:), r_hadv_total)


    ;;;hadv dynamic
    
    r_hadv_dynamic =   regCoef_n(ispan(1,yearnum,1), hadv_dynamic_mean, 0, 0)  
    copy_VarCoords(hadv_dynamic_mean(0,:,:), r_hadv_dynamic)

    ;;;hadv_thermaldy
    r_hadv_thermaldy =   regCoef_n(ispan(1,yearnum,1), hadv_thermaldy_mean, 0, 0)    
    copy_VarCoords(hadv_thermaldy_mean(0,:,:), r_hadv_thermaldy)
    ;;;

    ;;;quad 
    r_hadv_quad =   regCoef_n(ispan(1,yearnum,1), hadv_quad_mean, 0, 0)    
    copy_VarCoords(hadv_quad_mean(0,:,:), r_hadv_quad)
    
    ;;;;eddy 
    r_hadv_eddy =   regCoef_n(ispan(1,yearnum,1), hadv_eddy_mean, 0, 0)    
    copy_VarCoords(hadv_eddy_mean(0,:,:), r_hadv_eddy)
    ;;;
   
    r_hadv_dynamic = (/r_hadv_dynamic*yearnum/)
    r_hadv_thermaldy = (/r_hadv_thermaldy*yearnum/)
    r_hadv_quad = (/r_hadv_quad*yearnum/)
    r_hadv_eddy = (/r_hadv_eddy*yearnum/)
    r_hadv_total = (/r_hadv_total*yearnum/)


;;;;

print("*******显著性检验*****************")

;;;t test 
r_hadv_dynamic&lat@units = "degrees-north"
r_hadv_dynamic&lon@units = "degrees-east"
  prob_p1= rc_test(r_hadv_dynamic)
  copy_VarCoords(r_hadv_dynamic, prob_p1)

  printVarSummary(prob_p1)
  prob_p2= rc_test(r_hadv_thermaldy)

  copy_VarCoords(r_hadv_dynamic, prob_p2)

   prob_p3= rc_test(r_hadv_quad)
  copy_VarCoords(r_hadv_dynamic, prob_p3)
  
  prob_p4= rc_test(r_hadv_eddy)
  copy_VarCoords(r_hadv_dynamic, prob_p4)


  prob_p5= rc_test(r_hadv_total)
  copy_VarCoords(r_hadv_dynamic, prob_p5)




print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "海洋的平流项分解的动力与热力各项的趋势")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
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
  res@cnMinLevelValF       =   -1
  res@cnMaxLevelValF       =   1
  res@cnLevelSpacingF      =  0.1


  plot = new(5,graphic)
  res@gsnLeftString = "hadv_dynamic"
  plot(0) = gsn_csm_contour_map(wks, r_hadv_dynamic, res)
 
   res@gsnLeftString = "hadv_thermaldy"
  plot(1) = gsn_csm_contour_map(wks, r_hadv_thermaldy, res)
  
  res@gsnLeftString = "hadv_quad"
  plot(2) = gsn_csm_contour_map(wks, r_hadv_quad, res)
  
  res@gsnLeftString = "hadv_eddy"
  plot(3) = gsn_csm_contour_map(wks, r_hadv_eddy, res)
  
  res@gsnLeftString = "hadv_total"
  plot(4) = gsn_csm_contour_map(wks, r_hadv_total, res)
  

  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,5,1)

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
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   str2 = "0" + ispan(1,5,1)
   do i = 0,4
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
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(5,graphic)

    
      iplot(0) = gsn_csm_contour(wks,prob_p1,res_pro)
      iplot(1) = gsn_csm_contour(wks,prob_p2,res_pro)
      iplot(2) = gsn_csm_contour(wks,prob_p3,res_pro)
      iplot(3) = gsn_csm_contour(wks,prob_p4,res_pro)
      iplot(4) = gsn_csm_contour(wks,prob_p5,res_pro)
      

  

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
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
  pres@gsnPanelLabelBar  = True

  gsn_panel(wks, plot, (/3,2/), pres)
 


end