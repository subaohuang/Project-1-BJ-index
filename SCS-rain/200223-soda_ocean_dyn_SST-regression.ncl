;;********************************************************
;;计算海洋动力过程项,SST倾向以及SST值.
;;的相关系数场，
;;NCEP-NCAR的精度为1.87×1.9,soda数据已经插值到1°乘1°
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
 

;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)


    ssta_annual = month_to_annual(ssta(:,{lat1:lat2},{lon1:lon2}), 1)
    printVarSummary(ssta_annual)

   
    ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
    ;;;sst 倾向值
     dsst = center_finite_diff_n(ssta, 1, False, 0, 0) 
     copy_VarCoords(ssta, dsst)
     dsst_annual = month_to_annual(dsst(:,{lat1:lat2},{lon1:lon2}), 1)
     dsst_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),dsst_annual,False,True,0)/)

    ;;;
    delete(date)
    delete(timeind)
;;;

;;get region mean as index     
   ; lat1 = 0
   ; lat2 = 15
   ; lon1 = 110
   ; lon2 = 155

   ssta_index = dim_avg_n_Wrap(ssta_annual(:,{0:15},{110:155}), (/1,2/))
   dsst_index = dim_avg_n_Wrap(dsst_annual(:,{0:15},{110:155}), (/1,2/))

   ;;;normalized 
   
   ssta_index = (/dim_standardize_n(ssta_index, 1, 0)/)
   dsst_index = (/dim_standardize_n(dsst_index, 1, 0)/)
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


    ;;;calculate Horizontal and Vertical advction

    gradLatLon  = grad_latlon_cfd (T, T&lat, T&lon, True, False)
                             
    T_grad_lat  = gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list' 
    T_grad_lon  = gradLatLon[1]
    delete(gradLatLon)
    
    hadv = u
    hadv = u*T_grad_lon + v*T_grad_lat
    hadv = (/-1*hadv/)

    
    T_grad_z = center_finite_diff_n(T, T&depth,False,0,1)
    
    vadv = wt
    vadv = wt*T_grad_z


  

    ;;;;垂直积分平均，在混合层内进行
    mlt_clm = dim_avg_n_Wrap(mlt, 0)

    z =  vadv&depth
    z&depth = z
    print(z)    
    dz = dpres_plevel(z, dim_avg_n_Wrap(mlt, 0), 0, 0)
    copy_VarCoords(u(1,:,:,:), dz)
     
    dz_4d = conform_dims(dimsizes(u), dz, (/1,2,3/))

    hadv_mean = dim_sum_n_Wrap(hadv*dz_4d,1)/sum(dz)
    vadv_mean  = dim_sum_n_Wrap(vadv*dz_4d,1)/sum(dz)
   
    copy_VarCoords(u(:,1,:,:), hadv_mean)
    copy_VarCoords(u(:,1,:,:), vadv_mean)

    ;;;年平均
       hadv_mean = (/rmMonAnnCycTLL(hadv_mean)/)
       vadv_mean = (/rmMonAnnCycTLL(vadv_mean)/)

       hadv_mean_annual = month_to_annual(hadv_mean, 1)
       vadv_mean_annual = month_to_annual(vadv_mean, 1)

       hadv_mean_annual = (/dtrend_msg_n(ispan(1,yearnum,1),hadv_mean_annual,False,True,0)/)
       vadv_mean_annual = (/dtrend_msg_n(ispan(1,yearnum,1),vadv_mean_annual,False,True,0)/)
    ;;;



  
 
print("*******计算回归系数*****************")
;;;calculate correlation coef

    ;;;hadv correlation with sst
    r_sst_hadv =   regCoef_n(ssta_index, hadv_mean_annual, 0, 0)
    r_dsst_hadv =   regCoef_n(dsst_index, hadv_mean_annual, 0, 0)
      
    copy_VarCoords(ssta_annual(0,:,:), r_sst_hadv)
    copy_VarCoords(dsst_annual(0,:,:), r_dsst_hadv)

    
  
    ;;;hadv correlation with sst
    r_sst_vadv =   regCoef_n(ssta_index, vadv_mean_annual, 0, 0)
    r_dsst_vadv =   regCoef_n(dsst_index, vadv_mean_annual, 0, 0)
      
    copy_VarCoords(ssta_annual(0,:,:), r_sst_vadv)
    copy_VarCoords(dsst_annual(0,:,:), r_dsst_vadv)
    ;;;
   
    ;;;
    r_sst_hadv = r_sst_hadv*10^11
    r_dsst_hadv = r_dsst_hadv*10^11
    r_sst_vadv = r_sst_vadv*10^11
    r_dsst_vadv = r_dsst_vadv*10^11
   
   printMinMax(r_sst_hadv, 1)
       


;;;;

print("*******显著性检验*****************")

;;;t test 
  prob_p1= rc_test(r_sst_hadv)
  prob_p2= rc_test(r_dsst_hadv)

  copy_VarCoords(r_sst_hadv, prob_p1)
  copy_VarCoords(r_dsst_hadv, prob_p2)

  printVarSummary(prob_p1)

  prob_p3= rc_test(r_sst_vadv)
  prob_p4= rc_test(r_dsst_vadv)

  copy_VarCoords(r_sst_vadv, prob_p3)
  copy_VarCoords(r_dsst_vadv, prob_p4)

  ; prob_p1 = where(prob_p1.eq.0, 2, prob_p1)
  ; prob_p2 = where(prob_p2.eq.0, 2, prob_p2)

  ; prob_p3 = where(prob_p3.eq.0, 2, prob_p3)
  ; prob_p4 = where(prob_p4.eq.0, 2, prob_p4)

;;;;


print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "海洋的平流项和垂直扩散项与DSST及SST的回归系数")

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
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.1
 

  plot = new(4,graphic)
  res@gsnLeftString = "SST&Horizontal advction"
  plot(0) = gsn_csm_contour_map(wks, r_sst_hadv, res)
  res@gsnLeftString = "DSST&Horizontal advction"
  plot(1) = gsn_csm_contour_map(wks, r_dsst_hadv, res)
  
  res@gsnLeftString = "SST&Vertical advction"
  plot(2) = gsn_csm_contour_map(wks, r_sst_vadv, res)
  res@gsnLeftString = "DSST&Vertical advction"
  plot(3) = gsn_csm_contour_map(wks, r_dsst_vadv, res)
 
  
  


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,4,1)

   do i = 0,3
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
   str2 = "0" + ispan(1,4,1)
   do i = 0,3
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

    iplot=new(4,graphic)

    
      iplot(0) = gsn_csm_contour(wks,prob_p1,res_pro)
      iplot(1) = gsn_csm_contour(wks,prob_p2,res_pro)
      iplot(2) = gsn_csm_contour(wks,prob_p3,res_pro)
      iplot(3) = gsn_csm_contour(wks,prob_p4,res_pro)
     

  

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


    do i = 0,3
     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    end do

    ;;;  
    do i = 0,3
      overlay(plot(i), iplot(i))
    end do 

  pres = True
  pres@gsnPanelLabelBar  = True

  gsn_panel(wks, plot, (/2,2/), pres)
 


end