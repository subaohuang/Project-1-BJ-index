;;;画ERA 的0-12小时forecast月平均的通量数据的趋势
;;
;;
;;                                         lnshuheng 2020.03.10
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
    endmon = 201312

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
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160
;;;;

 

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
    ssta_annual =(/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)

    ;;;sst 倾向值
     dsst = center_finite_diff_n(ssta, 1, False, 0, 0) 
     copy_VarCoords(ssta, dsst)
     dsst_annual = month_to_annual(dsst(:,{lat1:lat2},{lon1:lon2}), 1)
     dsst_annual =(/dtrend_msg_n(ispan(1,yearnum,1),dsst_annual,False,True,0)/)

    ;;;
    delete(date)
    delete(timeind)
;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
     
    f2=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")
    f6=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_mlt_mn_1980-2015_chazhi.nc", "r") 
    


    indyear = (2013-1980+1)*12-1
    u_r = f2->u(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    v_r = f3->v(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    wt_r = f4->wt(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    T_r = f5->temp(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    mlt_r = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})
    
    u = tofloat(u_r)
    v = tofloat(v_r)
    wt = tofloat(wt_r)
    T = tofloat(T_r)
    mlt = tofloat(mlt_r)

   copy_VarCoords(u_r, u) 
   copy_VarCoords(v_r, v)
   copy_VarCoords(wt_r, wt)
   copy_VarCoords(T_r, T)
   copy_VarCoords(mlt_r, mlt)

    delete(u_r)
    delete(v_r)
    delete(wt_r)
    delete(mlt_r)
    delete(T_r)


 
    ;;;calculate Horizontal and Vertical advction

     
    ;;;;垂直积分平均，在混合层内进行
 

 
     z = u&depth 
    ; dz_4d = conform_dims(dimsizes(u), dz, (/1,2,3/))
      dz_4d = new(dimsizes(u),float)
      do i = 0,dimsizes(mlt&time)-1
        dz_4d(i,:,:,:) = dpres_plevel(z, mlt(i,:,:), 0, 0)
      end do 
    ;;;1混合层平均的
     T_mean = T(:,1,:,:)
     u_mean = u(:,1,:,:)
     v_mean = v(:,1,:,:)

     T_mean = dim_sum_n_Wrap(T*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
     u_mean = dim_sum_n_Wrap(u*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)
     v_mean = dim_sum_n_Wrap(v*dz_4d,1)/dim_sum_n_Wrap(dz_4d, 1)

    ;;;

    gradLatLon  = grad_latlon_cfd (T_mean, T_mean&lat, T_mean&lon, True, False)
                             
    T_grad_lat  = gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list' 
    T_grad_lon  = gradLatLon[1]
    delete(gradLatLon)
    
    hadv = u_mean
    hadv = u_mean*T_grad_lon + v_mean*T_grad_lat
    hadv = (/-1*hadv/)

    ;;;垂直项


     ;;;将wt进行垂直方向差值，为了得到更加精确的结果
        

      ;;;mlt变化引起温度的变化，可以不考虑
          x1 = 0
          x2 = 15
          y1 = 110
          y2 = 155
          mlt_scs = mlt(:,{x1:x2},{y1:y2})
          printVarSummary(mlt_scs)
          
        ;   do i = 0,dimsizes(mlt&time)-1
        ;   mlt_scs(i,:,:) = where(ismissing(mlt_scs(i,:,:)),dim_avg_n_Wrap(mlt_scs, 0), mlt_scs(i,:,:))
        ; end do 
        mlt_scs = where(ismissing(mlt_scs),avg(mlt_scs), mlt_scs) 
          
          printMinMax(mlt_scs, 1)
          dim  =dimsizes(mlt_scs)


    f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/SODA3.4.2_wt_T_in_MLD.nc", "r")
       
       wt_mlt = f1->wt_mlt
       T_mlt = f1->T_mlt


       vadv  = mlt_scs
       vadv = wt_mlt*(T_mean(:,{x1:x2},{y1:y2}) - T_mlt)/mlt_scs

    vadv = (/-1*vadv/)
    printVarSummary(vadv)

    ;;;convert to units degree ℃/year
      hadv =365*86400*hadv
      vadv =365*86400*vadv
    ;;;
    

   




    hadv_mean = hadv
    vadv_mean  = vadv
   
    copy_VarCoords(u(:,1,:,:), hadv_mean)
    copy_VarCoords(u(:,1,{x1:x2},{y1:y2}), vadv_mean)

    hadv_mean&lat@units = "degrees_north"
    hadv_mean&lon@units = "degrees_east"
    
    vadv_mean&lat@units = "degrees_north"
    vadv_mean&lon@units = "degrees_east"
 


    ;;;年平均
       ; hadv_mean = (/rmMonAnnCycTLL(hadv_mean)/)
       ; vadv_mean = (/rmMonAnnCycTLL(vadv_mean)/)

       hadv_mean_annual = month_to_annual(hadv_mean, 1)
       vadv_mean_annual = month_to_annual(vadv_mean, 1)
       total_ocean_annual = vadv_mean_annual
       total_ocean_annual = hadv_mean_annual(:,{x1:x2},{y1:y2}) + vadv_mean_annual

       total_ocean_annual&lat@units = "degrees_north"
       total_ocean_annual&lon@units = "degrees_east"
 

    ;;;

;;;read NCEP-NCAR heat flux data
    
 
    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/ERAI-heatflux-1degree.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)
  
    
    lhtfl = f3->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f3->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f3->nlwrs(timeind,{lat1:lat2},{lon1:lon2})

    lhtfl = (/lhtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    shtfl = (/shtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nswrs = (/nswrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nlwrs = (/nlwrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2


    ;;;nswrs_pen 
     R = 0.58
     L1 = 0.35
     L2 = 23
     nswrs_pen = nswrs
     nswrs_pen =  nswrs*(R*exp(-1*mlt/L1)+ (1-R)*exp(-1*mlt/L2))
    ;;;

    ;;;convert to units degree ℃/year
      density = 1025  ;;kg/m3
      cp  = 3850         ;;j/kg℃
      coef  = mlt*density*cp

      lhtfl = 365*86400*lhtfl/coef
      shtfl = 365*86400*shtfl/coef
      nswrs = 365*86400*nswrs/coef
      nlwrs = 365*86400*nlwrs/coef
      nswrs_pen = 365*86400*nswrs_pen/coef
      
    ;;;

    ; lhtfl  = (/rmMonAnnCycTLL(lhtfl)/)
    ; shtfl  = (/rmMonAnnCycTLL(shtfl)/)
    ; nswrs  = (/rmMonAnnCycTLL(nswrs)/)
    ; nlwrs  = (/rmMonAnnCycTLL(nlwrs)/)
    ; nswrs_pen  = (/rmMonAnnCycTLL(nswrs_pen)/)


    lhtfl_annual = month_to_annual(lhtfl, 1)
    shtfl_annual = month_to_annual(shtfl, 1)
    nswrs_annual = month_to_annual(nswrs, 1)
    nlwrs_annual = month_to_annual(nlwrs, 1)
    nswrs_pen_annual = month_to_annual(nswrs_pen, 1)
    total_annual = lhtfl_annual
    total_annual  = lhtfl_annual + shtfl_annual + nswrs_annual + nlwrs_annual -nswrs_pen_annual


       ; printMinMax(total_ocean_annual, 1)
       ; printMinMax(total_annual, 1)
       ; printMinMax(hadv_mean_annual, 1)
       ; printMinMax(vadv_mean_annual, 1)

 total_ocean_mean = dim_avg_n_Wrap(total_ocean_annual, 0)
 ; total_mean = dim_avg_n_Wrap(total_annual, 0)
 ; total = total_ocean_mean + total_mean

;;;




print("*******开始计算趋势系数*****************")
;;;calculate correlation coef

    ;;;大气
    r_lh =  regCoef_n(ispan(1, yearnum, 1), lhtfl_annual, 0, 0)
    r_sh   = regCoef_n(ispan(1, yearnum, 1), shtfl_annual, 0, 0)

    r_ns =  regCoef_n(ispan(1, yearnum, 1), nswrs_annual, 0, 0)
    r_nl   = regCoef_n(ispan(1, yearnum, 1), nlwrs_annual, 0, 0)
    r_total = regCoef_n(ispan(1, yearnum, 1), total_annual, 0, 0)

    copy_VarCoords(ssta_annual(0,:,:), r_lh)
    copy_VarCoords(ssta_annual(0,:,:), r_sh)
    copy_VarCoords(ssta_annual(0,:,:), r_ns)
    copy_VarCoords(ssta_annual(0,:,:), r_nl)
    copy_VarCoords(ssta_annual(0,:,:), r_total)
    ;;;

    ;;;海洋
    r_hadv =  regCoef_n(ispan(1, yearnum, 1), hadv_mean_annual, 0, 0)
    r_vadv   = regCoef_n(ispan(1, yearnum, 1), vadv_mean_annual, 0, 0)
    r_total_ocean   = regCoef_n(ispan(1, yearnum, 1), total_ocean_annual, 0, 0)    

    copy_VarCoords(hadv_mean(0,:,:), r_hadv)
    copy_VarCoords(vadv_mean(0,:,:), r_vadv)
    copy_VarCoords(vadv_mean(0,:,:), r_total_ocean)
    ;;;;;

    r_lh = (/r_lh*yearnum/)
    r_sh = (/r_sh*yearnum/)
    r_ns = (/r_ns*yearnum/)
    r_nl = (/r_nl*yearnum/)
    r_total = (/r_total*yearnum/)
    
    r_hadv = (/r_hadv*yearnum/)
    r_vadv = (/r_vadv*yearnum/)
    r_total_ocean = (/r_total_ocean*yearnum/)

    ; r_total_mean = avg(dim_avg_n_Wrap(r_total({0:15},{110:155}), 1))
    ; r_total_ocean_mean = avg(dim_avg_n_Wrap(r_total_ocean({0:15},{110:155}), 1))
; 
; print("r_total="+r_total_mean)
; print(data)t("r_total_ocean="+r_total_ocean_mean)




;;;;

print("*******显著性检验*****************")

;;;t test 
   prob_lh= rc_test(r_lh)
   prob_sh= rc_test(r_sh)
   prob_nl= rc_test(r_nl)
   prob_ns= rc_test(r_ns)
   prob_total= rc_test(r_total)

   prob_hadv= rc_test(r_hadv)
   prob_vadv= rc_test(r_vadv)
   prob_total_ocean = rc_test(r_total_ocean)
   


  copy_VarCoords(r_lh, prob_lh)
  copy_VarCoords(r_lh, prob_sh)
  copy_VarCoords(r_lh, prob_nl)
  copy_VarCoords(r_lh, prob_ns)
  copy_VarCoords(r_total, prob_total)

  copy_VarCoords(r_hadv, prob_hadv)
  copy_VarCoords(r_vadv, prob_vadv)
  copy_VarCoords(r_total_ocean, prob_total_ocean)
  
  
 ;  printVarSummary(prob_p1)
 ; prob_p1 = where(prob_p1.eq.0, 2, prob_p1)
 ; prob_p2 = where(prob_p2.eq.0, 2, prob_p2)
 ; prob_p3 = where(prob_p3.eq.0, 2, prob_p3)
 ; prob_p4 = where(prob_p4.eq.0, 2, prob_p4)
 ; prob_p5 = where(prob_p5.eq.0, 2, prob_p5)
 ; prob_p6 = where(prob_p6.eq.0, 2, prob_p6)
 ; prob_p7 = where(prob_p7.eq.0, 2, prob_p7)
 ; prob_p8 = where(prob_p8.eq.0, 2, prob_p8)
;;;;


print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "所有通量的线性趋势-ERAI-newthod")

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
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

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
  res@cnMinLevelValF       =   -8
  res@cnMaxLevelValF       =   8
  res@cnLevelSpacingF      =  1
 

  plot = new(8,graphic)

  res@gsnLeftString = "Net Heat Flux"
  plot(0) = gsn_csm_contour_map(wks, r_total({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "Ocean advction"
  plot(1) = gsn_csm_contour_map(wks, r_total_ocean({lat1:lat2},{lon1:lon2}), res)
  
  res@gsnLeftString = "Latent heat flux"
  plot(2) = gsn_csm_contour_map(wks, r_lh({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "Sensible heat flux"
  plot(3) = gsn_csm_contour_map(wks, r_sh({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "Net Longwave radiation"
  plot(4) = gsn_csm_contour_map(wks, r_nl({lat1:lat2},{lon1:lon2}), res)
  res@gsnLeftString = "Net Shortwave radiation"
  plot(5) = gsn_csm_contour_map(wks, r_ns({lat1:lat2},{lon1:lon2}), res)

  res@gsnLeftString = "Ocean horizonal advction"
  plot(6) = gsn_csm_contour_map(wks, r_hadv({lat1:lat2},{lon1:lon2}), res)
  
  res@gsnLeftString = "Ocean Vertical advction"
  plot(7) = gsn_csm_contour_map(wks, r_vadv, res)

  


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,8,1)

   do i = 0,7
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
   str2 = "0" + ispan(1,8,1)
   do i = 0,7
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
    res_pro@cnLevelSpacingF     = 0.1      ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(8,graphic)

    iplot(0) = gsn_csm_contour(wks,prob_total({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_total_ocean({lat1:lat2},{lon1:lon2}) ,res_pro)
    iplot(2) = gsn_csm_contour(wks,prob_lh({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(3) = gsn_csm_contour(wks,prob_sh({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(4) = gsn_csm_contour(wks,prob_nl({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(5) = gsn_csm_contour(wks,prob_ns({lat1:lat2},{lon1:lon2}),res_pro)

    iplot(6) = gsn_csm_contour(wks,prob_hadv({lat1:lat2},{lon1:lon2}),res_pro)
    iplot(7) = gsn_csm_contour(wks,prob_vadv,res_pro)

    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


    do i = 0,7
     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    end do

    ;;;  
    do i = 0,7
      overlay(plot(i), iplot(i))
    end do 

  pres = True
  pres@gsnPanelLabelBar  = True
  pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f","g","h"/)

  gsn_panel(wks, plot, (/4,2/), pres)
 


end