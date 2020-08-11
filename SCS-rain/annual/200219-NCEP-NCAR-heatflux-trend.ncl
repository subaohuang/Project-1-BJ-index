;;********************************************************
;;计算SST倾向以及SST值与净长波辐射，短波辐射，感热，潜热以及平流和对流,
;;其中热通量的方向都是向上为正，向下为负。。
;;;的回归场，热通量数据是NCEP-NCAR数据，海洋数据是SODA2.2.4,SST是哈哈得来SST数据
;;由于数据之间的精度不同，需要进行插值预处理,SODA精度是0.5×0.5,HADISST是1×1
;;NCEP-NCAR的精度为1.87×1.9
;;s
;;                                         lnshuheng 2020.02.18
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
    startmon =197901
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

    timeind:=ind(date.le.201312.and.date.ge.197901)

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






;;;read NCEP-NCAR heat flux data
    
 
    fils2 = "/home/ys17-19/lsh/Project/SCS-rain/NCEP-NCAR-heatflux-1degree.nc"
    f2 = addfile(fils2, "r")

    date:=cd_calendar(f2->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)
  
    
    lhtfl = f2->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f2->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f2->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f2->nlwrs(timeind,{lat1:lat2},{lon1:lon2})

    lhtfl  = (/rmMonAnnCycTLL(lhtfl)/)
    shtfl  = (/rmMonAnnCycTLL(shtfl)/)
    nswrs  = (/rmMonAnnCycTLL(nswrs)/)
    nlwrs  = (/rmMonAnnCycTLL(nlwrs)/)

    ;;;转化成方向向下

    lhtfl = (/-1*lhtfl/)
    shtfl = (/-1*shtfl/)
    nswrs = (/-1*nswrs/)
    nlwrs = (/-1*nlwrs/)



    lhtfl_annual = month_to_annual(lhtfl, 1)
    shtfl_annual = month_to_annual(shtfl, 1)
    nswrs_annual = month_to_annual(nswrs, 1)
    nlwrs_annual = month_to_annual(nlwrs, 1)
    total_annual = lhtfl_annual
    total_annual  = lhtfl_annual + shtfl_annual + nswrs_annual + nlwrs_annual



;;;;
 printVarSummary(dsst_index)

print("*******开始计算回归系数*****************")
;;;calculate correlation coef

    ;;;correlation with sst
    r_sst_lh =  regCoef_n(ispan(1, yearnum, 1), lhtfl_annual, 0, 0)
    r_sst_sh   = regCoef_n(ispan(1, yearnum, 1), shtfl_annual, 0, 0)

    r_sst_ns =  regCoef_n(ispan(1, yearnum, 1), nswrs_annual, 0, 0)
    r_sst_nl   = regCoef_n(ispan(1, yearnum, 1), nlwrs_annual, 0, 0)
    r_sst_total = regCoef_n(ispan(1, yearnum, 1), total_annual, 0, 0)
    ;;;;


    copy_VarCoords(ssta_annual(0,:,:), r_sst_lh)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_sh)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_ns)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_nl)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_total)
    
    r_sst_lh = (/r_sst_lh*yearnum/)
    r_sst_sh = (/r_sst_sh*yearnum/)
    r_sst_ns = (/r_sst_ns*yearnum/)
    r_sst_nl = (/r_sst_nl*yearnum/)
    r_sst_total = (/r_sst_total*yearnum/)
 


    ;;;;;
;;;;

print("*******显著性检验*****************")

;;;t test 
   prob_p1= rc_test(r_sst_lh)
   prob_p2= rc_test(r_sst_sh)
   prob_p3= rc_test(r_sst_nl)
   prob_p4= rc_test(r_sst_ns)
   prob_p5= rc_test(r_sst_total)

 
  
  copy_VarCoords(r_sst_lh, prob_p1)
  copy_VarCoords(r_sst_lh, prob_p2)
  copy_VarCoords(r_sst_lh, prob_p3)
  copy_VarCoords(r_sst_lh, prob_p4)
  copy_VarCoords(r_sst_total, prob_p5)
  
  
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
  wks   = gsn_open_wks("pdf", "NCEP-NACR-Heatflux的线性趋势")

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
  res@cnMinLevelValF       =   -8
  res@cnMaxLevelValF       =   8
  res@cnLevelSpacingF      =  1
 

  plot = new(5,graphic)
  res@gsnLeftString = "Latent heating"
  plot(0) = gsn_csm_contour_map(wks, r_sst_lh, res)
  res@gsnLeftString = "Sensible heating"
  plot(1) = gsn_csm_contour_map(wks, r_sst_sh, res)
  res@gsnLeftString = "Net Longwave radiation"
  plot(2) = gsn_csm_contour_map(wks, r_sst_nl, res)
  res@gsnLeftString = "Net Shortwave radiation"
  plot(3) = gsn_csm_contour_map(wks, r_sst_ns, res)
  res@gsnLeftString = "Total Heat Flux"
  plot(4) = gsn_csm_contour_map(wks, r_sst_total, res)
  


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
   str2 = "0" + ispan(1,8,1)
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
    res_pro@cnLevelSpacingF     = 0.1      ; set contour spacing

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