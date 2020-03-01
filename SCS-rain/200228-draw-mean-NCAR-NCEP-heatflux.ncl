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
  

    timeind:=ind(date.le.201312.and.date.ge.198001)
    

    lhtfl = f2->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f2->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f2->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f2->nlwrs(timeind,{lat1:lat2},{lon1:lon2})
    

    ; fils1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/hflx/lhtfl.sfc.mon.mean.nc"
    ; fils2 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/hflx/shtfl.sfc.mon.mean.nc"
    ; fils3 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/hflx/nswrs.sfc.mon.mean.nc"
    ; fils4 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/hflx/nlwrs.sfc.mon.mean.nc"


    ; f1 = addfile(fils1, "r")
    ; f2 = addfile(fils2, "r")
    ; f3 = addfile(fils3, "r")
    ; f4 = addfile(fils4, "r")
    f6=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_mlt_mn_1980-2015_chazhi.nc", "r") 


    indyear = (2013-1980+1)*12-1
    mlt = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})
   


    ;;;nswrs_pen 
     R = 0.58
     L1 = 0.35
     L2 = 23
     nswrs_pen = nswrs
     nswrs_pen =  nswrs*(R*exp(-mlt/L1)+ (1-R)*exp(-mlt/L2))
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

    ;;;转化成方向向下

    lhtfl = (/-1*lhtfl/)
    shtfl = (/-1*shtfl/)
    nswrs = (/-1*nswrs/)
    nlwrs = (/-1*nlwrs/)
    nswrs_pen = (/-1*nswrs_pen/)



    lhtfl_annual = month_to_annual(lhtfl, 1)
    shtfl_annual = month_to_annual(shtfl, 1)
    nswrs_annual = month_to_annual(nswrs, 1)
    nlwrs_annual = month_to_annual(nlwrs, 1)
    nswrs_pen_annual = month_to_annual(nswrs_pen, 1)
    total_annual = lhtfl_annual
    total_annual  = lhtfl_annual + shtfl_annual + nswrs_annual + nlwrs_annual -nswrs_pen_annual
    printMinMax(total_annual, 1)
    exit()



;;;;
 printVarSummary(dsst_index)

print("*******开始计算平均*****************")
;;;calculate correlation coef

    ;;;correlation with sst
    r_sst_lh =  dim_avg_n_Wrap(lhtfl_annual, 0)
    r_sst_sh   = dim_avg_n_Wrap(shtfl_annual, 0)

    r_sst_ns =  dim_avg_n_Wrap(nswrs_annual, 0)
    r_sst_nl   = dim_avg_n_Wrap(nlwrs_annual, 0)
    r_sst_total = dim_avg_n_Wrap(total_annual, 0)
    ;;;;


    copy_VarCoords(ssta_annual(0,:,:), r_sst_lh)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_sh)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_ns)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_nl)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_total)
    



    ;;;;;
;;;;

print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "NCEP-NACR-Heatflux的long-term-mean-摄氏度")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
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
  res@cnMinLevelValF       =   -50
  res@cnMaxLevelValF       =   50
  res@cnLevelSpacingF      =  10
 

  plot = new(5,graphic)
  res@gsnLeftString = "Latent heating"
  plot(0) = gsn_csm_contour_map(wks, r_sst_lh, res)
  res@gsnLeftString = "Sensible heating"
  plot(1) = gsn_csm_contour_map(wks, r_sst_sh, res)
  res@gsnLeftString = "Net Longwave radiation"
  plot(2) = gsn_csm_contour_map(wks, r_sst_nl, res)
  res@gsnLeftString = "Net Shortwave radiation"
  plot(3) = gsn_csm_contour_map(wks, r_sst_ns, res)

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -8
  res@cnMaxLevelValF       =   8
  res@cnLevelSpacingF      =  2
 
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

   
  pres = True
  pres@gsnPanelLabelBar  = False

  gsn_panel(wks, plot, (/3,2/), pres)
 


end