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
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
 







;;;read NCEP-NCAR heat flux data
    
 

    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/annual/ERAI-heatflux-1degree.nc"
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



    lhtfl_annual = month_to_annual(lhtfl, 1)
    shtfl_annual = month_to_annual(shtfl, 1)
    nswrs_annual = month_to_annual(nswrs, 1)
    nlwrs_annual = month_to_annual(nlwrs, 1)
   
   lhtfl_clm = dim_avg_n_Wrap(lhtfl_annual, 0)
     total_annual = lhtfl_annual
    total_annual  = lhtfl_annual + shtfl_annual + nswrs_annual + nlwrs_annual









print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "ERA-Heatflux的气候态")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  ; res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
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

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -8
  ; res@cnMaxLevelValF       =   8
  ; res@cnLevelSpacingF      =  1
 

  plot = new(1,graphic)
  res@gsnLeftString = "Latent heating"
  plot(0) = gsn_csm_contour_map(wks, lhtfl_clm, res)

  draw(plot)
  frame(wks)
 


end