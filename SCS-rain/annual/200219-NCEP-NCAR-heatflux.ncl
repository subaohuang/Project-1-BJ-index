;;********************************************************
;;NCEP-NCAR heat-fliu的气候态，global region
;;
;;                                         lnshuheng 2020.02.19
;**********************************************************

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
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
 

;;;read NCEP-NCAR heat flux data
    
    fils2 = "/home/ys17-19/lsh/Project/SCS-rain/NCEP-NCAR-heatflux-1degree.nc"
    f2 = addfile(fils2, "r")

    date:=cd_calendar(f2->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)
  
    
    lhtfl = f2->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f2->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f2->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f2->nlwrs(timeind,{lat1:lat2},{lon1:lon2})
    
    lhtfl&lat@units = "degrees_north"
    lhtfl&lon@units = "degrees_east"

    copy_VarCoords(lhtfl, shtfl)
    copy_VarCoords(lhtfl, nswrs)
    copy_VarCoords(lhtfl, nlwrs)

    printVarSummary(lhtfl)


    lhtfl_clm = dim_avg_n_Wrap(lhtfl, 0)
    shtfl_clm = dim_avg_n_Wrap(shtfl, 0)
    nswrs_clm = dim_avg_n_Wrap(nswrs, 0)
    nlwrs_clm = dim_avg_n_Wrap(nlwrs, 0)
;;;;



print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "NCEP-NACR-Heatflux气候态")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlueWhiteOrangeRed"
  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -90
  res@mpMaxLatF             = 90
  res@mpMinLonF             =0
  res@mpMaxLonF             = 360


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
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.1
 

  plot = new(4,graphic)
  res@gsnLeftString = "Latent heating"
  plot(0) = gsn_csm_contour_map(wks, lhtfl_clm, res)
  res@gsnLeftString = "Sensible heating"
  plot(1) = gsn_csm_contour_map(wks, shtfl_clm, res)
  res@gsnLeftString = "Net Longwave radiation"
  plot(2) = gsn_csm_contour_map(wks, nlwrs_clm, res)
  res@gsnLeftString = "Net Shortwave radiation"
  plot(3) = gsn_csm_contour_map(wks, nswrs_clm, res)

  
  


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,8,1)

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

  
  pres = True
  pres@lbLabelBarOn  = False

  gsn_panel(wks, plot, (/2,2/), pres)
 


end