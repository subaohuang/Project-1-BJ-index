

;;**************************************************************
;;;main program
;;***************************************************************
begin


;;;读入HadiSST海温数据
    diri="/mnt/f/share/SODA/SODA3.4.2_potensial_density_mld.nc"
    
    f1=addfile(diri, "r")
    potensial_density = f1->phro_mld
    lat =f1->lat
    lon =f1->lon

    d_r =potensial_density(time|:,lat|:,lon|:)
    d =tofloat(d_r)
    copy_VarMeta(d_r, d)
    printVarSummary(d)


    
    d = where(d.eq."-nan", 10^9, d)
    d@_FillValue = 10^9
    d&lat@units = "degrees_north"
    d&lon@units = "degrees_east"


    

    density_clm = dim_avg_n_Wrap(d, 0)
    
    density_clm = (/density_clm - 10^3/)
    printVarSummary(density_clm)
    

    


  wks = gsn_open_wks("pdf", "/mnt/f/share/SODA/密度的气候态")
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="matlab_jet"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   20.0
  res@cnMaxLevelValF       =   30.0
  res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(lat)
  res@mpMaxLatF             = max(lat)
  res@mpMinLonF             = min(lon)
  res@mpMaxLonF             = max(lon)
  
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

  res@mpFillOn = False

  plot =gsn_csm_contour(wks,density_clm,res)
  draw(plot)
  frame(wks)


end 

;;;