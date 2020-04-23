

;;**************************************************************
;;;main program
;;***************************************************************
begin


;;;读入HadiSST海温数据
    diri="/mnt/f/share/SODA/SODA3.4.2_temp_5m.nc"
    
    f1=addfile(diri, "r")
    salt = f1->temp
    lat =f1->lat
    lon =f1->lon

    
    
    
    salt&lat@units = "degrees_north"
    salt&lon@units = "degrees_east"




    density_clm = dim_avg_n_Wrap(salt, 0)
    
    
    printVarSummary(density_clm)
    printMinMax(density_clm,1)

    


  wks = gsn_open_wks("pdf", "/mnt/f/share/SODA/温度的气候态")
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  ; res@cnFillPalette="MPL_RdBu"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   30.0
  res@cnMaxLevelValF       =   34.0
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