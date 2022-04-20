

;;**************************************************************
;;;main program
;;***************************************************************
begin


;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc"
    
    f1=addfile(diri, "r")
    salt_r = f1->wt(:,0,:,:)

    time = f1->time
    

    date = cd_calendar(time, 1)
    print(date)
    ; exit()

    salt = tofloat(salt_r)
    copy_VarMeta(salt_r, salt)
    printVarSummary(salt)
   


    ; salt  = where(salt.eq.salt@_FillValue, 10^9, salt)
    ; salt@_FillValue = 10^9
    ; salt@missing_value = 10^9
    ; salt = where(salt.eq.0, salt@_FillValue,salt)
  
    printVarSummary(salt)
   
    salt!1 = "lat"
    salt!2 = "lon"
    salt&lat = f1->lat
    salt&lon = f1->lon
  
    
    salt&lat@units = "degrees_north"
    salt&lon@units = "degrees_east"


    density_clm = dim_avg_n_Wrap(salt, 0)
    
    ; print(salt(1,:,:))
    printVarSummary(density_clm)
   
    ; printMinMax(salt(400,:,:),1)
 
    


  wks = gsn_open_wks("pdf", "垂直速度的气候态")
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  ; res@cnFillPalette="MPL_RdBu"
  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   30.0
  ; res@cnMaxLevelValF       =   37.0
  ; res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(salt&lat)
  res@mpMaxLatF             = max( salt&lat)
  res@mpMinLonF             = min(salt&lon)
  res@mpMaxLonF             = max(salt&lon)
  
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