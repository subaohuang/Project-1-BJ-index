begin
  diri="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"giss_svd.nc")
  f1=addfile(fils1, "r")
  land_svd=f1->land_svd
  sst_svd=f1->sst_svd
  land_svd_timef=f1->land_svd_timef
  sst_svd_timef=f1->sst_svd_timef
  svd_var=f1->svd_var
 
  lon=f1->lon
  lat=f1->lat
  
  printVarSummary(land_svd)
  printVarSummary(sst_svd)
  
  land_svd&lat@units="degrees_north"
  land_svd&lon@units="degrees_east"
  land_svd@_FillValue=999.0
  sst_svd@_FillValue=999.0
  copy_VarCoords(land_svd, sst_svd)
 
  land_svd_new=land_svd(lat|:,lon|:,lev|:)
  sst_svd_new=sst_svd(lat|:,lon|:,lev|:)

  printVarSummary(sst_svd_new)

  corr0=escorc(land_svd_timef(:,0),sst_svd_timef(:,0))
  corr1=escorc(land_svd_timef(:,1),sst_svd_timef(:,1))
  
  mode=1

;;;画图
  ; wks=gsn_open_wks("pdf", ""+diri+"giss_sst_svd1&2"+(/mode+1/))
  wks=gsn_open_wks("pdf", ""+diri+"giss_sst&land_svd2")
  res = True
  res@gsnAddCyclic=False
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnLeftString="SVD"+(/mode+1/)
  res@gsnRightString=""+svd_var(mode)
  ; res@cnFillPalette  = "temp_diff_18lev"
  res@cnFillPalette="BlueRed"
  ; res@cnFillPalette="precip_diff_12lev"
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False
 
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.02
  res@cnMaxLevelValF       =   0.02
  res@cnLevelSpacingF      =  0.0025

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  ; res@cnMinLevelValF       =   -0.40
  ; res@cnMaxLevelValF       =   0.40
  ; res@cnLevelSpacingF      =  0.08

  res@lbLabelBarOn=True
  res@lbOrientation   ="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(land_svd&lat)
  res@mpMaxLatF             = max(land_svd&lat)
  res@mpMinLonF             = min(land_svd&lon)
  res@mpMaxLonF             = max(land_svd&lon)
  res@mpCenterLonF          =avg(land_svd&lon)
  
  res@tmXTOn=False
  res@tmYROn=False
  res@tmYLMode="Explicit"
  res@tmYLValues=ispan(-60,60,30)
  res@tmYLMinorValues=ispan(-60,60,10)
  res@tmYLMinorOn=True
  res@tmYLLabels=(/"60~S~o~N~S","30~S~o~N~S","0","30~S~o~N~N","60~S~o~N~N"/)
  res@tmYLLabelFontHeightF=0.01

  res@tmXBMode="Explicit"
  res@tmXBValues=fspan(0.5,358,13)
  res@tmXBMinorOn=False  
  res@tmXBLabelFontHeightF=0.01
  res@tmXBLabelFontThicknessF=0.02
  res@tmXBLabels=(/"0","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
            "180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0"/)
  
  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.2
  ; res@vpWidthF=0.4


  plot=new(3,graphic)


  plot(0)=gsn_csm_contour_map(wks, sst_svd_new(:,:,mode), res)
  ; res@cnFillPalette="precip_diff_12lev"
 
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.04
  res@cnMaxLevelValF       =   0.04
  res@cnLevelSpacingF      =  0.005

  plot(1)=gsn_csm_contour_map(wks, land_svd_new(:,:,mode), res)
  
  lineres=True 
  lineres@gsnDraw=False
  lineres@gsnFrame=False
  lineres@gsnYRefLine=0.0
  lineres@xyLineColors=(/"red","blue"/)
  lineres@xyDashPatterns=(/0,0/)
 
  lineres@tmXTOn=False
  lineres@tmYROn=False
  lineres@tmXBMode="Explicit"
  lineres@tmXBMinorValues=ispan(1,55,1)
  lineres@tmXBValues=ispan(1,55,10)
  lineres@tmXBLabels=ispan(1958,2012,10)
  lineres@tmYLLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  ; lineres@vpHeightF=0.18
  ; lineres@vpWidthF=0.46
  lineres@vpHeightF=0.25
  lineres@vpWidthF=0.75
  
  lineres@pmLegendDisplayMode="Always"
  lineres@xyExplicitLegendLabels=(/"Land","SST"/)
  lineres@lgJustification        = "TopRight"
  lineres@pmLegendOrthogonalPosF = -0.7     ; Move into plot
  lineres@pmLegendParallelPosF   = 1.1      ; Move to right
  lineres@pmLegendWidthF         = 0.3       ; Decrease width
  lineres@pmLegendHeightF        = 0.1       ; Decrease height
  lineres@lgBoxMinorExtentF      = 0.1       ; Shorten the legend lines
  lineres@lgLabelFontHeightF     = 0.015     ; Change the font size
  lineres@lgPerimOn=False
  lineres@gsnRightString="corr="+corr1

  plot(2)=gsn_csm_xy(wks, ispan(1,55,1), (/land_svd_timef(:,mode),sst_svd_timef(:,mode)/), \
                     lineres)
  


  ; plot(3)=gsn_csm_xy(wks, ispan(1,55,1), sst_svd_timef(:,1), \
  ;                    lineres)


  pres=True

  pres@gsnPanelLabelBar=False
  pres@gsnPanelDebug=True
  gsn_panel(wks, plot, (/3,1/), pres)



end


