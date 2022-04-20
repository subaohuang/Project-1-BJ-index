begin
;;;读入HadiSST海温数据

  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201012.and.date.ge.195001)
  sst=lonFlip(f1->sst(timeind,{20:-20},:))
  delete(sst&time)
  sst&time=ispan(1,732,1)
  
  printVarSummary(sst)

  sst_annual=month_to_annual(sst, 1)
  printVarSummary(sst_annual)

  sst_sm=runave_n_Wrap(sst_annual, 9, 0, 0)

  sst_std=dim_stddev_n_Wrap(sst_annual, 0)

  sst_sm_std=dim_stddev_n_Wrap(sst_sm, 0)
;;;



;;;画图
  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/sst的方差变化")
  gsn_define_colormap(wks, "MPL_YlOrBr")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnLeftString="sst std"
  res@gsnRightString="1950-2010"
  
  res@cnLineColors="Transparent"  ;
  res@lbLabelBarOn=True
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  res@cnFillMode="RasterFill"
  res@cnRasterSmoothingOn=True
  res@cnFillDrawOrder="PreDraw"

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   0.08
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.04

  res@tmXTOn=False
  res@tmYROn=False
  
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(sst_std&latitude)
  res@mpMaxLatF             = max(sst_std&latitude)
  res@mpMinLonF             = min(sst_std&longitude)
  res@mpMaxLonF             = max(sst_std&longitude)
  res@mpCenterLonF          =avg(sst_std&longitude)
  res@mpFillDrawOrder="PreDraw"
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.8
  res@vpHeightF = 0.3


  plot1=gsn_csm_contour_map(wks, sst_std, res)
  ; draw(plot1)
  plot2=gsn_csm_contour_map(wks, sst_sm_std, res)

;;;添加区域
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgx=(/80,170,170,80,80/)
   pgy=(/-5,-5,5,5,-5/)

   dum = gsn_add_polyline(wks, plot2, pgx, pgy, pgres)

   draw(plot2)
   frame(wks)
  end
;;;