begin
;;;读入HadiSST海温数据

  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201212.and.date.ge.195801)
  sst=lonFlip(f1->sst(timeind,{20:-20},:))
  sst=(/rmMonAnnCycTLL(sst)/)
  delete(sst&time)
  sst&time=date(timeind)
;;;
 sst_djf=month_to_season(sst, "DJF")
 sst_djf&time=ispan(1957,2011,1)
  printVarSummary(sst_djf)
printVarSummary(sst_djf&time)

  ; canonical_yr=(/1965,1972,1976,1982,1997/)
 canonical_yr=1998
 caind=get1Dindex(sst_djf&time,canonical_yr)
 ; ss_djf=dim_avg_n((/sst_djf(caind,:,:)/),0)
 ss_djf=sst_djf(caind,:,:)
 copy_VarCoords(sst_djf(1,:,:), ss_djf)


 sst_1998=dim_avg_n_Wrap(sst({199812:199902},:,:),0)

;;;画图
  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/hadisst")
  ; gsn_define_colormap(wks, "MPL_YlOrBr")
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

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  ; res@cnMinLevelValF       =   0.08
  ; res@cnMaxLevelValF       =   0.4
  ; res@cnLevelSpacingF      =  0.04

  res@tmXTOn=False
  res@tmYROn=False
  
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(sst_djf&latitude)
  res@mpMaxLatF             = max(sst_djf&latitude)
  res@mpMinLonF             = min(sst_djf&longitude)
  res@mpMaxLonF             = max(sst_djf&longitude)
  res@mpCenterLonF          =avg(sst_djf&longitude)
  res@mpFillDrawOrder="PreDraw"
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.8
  res@vpHeightF = 0.3


  plot1=gsn_csm_contour_map(wks, ss_djf, res)


;;;添加区域
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgx=(/80,170,170,80,80/)
   pgy=(/-5,-5,5,5,-5/)

   dum = gsn_add_polyline(wks, plot1, pgx, pgy, pgres)

   draw(plot1)
   frame(wks)
  end
;;;