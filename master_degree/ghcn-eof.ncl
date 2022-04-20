begin
;;;读入hadisst
  diri="/mnt/f/share/NOAA/GHCN_CAMS/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"air.mon.mean.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201212.and.date.ge.195801)
  sst=f1->air(timeind,{60:-60},:)
  sst=(/rmMonAnnCycTLL(sst)/)
  sst_annual=month_to_annual(sst, 1)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  printVarSummary(sst_annual)
  sst_sm=runave_n_Wrap(sst_annual, 3, 1, 0)
  printVarSummary(sst_sm)
  
   sst_sm=(/dim_standardize_n_Wrap(sst_sm, 0, 0)/)


  ;;插值到2.5乘2.5网格
  delete(diri)  
  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)
  delete(sst)
  

  optEOF      = True
  optEOF@jopt = 0

  neval=30
  ev_sst=eofunc_n_Wrap(sst_sm,neval,True,0)    ;;;得到的标准化后的场，要使得不标准化
  printVarSummary(ev_sst)
  
  
  do ne=0,neval-1
     ev_sst(ne,:,:) = ev_sst(ne,:,:)*sqrt(ev_sst@eval(ne)) ;;;使得量级能够相同
  end do

  ev_timef=eofunc_ts_n_Wrap(sst_sm, ev_sst, False,0)
  printVarSummary(ev_timef)
  ev_timef_std=new((/2,dimsizes(ev_timef&year)/), double)
  eval=ev_sst@eval
  tsmean=ev_timef@ts_mean
  ev_timef_std(0,:)=ev_timef(0,:)/eval(0)
  ev_timef_std(1,:)=ev_timef(1,:)/eval(1)


mode=1
;;;画图
 wks=gsn_open_wks("pdf", ""+diri_pi+"GHCH_CAM_eof"+(/mode+1/))
  res = True
  res@gsnAddCyclic=False
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnLeftString="EOF"+(/mode+1/)
  res@gsnRightString=""+ev_sst@pcvar(mode)+""
  ; res@cnFillPalette  = "temp_diff_18lev"
  res@cnFillPalette="BlueRed"
  ; res@cnFillPalette="precip_diff_12lev"
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False
 
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  ; res@cnMinLevelValF       =   -0.40
  ; res@cnMaxLevelValF       =   0.40
  ; res@cnLevelSpacingF      =  0.08

  res@lbLabelBarOn=True
  res@lbOrientation   ="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ev_sst&lat)
  res@mpMaxLatF             = max(ev_sst&lat)
  res@mpMinLonF             = min(ev_sst&lon)
  res@mpMaxLonF             = max(ev_sst&lon)
  res@mpCenterLonF          =avg(ev_sst&lon)
  
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


  plot=new(2,graphic)


  plot(0)=gsn_csm_contour_map(wks, ev_sst(mode,:,:), res)
  res@cnFillPalette="precip_diff_12lev"
  ; plot(1)=gsn_csm_contour_map(wks, ev_sst(0,:,:,1), res)
  
  lineres=True 
  lineres@gsnDraw=False
  lineres@gsnFrame=False
  lineres@gsnYRefLine=0.0
  lineres@xyLineColors="blue"
  lineres@xyDashPatterns=0
 
  lineres@tmXTOn=False
  lineres@tmYROn=False
  lineres@tmXBMode="Explicit"
  lineres@tmXBMinorValues=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  lineres@tmXBValues=ispan(1,dimsizes(sst_annual(:,1,1)),10)
  lineres@tmXBLabels=ispan(1958,2012,10)
  lineres@tmYLLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  lineres@vpHeightF=0.25
  lineres@vpWidthF=0.75
  

  plot(1)=gsn_csm_xy(wks, ispan(1,dimsizes(sst_annual(:,1,1)),1), ev_timef_std(mode,:), \
                     lineres)
  
  pres=True
  pres@gsnPanelLabelBar=False
  pres@gsnPanelDebug=True
  gsn_panel(wks, plot, (/2,1/), pres)



end