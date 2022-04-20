    ;;;Icoads 云量
    begin
    model1="FGOALS-s2"
    model2="inmcm4"
    diri_input="/mnt/f/share/CMIP5/AMIP/ts_new/"

    f1=addfile(diri_input+"ts_"+model1+"_amip_r1i1p1_197901-200812.nc", "r")
    f2=addfile(diri_input+"ts_"+model2+"_amip_r1i1p1_197901-200812.nc", "r")
    
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    ts1=f1->ts(timeind,:,:)
    ts2=f2->ts(timeind,:,:)
    
    ts1 = (/ts1-273.15/)
    ts2 = (/ts2-273.15/)

    delete(ts1&time)
    delete(ts2&time)
    ts1&time=date(timeind)
    ts2&time=date(timeind)
    printVarSummary(ts1)
    
    ts_pi1=lonPivot(ts1, 30)
    ts_pi2=lonPivot(ts2, 30)
    ts_djf_r1=month_to_season(ts_pi1(time|:,lat|:,lon|:), "DJF")
    ts_djf_r2=month_to_season(ts_pi2(time|:,lat|:,lon|:), "DJF")


    ;;注意点
    ts_djf1=ts_djf_r1(0:,:,:)
    delete(ts_djf1&time)
    ts_djf1&time=ispan(1979,2008,1)

    ts_djf2=ts_djf_r2(0:,:,:)
    delete(ts_djf2&time)
    ts_djf2&time=ispan(1979,2008,1)
    
    ts_djf1_clm=dim_avg_n_Wrap(ts_djf1, 0)
    ts_djf2_clm=dim_avg_n_Wrap(ts_djf2, 0)

    ts_djff=ts_djf1_clm
    ts_djff=ts_djf2_clm - ts_djf1_clm


  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/amip冬季气温的气候态比较")


  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftString="(a) 1979-1994 mean" 
  res@cnFillPalette="NCV_blu_red"
  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  .1
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  ; res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             =30
  res@mpMaxLonF             = 390
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  res@mpShapeMode="FreeAspect"
  res@vpHeightF=0.3
  res@vpWidthF=0.7  
  plot=new(3,graphic)
  ; plot(0)=gsn_csm_contour_map(wks, ts_djf1_clm, res)
  ; res@gsnLeftString="(a) 1995-2012 mean" 
  ; plot(1)=gsn_csm_contour_map(wks, ts_djf2_clm, res)
  ;   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  .1
  plot(2)=gsn_csm_contour_map(wks, ts_djff, res)
draw(plot(2))
frame(wks)

  ; pres=True
  ; pres@gsnPanelLabelBar=True
  ; gsn_panel(wks, plot, (/3,1/), pres)

end 



  
