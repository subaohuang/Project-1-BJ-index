 begin
   diri="/mnt/f/share/NOAA/"
   fils1=systemfunc("ls "+diri+"olr.mon.mean.nc")
    ; fils1=systemfunc("ls "+diri+"/MERRA/MERRA/msf_time_MERRA_global.nc")
    ; fils1=systemfunc("ls "+diri+"MRE/msf_time_mre_global.nc")
    year1=1979
    year2=2012
    yearend=201212
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    olr_r=f1->olr(timeind,:,:)
    olr=(/olr_r@add_offset+olr_r@scale_factor*olr_r/)
    copy_VarCoords(olr_r, olr)
    copy_VarAtts(olr_r, olr)
    delete(olr&time)
    olr&time=date(timeind)
    printVarSummary(olr)
    
   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    olr_jja_r=dim_avg_n_Wrap(month_to_season(olr(time|:,lat|:,lon|:), "JJA"),0)
    olr_djf_r=dim_avg_n_Wrap(month_to_season(olr(time|:,lat|:,lon|:), "MAM"),0)

    ; olr_djf_r=(/month_to_annual(olr_pi(time|:,lat|:,lon|:), 1)/)
    ;;;

 wks   = gsn_open_wks("pdf", "/mnt/f/share/era/data/OLR夏季与冬季的气候态")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(2,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnLeftString="(a) OLR"
  res@cnFillPalette="NCV_blu_red"
  ; res@cnFillOn=True
  ; res@cnLinesOn=False
  res@cnFillOn=False
  res@cnLinesOn=True

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -8.0
  ; res@cnMaxLevelValF       =   8.0
  ; res@cnLevelSpacingF      =  1.0

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             =0
  res@mpMaxLonF             = 360
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

  res@mpShapeMode="FreeAspect"
  res@vpHeightF=0.3
  res@vpWidthF=0.7
  plot(0)=gsn_csm_contour_map(wks, olr_jja_r, res)
  plot(1)=gsn_csm_contour_map(wks, olr_djf_r, res)
 

  pres=True
  gsn_panel(wks, plot, (/2,1/), pres)
end