function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

begin
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    p=f2->precip(timeind,:,:)
    p=(/rmMonAnnCycTLL(p)/)
    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
    p_annual = month_to_annual(p(time|:,lat|:,lon|:), 1)
    ;;注意点
     
 
    ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)

    sst=f1->sst(timeind,{20:60},{-70:-20})
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_annual = month_to_annual(ssta, 1)

    sst_index = dim_avg_n_Wrap(ssta_annual, (/1,2/))
    printVarSummary(sst_index)
    ;;;

    ;;;;
    f4 = addfile("/home/yangsong3/data-observation/ERA-interim/monthly/slp.mon.mean.nc","r")
    date :=cd_calendar(f4->time, 1)
    timeind :=ind(date.ge.197901.and.date.le.201312)
    z = short2flt(f4->msl(timeind,:,:))   
    print("=====finish read data==========")
    ;;;

    z = (/rmMonAnnCycTLL(z)/)
    z=(/dtrend_msg_n(ispan(1,dimsizes(z&time),1),z,False,True,0)/)
    
    z_annual = month_to_annual(z, 1)


     
    ;;;回归pre
     p_rc = regCoef_n(sst_index, p_annual, 0, 0)
     copy_VarCoords(p_annual(1,:,:), p_rc)

      z_rc = regCoef_n(sst_index, z_annual, 0, 0)
     copy_VarCoords(z_annual(1,:,:), z_rc)
    ;;
    
    ;;检验
    prob_p = rc_test(p_rc)
    p_rc = where(prob_p.le.0.05, p_rc, p_rc@_FillValue)

    ;;;
wks = gsn_open_wks("pdf", "降水回归大西洋的sst")

  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnLeftString="(a) OLR"
  res@cnFillPalette="NCV_blu_red"
  res@cnFillOn=True
  res@cnLinesOn=False

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

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.7

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -2.5
  ; res@cnMaxLevelValF       =   2.5
  ; res@cnLevelSpacingF      =  0.25

  res@cnFillPalette="MPL_BrBG"
  res@gsnLeftString="(b)GPCP Precipitation"
  plot = gsn_csm_contour_map(wks, z_rc, res)

   draw(plot)
   frame(wks)
end



    ;;;




    

