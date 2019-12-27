;;********************************************************
;;回归ERA的u,v,slp,z到降水的变化上
;;
;;lnshuheng 2019.12.24
;**********************************************************
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
    
    p = (/rmMonAnnCycTLL(p)/)
    p=(/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
    

    p_annual = month_to_annual(p, 1)
    

    ;;;计算区域的降水平均，得到一个index,先看降水由哪里因素决定
     pre_index =dim_avg_n_Wrap(p_annual(:,{0:15},{110:155}), (/1,2/))
     pre_index = (/dim_standardize(pre_index, 1)/)
    ;;;
    delete(f2)
    delete(timeind)
    delete(date)

     ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)

    sst=lonFlip(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_annual = month_to_annual(ssta, 1)

    ;;;;;回归pre
    
   ;;;回归
    sst_rc = regCoef_n(pre_index, ssta_annual, 0, 0)
    copy_VarCoords(ssta_annual(1,:,:), sst_rc)
    prob_sst = rc_test(sst_rc)
    sst_rc = where(prob_sst.le.0.05, sst_rc, sst_rc@_FillValue)

    

  wks = gsn_open_wks("pdf", "SST回归到pre")

  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnLeftString="(a) OLR"
  res@cnFillPalette="BlueWhiteOrangeRed"
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



  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.5
  res@cnMaxLevelValF       =   0.5
  res@cnLevelSpacingF      =  0.05
  res@gsnLeftString="(b)GPCP Precipitation"
  plot = gsn_csm_contour_map(wks, sst_rc, res)

   draw(plot)
   frame(wks)
end
