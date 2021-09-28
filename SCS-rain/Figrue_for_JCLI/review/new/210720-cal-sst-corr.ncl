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
;;**************************************************************
;;;main program
;;***************************************************************
begin
;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.195001)

    yearnum = 2015-1950+1 

    sst  = lonFlip(f1->sst(timeind,:,:))

    ; sst = (/rmMonAnnCycTLL(sst)/)
    
    sst_clm = clmMonTLL(sst)  ;;;12,:,:
  ;;;

  ;;;cal change 
    ssta = sst
    ssta = (/rmMonAnnCycTLL(ssta)/)

    dim = dimsizes(sst)
   
    ssta_12 = new((/yearnum,12,dim(1),dim(2)/),float)

    do i=0,11   
    ssta_12(:,i,:,:) = ssta(i::12,:,:)
    end do 
    
    rc_12 = regCoef_n(ispan(1,yearnum,1), ssta_12, 0, 0)
    
    corr = sst(1,:,:)
    corr = escorc_n(sst_clm, rc_12, 0, 0)

    prob = rtest(corr, 12, 0)

    corr=  where(ismissing(corr), 0, corr)
    corr = where(prob.le.0.05, corr, 0)
 
  wks = gsn_open_wks("pdf", "corr")
  
   ;;;
  res = True
  res@gsnAddCyclic = True
  res@gsnFrame=False
  res@gsnDraw=False
  res@cnFillPalette="BlueWhiteOrangeRed"
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=1
  res@pmLabelBarOrthogonalPosF = 0.17
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0.5
  res@mpMaxLonF             = 359.5
  res@lbOrientation="horizontal"
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@mpFillOn = False
 
  res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res@cnMinLevelValF      = -1.0        ; set min contour level
  res@cnMaxLevelValF      = 1.0        ; set max contour level
  res@cnLevelSpacingF     = 0.2       ; set contour spacing


  res@lbLabelFontHeightF = 0.02
  res@lbAutoManage = False


  plot = gsn_csm_contour_map(wks, corr, res)

 draw(plot)
 frame(wks)


end 