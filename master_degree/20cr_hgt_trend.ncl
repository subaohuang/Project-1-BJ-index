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
;************************************************
; 读入数据
;************************************************
	diri="/mnt/f/share/era/data/"
	fils1=systemfunc("ls "+diri+"msf_time_era_global.nc")
	f1=addfile(fils1, "r")
	fils2=systemfunc("ls /mnt/f/share/NOAA/hgt.mon.mean.nc") 
	f2=addfile(fils2, "r")
	date=cd_calendar(f2->time, 1)
	timeind=ind(date.le.201012.and.date.ge.195001)
	hgt=f2->hgt(timeind,:,{20:-20},:)
	delete(hgt&time)
	hgt&time=date(timeind)

	
   hgt_avg=hgt(:,{200},:,:)

	hgt_50to78=hgt_avg({195001:197812},:,:)
	hgt_79to10=hgt_avg({197901:201012},:,:)
	hgt_50to10=hgt_avg({195001:201012},:,:)

   printVarSummary(hgt_50to78)
 
;;;

;************************************************
; 计算回归
;************************************************
	rc1 =regCoef_n(ispan(1, dimsizes(hgt_50to78&time), 1),hgt_50to78,0,0)
	rc2 =regCoef_n(ispan(1, dimsizes(hgt_79to10&time), 1),hgt_79to10,0,0)
	rc3 =regCoef_n(ispan(1, dimsizes(hgt_50to10&time), 1),hgt_50to10,0,0)
     
   
	copy_VarCoords(hgt_50to78(0,:,:), rc1)
	copy_VarCoords(hgt_50to78(0,:,:), rc2)
	copy_VarCoords(hgt_50to78(0,:,:), rc3)

;************************************************
; 显著性检验
;************************************************

  rc1=(/rc1*120/) 
  rc2=(/rc2*120/) 
  rc3=(/rc3*120/) 
	prob1=rc_test(rc1)
	prob2=rc_test(rc2)
	prob3=rc_test(rc3)

   

;************************************************
;  画图
;************************************************

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/20cr位势高度的趋势")
  gsn_define_colormap(wks, "ViBlGrWhYeOrRe")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnLeftString="w trend"
  res@gsnRightString="1950-1978"
  
  res@cnLineColors="Transparent"  ;
  res@lbLabelBarOn=False
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  ; res@cnFillMode="RasterFill"
  ; res@cnRasterSmoothingOn=True
  
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  
  ; res@cnMinLevelValF       =   -16.0
  ; res@cnMaxLevelValF       =   16.0
  ; res@cnLevelSpacingF      =  2.0

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  
  res@cnMinLevelValF       =   -3.6
  res@cnMaxLevelValF       =   3.6
  res@cnLevelSpacingF      =  0.4


  res@tmXTOn=False
  res@tmYROn=False

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(rc1&lat)
  res@mpMaxLatF             = max(rc1&lat)
  res@mpMinLonF             = min(rc1&lon)
  res@mpMaxLonF             = max(rc1&lon)
  res@mpCenterLonF          =avg(rc1&lon)
  
    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.8
    res@vpHeightF = 0.2

	plot=new(3,graphic)
	plot(0)= gsn_csm_contour_map(wks,rc1,res)
    res@gsnRightString="1979-2010"
	plot(1)= gsn_csm_contour_map(wks,rc2,res)
    res@gsnRightString="1950-2010"
    ; res@cnMinLevelValF       =   -5.0
    ; res@cnMaxLevelValF       =   5.0
    ; res@cnLevelSpacingF      =  1.0

    plot(2)= gsn_csm_contour_map(wks,rc3,res)

 ;;;打点
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  
  iplot=new(3,graphic)

  iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
  iplot(1)=gsn_csm_contour(wks,prob2,res_pro)
  iplot(2)=gsn_csm_contour(wks,prob3,res_pro)

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

  do i=0,2 
  iplot(i) = gsn_contour_shade(iplot(i), 0,0.1, opt)
end do
;;;  
 do i=0,2
overlay(plot(i), iplot(i))
end do
pres=True
 pres@gsnPanelLabelBar=True

gsn_panel(wks, plot, (/3,1/), pres)
end
