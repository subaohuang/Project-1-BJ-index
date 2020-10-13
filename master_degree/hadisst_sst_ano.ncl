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
	fils2=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc") 
	f2=addfile(fils2, "r")
	date=cd_calendar(f2->time, 1)
	timeind=ind(date.le.201012.and.date.ge.195001)
  sst=lonFlip(f2->sst(timeind,{20:-20},:))
	delete(sst&time)
	sst&time=date(timeind)
  sst_avg=sst(:,:,:)
 

   
  ssta=rmMonAnnCycTLL(sst_avg) 
  printVarSummary(ssta)
  
   ssta=(/dtrend_msg_n(ispan(0,dimsizes(ssta&time)-1,1) ,ssta,False,1,0)/) ;;;;;看下去除长期趋势的情况

	ssta_50to78=dim_avg_n_Wrap(ssta({195001:197812},:,:),0)
	ssta_79to10=dim_avg_n_Wrap(ssta({197901:201012},:,:),0)
	ssta_50to10=dim_avg_n_Wrap(ssta({195001:201012},:,:),0)



;;;



;************************************************
;  画图
;************************************************

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/hadisst_异常场分阶段的情况")
  gsn_define_colormap(wks, "BlueWhiteOrangeRed")
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
  res@cnFillMode="RasterFill"
  res@cnRasterSmoothingOn=True
  
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  
  ; res@cnMinLevelValF       =   -0.2
  ; res@cnMaxLevelValF       =   0.2
  ; res@cnLevelSpacingF      =  0.04

  res@tmXTOn=False
  res@tmYROn=False

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ssta&latitude)
  res@mpMaxLatF             = max(ssta&latitude)
  res@mpMinLonF             = min(ssta&longitude)
  res@mpMaxLonF             = max(ssta&longitude)
  res@mpCenterLonF          =avg(ssta&longitude)
  
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.8
  res@vpHeightF = 0.2


	plot=new(3,graphic)
	plot(0)= gsn_csm_contour_map(wks,ssta_50to78,res)
    res@gsnRightString="1979-2010"
	plot(1)= gsn_csm_contour_map(wks,ssta_79to10,res)
    res@gsnRightString="1950-2010"
    ; res@cnMinLevelValF       =   -5.0
    ; res@cnMaxLevelValF       =   5.0
    ; res@cnLevelSpacingF      =  1.0

    plot(2)= gsn_csm_contour_map(wks,ssta_50to10,res)

;  ;;;打点
;   res_pro=True
;   res_pro@gsnDraw=False
;   res_pro@gsnFrame=False
;   res_pro@cnLevelSelectionMode = "ExplicitLevels"
;   res_pro@cnLinesOn=False
;   res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
;   res_pro@cnMinLevelValF      = 0.00        ; set min contour level
;   res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
;   res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

;   res_pro@cnLineLabelsOn=False
;   res_pro@tmYRMode = "Automatic"
;   res_pro@cnInfoLabelOn=False
  
;   iplot=new(3,graphic)

;   iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
;   iplot(1)=gsn_csm_contour(wks,prob2,res_pro)
;   iplot(2)=gsn_csm_contour(wks,prob3,res_pro)

;   opt                      = True
;   opt@gsnShadeFillType     = "pattern"
;   opt@gsnShadeMid     = 17         ; stipple
;   ; opt@gsnShadeHigh   = "green"
;   opt@gsnShadeFillScaleF   = 0.2         ; add extra density
;   opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

;   do i=0,2 
;   iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
; end do
; ;;;  
;  do i=0,2
; overlay(plot(i), iplot(i))
; end do
pres=True
 pres@gsnPanelLabelBar=True

gsn_panel(wks, plot, (/3,1/), pres)
end