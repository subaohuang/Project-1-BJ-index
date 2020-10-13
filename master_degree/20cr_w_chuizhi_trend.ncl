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
  fils2=systemfunc("ls /mnt/f/share/ncep/ncep1/omega.mon.mean.nc") 
  f2=addfile(fils2, "r")
  date=cd_calendar(f2->time, 1)
  timeind=ind(date.le.201012.and.date.ge.195001)
  w=f2->omega(timeind,:,{5:-5},:)
  delete(w&time)
  w&time=date(timeind)
  w=(/w*-1000/)
  
  w_avg=dim_avg_n_Wrap(w(:,:,{5:-5},:), 2)

  w_clm=dim_avg_n_Wrap(w_avg, 0)

  w_50to78=w_avg({195001:197812},:,:)
  w_79to10=w_avg({197901:201012},:,:)
  w_50to10=w_avg({195001:201012},:,:)

  printVarSummary(w_50to78)
 
;;;

;************************************************
; 计算回归
;************************************************
  rc1 =regCoef_n(ispan(1, dimsizes(w_50to78&time), 1),w_50to78,0,0)
  rc2 =regCoef_n(ispan(1, dimsizes(w_79to10&time), 1),w_79to10,0,0)
  rc3 =regCoef_n(ispan(1, dimsizes(w_50to10&time), 1),w_50to10,0,0)
     
   
  copy_VarCoords(w_50to78(0,:,:), rc1)
  copy_VarCoords(w_50to78(0,:,:), rc2)
  copy_VarCoords(w_50to78(0,:,:), rc3)

;************************************************
; 显著性检验
;************************************************
  prob1=rc_test(rc1)
  prob2=rc_test(rc2)
  prob3=rc_test(rc3)

    rc1=(/rc1*120/) 
    rc2=(/rc2*120/) 
    rc3=(/rc3*120/) 
 

;************************************************
; 画图
;************************************************
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/垂直速度的趋势&垂直分布图")
  gsn_define_colormap(wks,"GMT_polar")
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="w trend"
    res2@gsnRightString="1950-1978"
   
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -6.0
    res2@cnMaxLevelValF = 6.0
    res2@cnLevelSpacingF = 1.0
    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
  
    res2@lbOrientation="Horizontal"
    res2@lbLabelBarOn=False

    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    res2@tmXBValues=(/0,90,180,270,357.5/)
    res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(0,360,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7
   
   plot=new(3,graphic)

    plot(0)=gsn_csm_pres_hgt(wks,rc1,res2)
     res2@gsnRightString="1979-2010"
    plot(1)=gsn_csm_pres_hgt(wks,rc2,res2)
     res2@gsnRightString="1950-2010"
    plot(2)=gsn_csm_pres_hgt(wks,rc3,res2)
   res2@gsnRightString=""
;;;gsn_contour_shade
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@gsnRightString=""
  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False

   iplot=new(3,graphic)
   iplot(0)=gsn_csm_pres_hgt(wks,prob1,res_pro)
   iplot(1)=gsn_csm_pres_hgt(wks,prob2,res_pro)
   iplot(2)=gsn_csm_pres_hgt(wks,prob3,res_pro)
  

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

  iplot(0) = gsn_contour_shade(iplot(0), 0,0.05, opt)
  iplot(1) = gsn_contour_shade(iplot(1), 0,0.05, opt)
  iplot(2) = gsn_contour_shade(iplot(2), 0,0.05, opt)
 ;;;

 res2@cnLevelSelectionMode = "ExplicitLevels"
 res2@cnLinesOn=True
 res2@cnLevels = ispan(-50,50,10)
 res2@cnLineLabelsOn=True
 res2@cnLineLabelBackgroundColor="Transparent"
 res2@cnLineThicknessF = 2.
 res2@cnFillOn        = False
 res2@gsnContourNegLineDashPattern=1
 plot1=new(3,graphic)
 do i=0,2
 plot1(i) = gsn_csm_pres_hgt(wks, w_clm, res2)
 end do
;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 15.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=15
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 25)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
      
  ployres@gsLineDashPattern=0
  dum1=new(3,graphic)
  dum2=dum1
  dum3=dum1
  do i=0,2
  dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
  ployres@gsLineDashPattern=16
  dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
  ployres@gsLineDashPattern=0
  dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
 end do 

do i=0,2
 overlay(plot(i), iplot(i))
 overlay(plot(i),plot1(i))
end do
pres=True
pres@gsnPanelLabelBar=True
pres@pmLabelBarWidthF=0.8
gsn_panel(wks, plot, (/3,1/), pres)
end