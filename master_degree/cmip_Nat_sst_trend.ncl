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

   model="GFDL-CM3"
   yearstart=195801
   yearend=200512
   diri_output="/mnt/f/share/CMIP5/msf_history_natureforce/"
	
 
    fils1=systemfunc("ls /mnt/j/CMIP5/history_natureforce/sst/"+model+"/r1i1p1/tos_Omon_"+model+"_historicalNat_r1i1p1"+"*"+".nc")
    
    print(fils1)
 
    
    f1=addfiles(fils1, "join")



    time=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date=cd_calendar(time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)


     sst=f1[:]->tos(timeind,{30:-30},:)

  sst=(/sst-273.15/)
  printVarSummary(sst)
;;;

;************************************************
; 计算回归
;************************************************
	rc1 =regCoef_n(ispan(1, dimsizes(sst&time), 1),sst,0,0)
	copy_VarCoords(sst(0,:,:), rc1)
  printVarSummary(rc1)

;************************************************
; 显著性检验
;************************************************
	prob1=rc_test(rc1)
  rc1=(/rc1*120/) 
  copy_VarCoords(sst(1,:,:), rc1)

;************************************************
;  画图
;************************************************

  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/"+model+"_sst趋势变化")
   gsn_define_colormap(wks, "GMT_polar")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=True

  
  res@cnLineColors="Transparent"  ;
  res@lbLabelBarOn=False
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  res@cnFillMode="RasterFill"
  res@cnRasterSmoothingOn=True
  
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  
  ; res@cnMinLevelValF       =   -0.4
  ; res@cnMaxLevelValF       =   0.4
  ; res@cnLevelSpacingF      =  0.1

  res@tmXTOn=False
  res@tmYROn=False

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30.0
  res@mpMaxLatF             = 30.0
  res@mpMinLonF             = -180.0
  res@mpMaxLonF             = 180.0
  res@mpCenterLonF          =180.0
  
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.8
  res@vpHeightF = 0.2
  
  plot=new(1,graphic)
	plot(0)= gsn_csm_contour_map(wks,rc1,res)

 ;;;打点
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.5       ; set contour spacing

  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  
  iplot=new(1,graphic)

  iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
 

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger


  iplot(0) = gsn_contour_shade(iplot(0), 0,0.5, opt)

;;;  

overlay(plot(0), iplot(0))

draw(plot)
frame(wks)
end