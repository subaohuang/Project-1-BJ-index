;;;
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
;;;
function clm_ano(x[*]:numeric)  ;去除一维时间序列的年循环
begin

   dimx  = dimsizes (x)
   ntim  = dimx(0)

   nmos  = 12
   modCheck ("clm", ntim, nmos)   ; error check; Compute all 12 monthly averages first. 
   aveMonth = new((/nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do
   copy_VarAtts (x,aveMonth)
   aveMonth!0   = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var
   xAnom=x
   do yr=0,ntim-1,nmos
      xAnom(yr:yr+nmos-1) = (/ x(yr:yr+nmos-1)- aveMonth /)
   end do
   return (xAnom)
end
begin

;;;位置指数的回归
  ; diri1="/mnt/f/share/era/data/"
  ; fils2=systemfunc("ls "+diri1+"msf_era_location1979-2017.nc")
  ; f2=addfile(fils2, "r")
  ; ipindex=f2->locationxi_a(0:455)
  ; ipindex_de=dtrend_msg_n(ispan(0,455,1),clm_ano(ipindex),\
  ;                    False,False,0)
;;;
ipindex=asciiread("/mnt/f/share/era/data/ip_index_1950_2010.ascii", 732, "float")
; ipindex_de=dtrend(clm_ano(ipindex),True)
ipindex_de=dim_standardize_n(dtrend(clm_ano(ipindex),True),0,0)
ipindex_de!0="time"
ipindex_de&time=ispan(1,732,1)



;;;读入20cr的风场数据
  fils=systemfunc("ls /mnt/f/share/ncep/20CR/omega.mon.mean.nc") 
  f1=addfile(fils, "r")

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201012.and.date.ge.195001)
  w_r=f1->omega(timeind,:,{30:-30},:)
  delete(timeind)
  w_r&time=ispan(1,732,1)
  w=dim_avg_n_Wrap(w_r(:,{850:200},:,:),1)
  
  w=(/w*-100/)

 
  wa=dtrend_msg_n(ispan(0,dimsizes(w&time)-1,1),rmMonAnnCycTLL(w),\
                     False,False,0)
  ; wa=dtrend_msg_n(ispan(0,dimsizes(w&time)-1,1),runave_n_Wrap(rmMonAnnCycTLL(w),3,0,0),\
  ;                    False,False,0)
  
  copy_VarCoords(w, wa)

;;;
 

;;;回归w场

  wa_rc=regCoef_n(ipindex_de,wa,0,0)  
  copy_VarCoords(wa(0,:,:), wa_rc)

;;;



;;;显著性检验
  prob_wa=rc_test(wa_rc)
;;;

;;;画图
  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/ip指数回归的20cr垂直速度异常")
  gsn_define_colormap(wks, "BrownBlue12")
  gsn_define_colormap(wks, "ViBlGrWhYeOrRe")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
 
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnFillMode="RasterFill"
  res@cnRasterSmoothingOn=True
  
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2

  ; res@cnMinLevelValF       =   -0.3
  ; res@cnMaxLevelValF       =   0.3
  ; res@cnLevelSpacingF      =  0.1


  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(wa&lat)
  res@mpMaxLatF             = max(wa&lat)
  res@mpMinLonF             = min(wa&lon)
  res@mpMaxLonF             = max(wa&lon)
  res@mpCenterLonF          =avg(wa&lon)

  res@tmXTOn=False
  res@tmYROn=False

  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.9
  res@vpHeightF = 0.2
 
 plot= gsn_csm_contour_map(wks,wa_rc,res)
 


;;;打点
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True
 
  iplot=gsn_csm_contour(wks,prob_wa,res_pro)

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
  iplot = gsn_contour_shade(iplot, 0,0.05, opt)
;;;  
overlay(plot, iplot)
 draw(plot)
 frame(wks)


end