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
ipindex_de=dim_standardize_n(dtrend(clm_ano(ipindex),True),0,0)
ipindex_de!0="time"
ipindex_de&time=ispan(1,732,1)

;;;读入HadiSST海温数据
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)

  timeind=ind(date.le.201012.and.date.ge.195001)
  sst=lonFlip(f1->sst(timeind,{60:-40},:))
  delete(sst&time)
  delete(timeind)
  delete(date)
  sst&time=ispan(1,732,1)
  


  ssta=dtrend_msg_n(ispan(0,dimsizes(sst&time)-1,1),rmMonAnnCycTLL(sst),\
                     False,False,0)

 
  ; ssta=dtrend_msg_n(ispan(0,dimsizes(sst&time)-1,1),runave_n_Wrap(rmMonAnnCycTLL(sst),3,0,0),\
  ;                    False,False,0)


   copy_VarCoords(sst, ssta)
;;;

;;;读入20cr的风场数据
  fils2=systemfunc("ls /mnt/f/share/ncep/20CR/uwnd.mon.mean.nc") 
  fils3=systemfunc("ls /mnt/f/share/ncep/20CR/vwnd.mon.mean.nc") 
  f2=addfile(fils2, "r")
  f3=addfile(fils3, "r")
  date=cd_calendar(f2->time, 1)
  timeind=ind(date.le.201012.and.date.ge.195001)
  u=f2->uwnd(timeind,{1000},{30:-30},:)
  v=f3->vwnd(timeind,{1000},{30:-30},:)
  delete(timeind)
  u&time=ispan(1,732,1)
  v&time=ispan(1,732,1)
  
 
  ua=dtrend_msg_n(ispan(0,dimsizes(u&time)-1,1),rmMonAnnCycTLL(u),\
                     False,False,0)
  copy_VarCoords(u, ua)


  va=dtrend_msg_n(ispan(0,dimsizes(v&time)-1,1),rmMonAnnCycTLL(v),\
                     False,False,0)
  copy_VarCoords(v, va)
;;;


;;;回归u场,v场

  ssta_rc=regCoef_n(ipindex_de,ssta,0,0)  
  copy_VarCoords(ssta(0,:,:), ssta_rc)

  ua_rc=regCoef_n(ipindex_de,ua,0,0)  
  copy_VarCoords(ua(0,:,:), ua_rc)

  va_rc=regCoef_n(ipindex_de,va,0,0)  
  copy_VarCoords(va(0,:,:), va_rc)
;;;



;;;显著性检验
  prob_ssta=rc_test(ssta_rc)
  prob_ua=rc_test(ua_rc)
  prob_va=rc_test(va_rc)
;;;


;;;风场的检验

  ua_rc=where(prob_ua.le.0.05, ua_rc, ua_rc@_FillValue)
  va_rc=where(prob_va.le.0.05, va_rc, va_rc@_FillValue)
;;;





;;;画图
  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/ip指数回归的哈得来sst和20cr风场")
  gsn_define_colormap(wks, "BlRe")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnScalarContour        = True
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =True
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.50
  res@cnMaxLevelValF       =   0.50
  res@cnLevelSpacingF      =  0.1
  

  ; res@cnMinLevelValF       =   -0.02
  ; res@cnMaxLevelValF       =   0.02
  ; res@cnLevelSpacingF      =  0.005

  res@lbOrientation = "Horizontal"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ssta&latitude)
  res@mpMaxLatF             = max(ssta&latitude)
  res@mpMinLonF             = min(ssta&longitude)
  res@mpMaxLonF             = max(ssta&longitude)
  res@mpCenterLonF          =avg(ssta&longitude)
  
  res@tmXTOn=False
  res@tmYROn=False

  res@vcRefMagnitudeF          = 1.0           ; define vector ref mag
  res@vcRefLengthF             = 0.04         ; define length of vec ref
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  res@vcMinDistanceF           = 0.008          ; thin vectors
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.9
  res@vpHeightF = 0.2


 plot= gsn_csm_vector_scalar_map(wks,ua_rc,va_rc,ssta_rc,res)
 


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
 
  iplot=gsn_csm_contour(wks,prob_ssta,res_pro)

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