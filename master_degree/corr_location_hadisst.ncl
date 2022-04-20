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
;;;;
function  draw_corr(mode,wks)
begin
; ipindex1=asciiread("/mnt/f/share/MRE_msf_mean_cold_locationxi_index.ascii", 55, "float")
; ipindex2=asciiread("/mnt/f/share/MRE_msf_mean_cold_locationdong_index.ascii", 55, "float")

ipindex1=asciiread("/mnt/f/share/JRA-55_msf_cold_locationxi_index.ascii", 55, "float")
ipindex2=asciiread("/mnt/f/share/JRA-55_msf_cold_locationdong_index.ascii", 55, "float")

ipindex3=ipindex2- ipindex1


ipindex1=dim_rmvmean(ipindex1)
ipindex2=dim_rmvmean(ipindex2)
ipindex3=dim_rmvmean(ipindex3)

;;;读入HadiSST海温数据
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)

  timeind=ind(date.le.201312.and.date.ge.195801)
  

  sst=lonFlip(f1->sst(timeind,{30:-30},:))
  ssta=(/rmMonAnnCycTLL(sst)/)
  copy_VarCoords(sst, ssta)
  ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
  copy_VarCoords(ssta, ssta_de)
  ssta_son=month_to_season(ssta_de, "SON")
  ssta_djf=month_to_season(ssta_de, "DJF")

  ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
  copy_VarCoords(ssta_son(0:54,:,:), ssta_s)

;;;
delete(date)
delete(timeind)
;;;读入20cr的风场数据
  fils2=systemfunc("ls /mnt/f/share/ncep/20CR/uwnd.mon.mean.nc") 
  fils3=systemfunc("ls /mnt/f/share/ncep/20CR/vwnd.mon.mean.nc") 
  f2=addfile(fils2, "r")
  f3=addfile(fils3, "r")
  date=cd_calendar(f2->time, 1)
  timeind=ind(date.le.201212.and.date.ge.195801)
  u=f2->uwnd(timeind,{1000},{30:-30},:)
  v=f3->vwnd(timeind,{1000},{30:-30},:)
  delete(u&time)
  delete(v&time)
  u&time=date(timeind)
  v&time=date(timeind)
 printVarSummary(u)

  
  ua=dtrend_msg_n(ispan(0,dimsizes(u&time)-1,1),rmMonAnnCycTLL(u),\
                     False,False,0)
  copy_VarCoords(u, ua)


  va=dtrend_msg_n(ispan(0,dimsizes(v&time)-1,1),rmMonAnnCycTLL(v),\
                     False,False,0)
  copy_VarCoords(v, va)


  ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 ua_son=month_to_season(ua(time|:,lat|:,lon|:), "SON")
 ua_djf=month_to_season(ua(time|:,lat|:,lon|:), "DJF")

 ua_son&time=ispan(1958,2012,1)
 ua_djf&time=ispan(1958,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
  printVarSummary(ua)
    ua_djf(0:dimsizes(ua_djf&time)-2,:,:)=(/ua_djf(1:,:,:)/)
    ua_djf(dimsizes(ua_djf&time)-1,:,:)=(/ua({201212},:,:)/)
 ;;;

 va_son=month_to_season(va(time|:,lat|:,lon|:), "SON")
 va_djf=month_to_season(va(time|:,lat|:,lon|:), "DJF")

 va_son&time=ispan(1958,2012,1)
 va_djf&time=ispan(1958,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ;;注意点
    va_djf(0:dimsizes(va_djf&time)-2,:,:)=(/va_djf(1:,:,:)/)
    va_djf(dimsizes(va_djf&time)-1,:,:)=(/va({201212},:,:)/)
 ;;;
 
 ua_cold=ua_son
 ua_cold=(/(ua_son(:,:,:)+ua_djf(:,:,:))*0.5/)
 
 va_cold=va_son
 va_cold=(/(va_son(:,:,:)+va_djf(:,:,:))*0.5/)

; ;;;回归u场,v场

;   ssta_rc=regCoef_n(ipindex2,ssta_s,0,0)  
;   copy_VarCoords(ssta_s(0,:,:), ssta_rc)

;   ua_rc=regCoef_n(ipindex2,ua_cold,0,0)  
;   copy_VarCoords(ua(0,:,:), ua_rc)

;   va_rc=regCoef_n(ipindex2,va_cold,0,0)  
;   copy_VarCoords(va(0,:,:), va_rc)
; ;;;显著性检验
;   prob_ssta=rc_test(ssta_rc)
;   prob_ua=rc_test(ua_rc)
;   prob_va=rc_test(va_rc)
; ;;;
; ;;;

if mode.eq.0 then 
  ipindex=ipindex1
 end if 
 if mode.eq.1 then 
  ipindex=ipindex2
 end if 

if mode.eq.2 then 
  ipindex=ipindex3
 end if 
  ip_corr=escorc_n(ipindex1,ipindex2,0,0)
;;;相关性
  ssta_rc=escorc_n(ipindex,ssta_s,0,0)  
  copy_VarCoords(ssta_s(0,:,:), ssta_rc)

  ua_rc=escorc_n(ipindex,ua_cold,0,0)  
  copy_VarCoords(ua(0,:,:), ua_rc)

  va_rc=escorc_n(ipindex,va_cold,0,0)  
  copy_VarCoords(va(0,:,:), va_rc)
;;;

;;;计算sst和xi index以及东index的偏相关性



;;;显著性检验
  prob_ssta=rtest(ssta_rc,dimsizes(ssta_s&time),0)
  prob_ua=rtest(ua_rc,dimsizes(ua_cold&time),0)
  prob_va=rtest(va_rc,dimsizes(va_cold&time),0)
  prob_ssta=where(prob_ssta.eq.0, prob_ssta@_FillValue, prob_ssta)
  copy_VarCoords(ssta_rc, prob_ssta)
  copy_VarCoords(ua_rc, prob_ua)
  copy_VarCoords(va_rc, prob_va)
;;;



;;;风场的检验
  ua_rc=where(prob_ua.le.0.05, ua_rc, ua_rc@_FillValue)
  va_rc=where(prob_va.le.0.05, va_rc, va_rc@_FillValue)
  printVarSummary(ua_rc)
  printVarSummary(va_rc)
;;;



;;;画图
  gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    ; gsn_define_colormap(wks, "BlRe ")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnScalarContour        = True
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =True
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.7
  res@cnMaxLevelValF       =   0.7
  res@cnLevelSpacingF      =  0.1
  res@lbLabelBarOn=False
  if mode.eq.0 then
  res@gsnLeftString="Western edge&Hadisst"
  end if
   if mode.eq.1 then
  res@gsnLeftString="Eastern edge&Hadisst"
  end if
   if mode.eq.2 then
  res@gsnLeftString="Width&Hadisst"
  end if

  ; res@cnMinLevelValF       =   -0.1
  ; res@cnMaxLevelValF       =   0.1
  ; res@cnLevelSpacingF      =  0.01
  
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

  res@vcRefMagnitudeF          = 3.0          ; define vector ref mag
  res@vcRefLengthF             = 0.04         ; define length of vec ref
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  res@vcMinDistanceF           = 0.008          ; thin vectors
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.8
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
 return[/plot,res,ssta_rc,ip_corr,res_pro,opt/]
end
;;;主程序
begin
plot=new(4,graphic)
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/MRE的位置指数相关hadisst")
m1=draw_corr(0,wks)
m2=draw_corr(1,wks)
m3=draw_corr(2,wks)

plot(0)=m1[0]
plot(1)=m2[0]
plot(3)=m3[0]
printVarSummary(plot(0))


ssta_rc=m1[2]
ssta_rc=(/m1[2]-m2[2]/)
res=m1[1]
; delete(res@cnLevelSelectionMode)
; delete(res@cnMinLevelValF )
; delete(res@cnMaxLevelValF )
; delete(res@cnLevelSpacingF)
res@gsnLeftString="corr West- corr East"
plot(2)=gsn_csm_contour_map(wks,ssta_rc,res)

pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks, plot, (/2,2/), pres)

;;;画偏相关的图
 r_xi_sst=m1[2]
 r_dong_sst=m2[2]
 r_xi_dong=m1[3]
 printVarSummary(r_xi_sst)
;;
r_xi_sst_pa=r_xi_sst
r_dong_sst_pa=r_xi_sst
r_cha=r_xi_sst

r_xi_sst_pa=(/(r_xi_sst- r_xi_dong*r_dong_sst)/sqrt((1- r_xi_dong^2)*(1- r_dong_sst^2))/)
r_dong_sst_pa=(/(r_dong_sst- r_xi_dong*r_xi_sst)/sqrt((1- r_xi_dong^2)*(1- r_xi_sst^2))/)
r_cha=(/r_xi_sst_pa- r_dong_sst_pa/)

plot_pa=new(3,graphic)
res_pa=m1[1]
res_pa@cnFillPalette="BlueWhiteOrangeRed"
res_pa@gsnLeftString="Western edge&Hadisst"
res_pa@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
res_pa@cnMinLevelValF       =   -0.7
res_pa@cnMaxLevelValF       =   0.7
res_pa@cnLevelSpacingF      =  0.1
wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/JRA55的位置指数偏相关系数hadisst")
plot_pa(0)=gsn_csm_contour_map(wks1,r_xi_sst_pa,res_pa)
res_pa@gsnLeftString="Eastern edge&Hadisst"
plot_pa(1)=gsn_csm_contour_map(wks1,r_dong_sst_pa,res_pa)
res_pa@gsnLeftString="corr of Western - easten"
plot_pa(2)=gsn_csm_contour_map(wks1,r_cha,res_pa)



gsn_panel(wks1, plot_pa, (/3,1/), pres)
frame(wks1)


end
