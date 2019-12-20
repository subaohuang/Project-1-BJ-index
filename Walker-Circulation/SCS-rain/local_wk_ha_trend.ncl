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
    diri="/home/ys17-19/lsh/data/era/"
    fils1=systemfunc("ls "+diri+"era_m_zonal_meri_new.nc")
    f1=addfile(fils1, "r")
    m_zonal=f1->m_zonal
    m_meri=f1->m_meri
    
    m_zonal_clm = dim_avg_n_Wrap(m_zonal(:,{500},:,:), 0)
    m_meri_clm = dim_avg_n_Wrap(m_meri(:,{500},:,:), 0)

    m_zonal_500=rmMonAnnCycTLL(m_zonal(:,{500},{-60:60},:))
    m_meri_500=rmMonAnnCycTLL(m_meri(:,{500},{-60:60},:))
    
    m_zonal_annual=month_to_annual(m_zonal_500, 1)
    m_meri_annual=month_to_annual(m_meri_500, 1)
    
    dims=dimsizes(m_zonal_annual)
;;

 ;;;回归p场

  m_zonal_rc1=regCoef_n(ispan(1,dims(0),1),m_zonal_annual,0,0)  
  copy_VarCoords(m_zonal_annual(0,:,:), m_zonal_rc1)
  m_zonal_rc1=(/m_zonal_rc1*10/)
  m_zonal_rc2=regCoef_n(ispan(1,dims(0),1),m_meri_annual,0,0)  
  copy_VarCoords(m_zonal_annual(0,:,:), m_zonal_rc2)
  m_zonal_rc2=(/m_zonal_rc2*10/)

print("趋势计算结束")


;;;显著性检验
  prob_m_zonal1=rc_test(m_zonal_rc1)
  prob_m_zonal2=rc_test(m_zonal_rc2)
;;;
print("检验结束")


;;;画图 
  wks=gsn_open_wks("pdf", "ERA-in的local沃克及哈得来环流的趋势")

  ; gsn_define_colormap(wks, "BlRe")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.0008
  res@cnMaxLevelValF       =   0.0008
  res@cnLevelSpacingF      =  0.0001
  res@gsnLeftString = "local WC trend"
  res@cnFillPalette = "BlRe"
  ; res@cnMinLevelValF       =   -0.02
  ; res@cnMaxLevelValF       =   0.02
  ; res@cnLevelSpacingF      =  0.005

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(m_zonal_annual&lat)
  res@mpMaxLatF             = max(m_zonal_annual&lat)
  res@mpMinLonF             = min(m_zonal_annual&lon)
  res@mpMaxLonF             = max(m_zonal_annual&lon)
  res@mpCenterLonF          =avg(m_zonal_annual&lon)
   

  res@tmXTOn=False
  res@tmYROn=False

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpWidthF = 0.9
  ; res@vpHeightF = 0.2 
  ;;; *************************************************
  plot=new(2,graphic)

  plot(0)= gsn_csm_contour_map(wks,m_zonal_rc1,res)
  res@gsnLeftString = "local HC trend"
  plot(1)= gsn_csm_contour_map(wks,m_zonal_rc2,res)
  ;*********************************************************
   ;plot clm
     wks1=gsn_open_wks("pdf", "ERA-in的local沃克及哈得来环流的气候态")
  plot_clm = new(2,graphic)

  res@cnMinLevelValF       =   -0.005
  res@cnMaxLevelValF       =   0.005
  res@cnLevelSpacingF      =  0.0005
  res@cnFillPalette = "NCV_blue_red"
    res@gsnLeftString = "local WC climatology"
  plot_clm(0) = gsn_csm_contour_map(wks1,m_zonal_clm, res)
  res@gsnLeftString = "local HC climatology"
  plot_clm(1) = gsn_csm_contour_map(wks1,m_meri_clm ,res)
   pres=True
  gsn_panel(wks1, plot_clm, (/2,1/), pres)

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

 
  iplot=new(2,graphic)

  iplot(0)=gsn_csm_contour(wks,prob_m_zonal1,res_pro)
  iplot(1)=gsn_csm_contour(wks,prob_m_zonal2,res_pro)
  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 2.         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.005      ; make dots larger

  do i=0,1
  iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
end do
;;;  
  do i=0,1
  overlay(plot(i), iplot(i))
  end do
  pres=True
  gsn_panel(wks, plot, (/2,1/), pres)

end 

