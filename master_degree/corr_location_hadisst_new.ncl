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
function  draw_corr(fils1,fils2,LeftString,wks1,wks2)
begin
    diri="/mnt/f/share/"
    locationxi_r=asciiread(diri+fils1, 55, "float")
    locationdong_r=asciiread(diri+fils2, 55, "float")   
    width_r=locationdong_r-locationxi_r
     
    ipindex1  =dim_rmvmean(locationxi_r)
    ipindex2 =dim_rmvmean(locationdong_r)
    

    ipindex3=width_r

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
  ssta=(/rmMonAnnCycTLL(sst(:,:,{100:300}))/)
  copy_VarCoords(sst(:,:,{100:300}), ssta)
  ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
  copy_VarCoords(ssta, ssta_de)
  ssta_son=month_to_season(ssta_de, "SON")
  ssta_djf=month_to_season(ssta_de, "DJF")

  ssta_s=(/ssta_son(0:54,:,:)+ssta_djf(1:,:,:)/)*0.5
  copy_VarCoords(ssta_son(0:54,:,:), ssta_s)

    delete(date)
    delete(timeind)
;;;

;;;
  r_xi_dong=escorc_n(ipindex1,ipindex2,0,0)
;;;相关性
  r_xi_sst=escorc_n(ipindex1,ssta_s,0,0)  
  r_dong_sst=escorc_n(ipindex2,ssta_s,0,0)
  copy_VarCoords(ssta_s(0,:,:), r_xi_sst)
  copy_VarCoords(ssta_s(0,:,:), r_dong_sst)

  r_xi_sst_pa=(/(r_xi_sst- r_xi_dong*r_dong_sst)/sqrt((1- r_xi_dong^2)*(1- r_dong_sst^2))/)
  r_dong_sst_pa=(/(r_dong_sst- r_xi_dong*r_xi_sst)/sqrt((1- r_xi_dong^2)*(1- r_xi_sst^2))/)
  r_cha=(/r_xi_sst_pa- r_dong_sst_pa/)

  copy_VarCoords(r_xi_sst, r_xi_sst_pa)
  copy_VarCoords(r_xi_sst, r_dong_sst_pa)
  copy_VarCoords(r_xi_sst, r_cha)
;;;
 
;;;偏相关的显著性检验
  prob_r_xi_sst_pa=rtest(r_xi_sst_pa,dimsizes(ssta_s&time)-1,0)
  prob_r_dong_sst_pa=rtest(r_dong_sst_pa,dimsizes(ssta_s&time)-1,0)
 copy_VarCoords(r_xi_sst, prob_r_xi_sst_pa)
 copy_VarCoords(r_xi_sst, prob_r_dong_sst_pa)
  prob_r_xi_sst_pa=where(prob_r_xi_sst_pa.eq.0,999 , prob_r_xi_sst_pa)
   prob_r_dong_sst_pa=where(prob_r_dong_sst_pa.eq.0,999 , prob_r_dong_sst_pa)


;;;
 

;;;画图
  gsn_define_colormap(wks1, "BlueWhiteOrangeRed")
  gsn_define_colormap(wks2, "BlueWhiteOrangeRed")
    ; gsn_define_colormap(wks, "BlRe ")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnLeftString=LeftString
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =True
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.7
  res@cnMaxLevelValF       =   0.7
  res@cnLevelSpacingF      =  0.1
  res@lbLabelBarOn=False

  res@lbOrientation = "Horizontal"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ssta&latitude)
  res@mpMaxLatF             = max(ssta&latitude)
  res@mpMinLonF             = min(ssta&longitude)
  res@mpMaxLonF             = max(ssta&longitude)
  res@mpCenterLonF          =avg(ssta&longitude)
  
  res@tmXTOn=False
  res@tmYROn=False

  ; res@vcRefMagnitudeF          = 3.0          ; define vector ref mag
  ; res@vcRefLengthF             = 0.04         ; define length of vec ref
  ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcMinDistanceF           = 0.008          ; thin vectors
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角
  res@mpShapeMode = "FreeAspect"
  res@vpWidthF = 0.7
  res@vpHeightF = 0.3

 plot=new(2,graphic)
  plot(0)= gsn_csm_contour_map(wks1,r_xi_sst_pa,res)

  plot(1)= gsn_csm_contour_map(wks2,r_dong_sst_pa,res)
 
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
  
  iplot=new(2,graphic)

  iplot(0)=gsn_csm_contour(wks1,prob_r_xi_sst_pa,res_pro)
  iplot(1)=gsn_csm_contour(wks2,prob_r_dong_sst_pa,res_pro)

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

  iplot(0) = gsn_contour_shade(iplot(0), 0,0.05, opt)
  iplot(1)= gsn_contour_shade(iplot(1), 0,0.05, opt)
;;;  
do i=0,1
overlay(plot(i), iplot(i))
end do
 return [/plot(0),plot(1)/]
end


;******************************************************************
;;主程序
;;**************************************************************
begin

fils1="20CR_msf_cold_locationxi_index.ascii"
fils2="NCPE1_msf_cold_locationxi_index.ascii"
fils3="ERA_msf_cold_locationxi_index.ascii"
fils4="JRA-55_msf_cold_locationxi_index.ascii"
fils5="MRE_msf_cold_locationxi_index.ascii"
; fils5="MRE_msf_mean_cold_locationxi_index.ascii"


fils6="20CR_msf_cold_locationdong_index.ascii"
fils7="NCPE1_msf_cold_locationdong_index.ascii"
fils8="ERA_msf_cold_locationdong_index.ascii"
fils9="JRA-55_msf_cold_locationdong_index.ascii"
fils10="MRE_msf_cold_locationdong_index.ascii"

plot1=new(4,graphic)
plot2=new(4,graphic)


wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个数据西位置指数hadisst偏相关")
wks2=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个数据东位置指数hadisst偏相关")

m1=draw_corr(fils1,fils6,"20CR",wks1,wks2)
m2=draw_corr(fils2,fils7,"NCPE1",wks1,wks2)
m3=draw_corr(fils3,fils8,"ERA",wks1,wks2)
m4=draw_corr(fils4,fils9,"JRA-55",wks1,wks2)

plot1(0)=m1[0]
plot1(1)=m2[0]
plot1(2)=m3[0]
plot1(3)=m4[0]

plot2(0)=m1[1]
plot2(1)=m2[1]
plot2(2)=m3[1]
plot2(3)=m4[1]

pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks1, plot1, (/2,2/), pres)
frame(wks1)

gsn_panel(wks2, plot2, (/2,2/), pres)
frame(wks2)


end
