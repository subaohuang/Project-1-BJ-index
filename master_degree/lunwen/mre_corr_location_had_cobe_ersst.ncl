;;;
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 3
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end
;;;

;;;;
function  draw_corr(rea,wks1,wks2,wks3,LeftString)
begin
    diri="/mnt/f/share/rea_msf/"
    locationxi_r=asciiread(diri+rea+"_msf_djf_locationxi_node_index.ascii", 34, "float")
    locationdong_r=asciiread(diri+rea+"_msf_djf_locationdong_node_index.ascii", 34, "float")   
    width_r=locationdong_r-locationxi_r
     

    ipindex1  =dim_rmvmean(locationxi_r)
    ipindex2 =dim_rmvmean(locationdong_r)
    
    ipindex3=width_r
    ipindex3=dim_rmvmean(ipindex3)

  ipindex1=(/dtrend_msg_n(ispan(1,34,1),ipindex1,False,True,0)/)
  ipindex2=(/dtrend_msg_n(ispan(1,34,1),ipindex2,False,True,0)/)
  ipindex3=(/dtrend_msg_n(ispan(1,34,1),ipindex3,False,True,0)/)

;;;读入海温数据
  if  rea.eq."JRA55" then
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"COBE_SST.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201312.and.date.ge.197901)
  sst=f1->sst(timeind,{30:-30},:)
  sst!0="time"
  sst!1="lat"
  sst!2="lon"
  delete(date)
  delete(timeind)
  delete(f1)

  
  else if rea.eq."20CR".or.rea.eq."ERA"  then 

  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201312.and.date.ge.197901)
  sst=lonFlip(f1->sst(timeind,{30:-30},:))
  sst!0="time"
  sst!1="lat"
  sst!2="lon"
  delete(date)
  delete(timeind)
  delete(f1)

   else 
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201312.and.date.ge.197901)
  sst=f1->sst(timeind,{30:-30},:)
  sst!0="time"
  sst!1="lat"
  sst!2="lon"
  delete(date)
  delete(f1)

   end if 
 end if

;;;处理成异常的冬季平均场
  ssta=(/rmMonAnnCycTLL(sst(:,:,{100:300}))/)
  copy_VarCoords(sst(:,:,{100:300}), ssta)
  ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
  copy_VarCoords(ssta, ssta_de)
  ssta_djf=month_to_season(ssta_de, "DJF")
  ssta_son=month_to_season(ssta_de, "SON")


  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
  ; ssta_s=(/ssta_son(0:33,:,:)/)
  copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
;;;


 ; nino4=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))
 nino4=dim_avg_n(ssta_s(:,{5:-5},{160:210}), (/1,2/))
    
 nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
;;;
 
  print("*********"+rea+"***********")
   r_xi_dong=escorc_n(ipindex1,ipindex2,0,0)
  print("r_xi_dong="+r_xi_dong)
  
  r_width_xi=escorc_n(ipindex3, ipindex1,0,0)
  r_width_dong=escorc_n(ipindex3, ipindex2,0,0)
  
  print("r_width_xi="+r_width_xi)
  print("r_width_dong="+r_width_dong)

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
  
   siglvl = 0.05
  r_nino3_xi=escorc_n(ipindex1, nino3, 0, 0)
            ; a-priori specified sig level
  print("r_nino3_xi="+r_nino3_xi)
   pr     = rtest(r_nino3_xi, dimsizes(ipindex1), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)

  r_nino3_dong=escorc_n(ipindex2, nino3, 0, 0)
  print("r_nino3_dong="+r_nino3_dong)
   pr     = rtest(r_nino3_dong, dimsizes(ipindex1), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)

  r_nino3_width=escorc_n(ipindex3, nino3, 0, 0)
  print("r_nino3_width="+r_nino3_width)
   pr     = rtest(r_nino3_width, dimsizes(ipindex1), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)

  r_nino4_xi=escorc_n(ipindex1, nino4, 0, 0)
  print("r_nino4_xi="+r_nino4_xi)
    pr     = rtest(r_nino4_xi, dimsizes(ipindex1), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)
  r_nino4_dong=escorc_n(ipindex2, nino4, 0, 0)
  print("r_nino4_dong="+r_nino4_dong)
   pr     = rtest(r_nino4_dong, dimsizes(ipindex1), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)

  r_nino4_width=escorc_n(ipindex3, nino4, 0, 0)
  print("r_nino4_width="+r_nino4_width)
     pr     = rtest(r_nino4_width, dimsizes(ipindex1), 0)    ; ===> pr=0.046
  if (pr.lt.siglvl) then
      print("显著检验=1")
  else
      print("显著检验=0")
  end if
  delete(pr)
  
  r_nino4_xi_pa=(/(r_nino4_xi- r_xi_dong*r_nino4_dong)/sqrt((1- r_xi_dong^2)*(1- r_nino4_dong^2))/)
  r_nino4_dong_pa=(/(r_nino4_dong- r_xi_dong*r_nino4_xi)/sqrt((1- r_xi_dong^2)*(1- r_nino4_xi^2))/)

  r_nino3_xi_pa=(/(r_nino3_xi- r_xi_dong*r_nino3_dong)/sqrt((1- r_xi_dong^2)*(1- r_nino3_dong^2))/)
  r_nino3_dong_pa=(/(r_nino3_dong- r_xi_dong*r_nino3_xi)/sqrt((1- r_xi_dong^2)*(1- r_nino3_xi^2))/)
  

  print("r_nino4_xi_pa="+r_nino4_xi_pa)
  print("r_nino4_dong_pa="+r_nino4_dong_pa)
  print("r_nino3_xi_pa="+r_nino3_xi_pa)
  print("r_nino3_dong_pa="+r_nino3_dong_pa)


;;;
 
;;;偏相关的显著性检验
  prob_r_xi_sst_pa=rtest(r_xi_sst_pa,dimsizes(ssta_s&time)-1,0)
  prob_r_dong_sst_pa=rtest(r_dong_sst_pa,dimsizes(ssta_s&time)-1,0)
  copy_VarCoords(r_xi_sst, prob_r_xi_sst_pa)
  copy_VarCoords(r_xi_sst, prob_r_dong_sst_pa)
  prob_r_xi_sst_pa=where(prob_r_xi_sst_pa.eq.0,999 , prob_r_xi_sst_pa)
  prob_r_dong_sst_pa=where(prob_r_dong_sst_pa.eq.0,999 , prob_r_dong_sst_pa)
;;;

  r_width_sst=escorc_n(ipindex3,ssta_s,0,0)  
  copy_VarCoords(ssta_s(0,:,:), r_width_sst)

  prob_width_sst=rtest(r_width_sst,dimsizes(ssta_s&time),0)
  copy_VarCoords(r_width_sst, prob_width_sst)
  prob_width_sst=where(prob_width_sst.eq.0,999 , prob_width_sst)



;;;画图
  gsn_define_colormap(wks1, "BlueWhiteOrangeRed")
  gsn_define_colormap(wks2, "BlueWhiteOrangeRed")
   gsn_define_colormap(wks3, "BlueWhiteOrangeRed")
    ; gsn_define_colormap(wks, "BlRe ")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnLeftString=LeftString
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  res@lbLabelBarOn=False

  res@lbOrientation = "Horizontal"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ssta&lat)
  res@mpMaxLatF             = max(ssta&lat)
  res@mpMinLonF             = min(ssta&lon)
  res@mpMaxLonF             = max(ssta&lon)
  res@mpCenterLonF          =avg(ssta&lon)
  res@mpShapeMode   ="FreeAspect"
  res@vpHeightF=0.35
  res@vpWidthF=0.65

  res@tmXTOn=False
  res@tmYROn=False
  


  ; res@vcRefMagnitudeF          = 3.0          ; define vector ref mag
  ; res@vcRefLengthF             = 0.04         ; define length of vec ref
  ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcMinDistanceF           = 0.008          ; thin vectors
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角
  ; res@mpShapeMode = "FreeAspect"
  ; res@vpWidthF = 0.8
  ; res@vpHeightF = 0.2

 plot=new(3,graphic)
 plot_zero=new(3,graphic)
  plot(0)= gsn_csm_contour_map(wks1,r_xi_sst_pa,res)
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.1
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2
  
  plot(1)= gsn_csm_contour_map(wks2,r_dong_sst_pa,res)
    res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  plot(2)=gsn_csm_contour_map(wks3, r_width_sst, res)

  res_z=True
  res_z@gsnDraw=False
  res_z@gsnFrame=False
  res_z@cnFillOn=False
  res_z@cnLinesOn=True
  res_z@cnLevelSelectionMode="ExplicitLevels"
  res_z@cnLevels=0
  res_z@cnLineThicknessF=3.0
  res_z@cnInfoLabelOn=False
  res_z@cnLineLabelsOn=False

  plot_zero(0)=gsn_csm_contour(wks1,r_xi_sst_pa,res_z)
  plot_zero(1) =gsn_csm_contour(wks2,r_dong_sst_pa,res_z)
  plot_zero(2) =gsn_csm_contour(wks3,r_width_sst,res_z)
overlay(plot(0), plot_zero(0))
overlay(plot(1), plot_zero(1))
overlay(plot(2), plot_zero(2))


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
  
  iplot=new(3,graphic)

  iplot(0)=gsn_csm_contour(wks1,prob_r_xi_sst_pa,res_pro)
  iplot(1)=gsn_csm_contour(wks2,prob_r_dong_sst_pa,res_pro)
  iplot(2)=gsn_csm_contour(wks3,prob_width_sst,res_pro)
  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadnino4d     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

  iplot(0) = gsn_contour_shade(iplot(0), 0,0.05, opt)
  iplot(1)= gsn_contour_shade(iplot(1), 0,0.05, opt)
  iplot(2)= gsn_contour_shade(iplot(2), 0,0.05, opt)
;;;  
do i=0,2
overlay(plot(i), iplot(i))
end do
 return [/plot(0),plot(1),plot(2)/]
end


;******************************************************************
;;主程序
;;**************************************************************
begin

rea=(/"20CR","ERA","JRA55","MERRA","NCEP1","NCEP2","MRE"/)
LeftString=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2","(g)MRE"/)


plot1=new(6,graphic)
plot2=new(6,graphic)
plot3=new(6,graphic)

wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/ 6个数据西位置指数与对应偏相关")
wks2=gsn_open_wks("pdf", "/mnt/f/share/era/data/6个数据东位置指数与对应偏相关")
wks3=gsn_open_wks("pdf", "/mnt/f/share/era/data/6个数据宽度与海温的相关")

m=NewList("lifo")
do i=0,6
  ListPush(m, draw_corr(rea(i),wks1,wks2,wks3,LeftString(i)))
end do
plot_20cr=m[5]
plot_ERA=m[4]
plot_JRA55=m[3]
plot_MERRA=m[2]
plot_NCEP1=m[1]
plot_NCEP2=m[0]

 
 plot1(0)=plot_20cr[0]
 plot1(1)=plot_ERA[0]
 plot1(2)=plot_JRA55[0] 
 plot1(3)=plot_MERRA[0]
 plot1(4)=plot_NCEP1[0]
 plot1(5)=plot_NCEP2[0]

 
 plot2(0)=plot_20cr[1]
 plot2(1)=plot_ERA[1]
 plot2(2)=plot_JRA55[1] 
 plot2(3)=plot_MERRA[1]
 plot2(4)=plot_NCEP1[1]
 plot2(5)=plot_NCEP2[1]
  
  
 plot3(0)=plot_20cr[2]
 plot3(1)=plot_ERA[2]
 plot3(2)=plot_JRA55[2] 
 plot3(3)=plot_MERRA[2]
 plot3(4)=plot_NCEP1[2]
 plot3(5)=plot_NCEP2[2]
  

pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks1, plot1, (/3,2/), pres)
frame(wks1)

gsn_panel(wks2, plot2, (/3,2/), pres)
frame(wks2)

gsn_panel(wks3, plot3, (/3,2/), pres)
frame(wks3)

end
