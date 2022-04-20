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
begin
rea=(/"20CR","ERA","JRA55","MERRA","NCEP1","NCEP2"/)
LeftString=(/"(a)20CR","(b)ERAIM","(c)JRA55","(d)MERRA","(e)NCEP1","(f)NCEP2"/)

    diri="/mnt/f/share/rea_msf/"
    locationxi_r=new((/dimsizes(rea),34/),float)
    locationdong_r=new((/dimsizes(rea),34/),float)
    width_r=new((/dimsizes(rea),34/),float)
 do i=0,dimsizes(rea)-1
    locationxi_r(i,:)=asciiread(diri+rea(i)+"_msf_djf_locationxi_node_index.ascii", 34, "float")
    locationdong_r(i,:)=asciiread(diri+rea(i)+"_msf_djf_locationdong_node_index.ascii", 34, "float")   
    width_r(i,:)=locationdong_r(i,:)-locationxi_r(i,:)
end do
locationxi_r_mean=dim_avg_n(locationxi_r, 0)
locationdong_r_mean=dim_avg_n(locationdong_r, 0)
width_r_mean=dim_avg_n(width_r, 0)


ipindex1  =dim_rmvmean(locationxi_r_mean)
ipindex2 =dim_rmvmean(locationdong_r_mean)

ipindex3=width_r_mean
ipindex3=dim_rmvmean(ipindex3)

ipindex1=(/dtrend_msg_n(ispan(1,34,1),ipindex1,False,True,0)/)
ipindex2=(/dtrend_msg_n(ispan(1,34,1),ipindex2,False,True,0)/)
ipindex3=(/dtrend_msg_n(ispan(1,34,1),ipindex3,False,True,0)/)
  


;;;读入海温数据
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
;;;处理成异常的冬季平均场
  ; ssta=(/rmMonAnnCycTLL(sst(:,:,:))/)
  ; copy_VarCoords(sst(:,:,:), ssta)
  ssta=(/rmMonAnnCycTLL(sst(:,:,{0:}))/)
  copy_VarCoords(sst(:,:,{0:}),ssta)

  ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
  copy_VarCoords(ssta, ssta_de)
  ssta_djf=month_to_season(ssta_de, "DJF")
  ssta_son=month_to_season(ssta_de, "SON")


  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
  ; ssta_s=(/ssta_son(0:33,:,:)/)
  copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
;;;


 emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
 nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
;;;
  r_xi_dong=escorc_n(ipindex1,ipindex2,0,0)

  print("r_xi_dong="+r_xi_dong)
; 
   
    ; ;remove enso cycle
    ;  nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))  
    ;  a1=regCoef_n(nino3_4, ipindex1, 0, 0)
    ;  ipindex1_enso=a1*nino3_4+a1@yintercept
    ;  ipindex1=(/ipindex1- ipindex1_enso/)

    ;  a2=regCoef_n(nino3_4, ipindex2, 0, 0)
    ;  ipindex2_enso=a2*nino3_4+a2@yintercept
    ;  ipindex2=(/ipindex2- ipindex2_enso/)
    ;  a3=regCoef_n(nino3_4, ipindex3, 0, 0)
    ;  ipindex3_enso=a3*nino3_4+a3@yintercept
    ;  ipindex3=(/ipindex3- ipindex3_enso/)
   
    ; printVarSummary(nino3_4)
    ; printVarSummary(ssta_s)
    ;  a4=regCoef_n(nino3_4, ssta_s, 0, 0)
    ;  a4_comform=conform_dims(dimsizes(ssta_s), a4, (/1,2/))
    ; printVarSummary(a4)
    ;  l=reshape(a4@yintercept, dimsizes(ssta_s(0,:,:)))
    ;  l_comform=conform_dims(dimsizes(ssta_s), l, (/1,2/))
    ;  ssta_s_enso=ssta_s
    ;  do i=0,dimsizes(ssta_s&lat)-1
    ;   do j=0,dimsizes(ssta_s&lon)-1
    ;   ssta_s_enso(:,i,j)=a4(i,j)*nino3_4+l(i,j)
    ;  end do 
    ; end do
    ;  ssta_s=(/ssta_s-ssta_s_enso/)
    ; ;;;
   

  r_xi_sst=escorc_n(ipindex1,ssta_s,0,0)  
  r_dong_sst=escorc_n(ipindex2,ssta_s,0,0)
  copy_VarCoords(ssta_s(0,:,:), r_xi_sst)
  copy_VarCoords(ssta_s(0,:,:), r_dong_sst)
;;;



  r_xi_sst_pa=(/(r_xi_sst- r_xi_dong*r_dong_sst)/sqrt((1- r_xi_dong^2)*(1- r_dong_sst^2))/)
  r_dong_sst_pa=(/(r_dong_sst- r_xi_dong*r_xi_sst)/sqrt((1- r_xi_dong^2)*(1- r_xi_sst^2))/)
  r_cha=(/r_xi_sst_pa- r_dong_sst_pa/)

  copy_VarCoords(r_xi_sst, r_xi_sst_pa)
  copy_VarCoords(r_xi_sst, r_dong_sst_pa)
  copy_VarCoords(r_xi_sst, r_cha)
  
  r_emi_dong=escorc_n(ipindex2, emi, 0, 0)
  print("r_emi_dong="+r_emi_dong)

  r_emi_xi=escorc_n(ipindex1, emi, 0, 0)
  print("r_emi_xi="+r_emi_xi)
   
  r_nino3_dong=escorc_n(ipindex2, nino3, 0, 0)
  print("r_nino3_dong="+r_nino3_dong)

  r_nino3_xi=escorc_n(ipindex1, nino3, 0, 0)
  print("r_nino3_xi="+r_nino3_xi)
  
  
;;
 
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

   prob_r_xi_sst=rtest(r_xi_sst,dimsizes(ssta_s&time),0)
  copy_VarCoords(r_xi_sst, prob_r_xi_sst)
  prob_r_xi_sst=where(prob_r_xi_sst.eq.0,999 , prob_r_xi_sst)
  
  prob_r_dong_sst=rtest(r_dong_sst,dimsizes(ssta_s&time),0)
  copy_VarCoords(r_dong_sst, prob_r_dong_sst)
  prob_r_dong_sst=where(prob_r_dong_sst.eq.0,999 , prob_r_dong_sst)

  

;;;画图

  wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/集合平均数据位置指数与海温的相关-全球2")
  gsn_define_colormap(wks1, "BlueWhiteOrangeRed")
    ; gsn_define_colormap(wks, "BlRe ")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnLeftString="(a)Western Edge"
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
  ; res@mpMinLatF             = min(ssta&lat)
  ; res@mpMaxLatF             = max(ssta&lat)
  ; res@mpMinLonF             = min(ssta&lon)
  ; res@mpMaxLonF             = max(ssta&lon)
  ; res@mpCenterLonF          =avg(ssta&lon)
    res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  res@mpCenterLonF          =180
  res@mpShapeMode   ="FreeAspect"
  res@vpHeightF=0.35
  res@vpWidthF=0.65

  res@tmXTOn=False
  res@tmYROn=False
  

  plot=new(3,graphic)
  plot_zero=new(3,graphic)
  plot(0)= gsn_csm_contour_map(wks1,r_xi_sst,res)
  ; res@cnMinLevelValF       =   -0.6
  ; res@cnMaxLevelValF       =   0.6
  ; res@cnLevelSpacingF      =  0.1
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
    res@gsnLeftString="(b)Eastern Edge"
  plot(1)= gsn_csm_contour_map(wks1,r_dong_sst,res)

  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  res@gsnLeftString="(c)Width"
  plot(2)=gsn_csm_contour_map(wks1, r_width_sst, res)

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

  plot_zero(0)=gsn_csm_contour(wks1,r_xi_sst,res_z)
  plot_zero(1) =gsn_csm_contour(wks1,r_dong_sst,res_z)
  plot_zero(2) =gsn_csm_contour(wks1,r_width_sst,res_z)
  overlay(plot(0), plot_zero(0))
  overlay(plot(1), plot_zero(1))
  overlay(plot(2), plot_zero(2))



;;;打点
  res_pro=True
  res_pro@gsnAddCyclic=False
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False

  iplot=new(3,graphic)

  iplot(0)=gsn_csm_contour(wks1,prob_r_xi_sst,res_pro)
  iplot(1)=gsn_csm_contour(wks1,prob_r_dong_sst,res_pro)
  iplot(2)=gsn_csm_contour(wks1,prob_width_sst,res_pro)

 printVarSummary(prob_r_xi_sst)
 printVarSummary(prob_width_sst)

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeMid   = "green"
  opt@gsnShadeFillScaleF   = 1.0        ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0    ; make dots larger
  
  iplot(0) = gsn_contour_shade(iplot(0), 0,0.05, opt)
  iplot(1)= gsn_contour_shade(iplot(1), 0,0.05, opt)
  iplot(2)= gsn_contour_shade(iplot(2), 0,0.05, opt)

;;;  
do i=0,2
overlay(plot(i), iplot(i))
end do





pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks1, plot, (/3,1/), pres)

end

















