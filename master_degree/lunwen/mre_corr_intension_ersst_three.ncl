;********************************************************
;計算三个环流的强度变化
;*******************************************************
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
undef("level_mean")
function level_mean(msf_time)   ;;;;select 850hpa-200hpa 
begin
; msf_time=(/msf_time/10^9/)
plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
; dp=(/50/)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end
;;; 
function  ploy(locationxi_rea)
begin
dim=dimsizes(locationxi_rea)
dim1=dim(1)
mnmx      = new ((/3,dim(1)/), typeof(locationxi_rea))
mnmx(0,:) = dim_avg_n_Wrap(locationxi_rea, 0)
mnmx(1,0:dim1-1) = dim_min_n_Wrap(locationxi_rea, 0)
mnmx(2,0:dim1-1) = dim_max_n_Wrap(locationxi_rea, 0)
return mnmx
end 
;;;
begin
    rea= (/"20cr","era","jra55","MERRA","ncep1","ncep2"/)

    id_inten=new((/6,34/),float)
    pa_inten=new((/6,34/),float)
    al_inten=new((/6,34/),float)


    do i=0,dimsizes(rea)-1
    diri="/mnt/f/share/rea_msf/"
    LeftString=rea(i)
    f1=addfile(diri+"msf_time_"+LeftString+"_global.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    
    year1=1979
    year2=2012
    yearend=201212
    msf_time=f1->msf_time(:,:,timeind)
    
    delete(f1)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    delete(timeind)
    printVarSummary(msf_time)
    msf_time=(/msf_time/10^9/)

    msf_mean=msf_time
  
    delete(msf_time)

    ; msf_mean=(/dtrend_msg_n(ispan(1,dimsizes(msf_mean&time),1),msf_mean,False,True,2)/)
      
      

    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")

    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    delete(msf_mean)
    ;;;
    inten_id= dim_avg_n_Wrap(level_mean(msf_djf(:,:,{70:100})),1)
    inten_pa=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{150:210})),1)
    inten_al=dim_avg_n_Wrap(level_mean(msf_djf(:,:,{320:340})),1)
    delete(msf_djf)

    inten_id=(/dim_rmvmean(inten_id)/)
    inten_pa=(/dim_rmvmean(inten_pa)/)
    inten_al=(/dim_rmvmean(inten_al)/)
    ;;;
    
   r_id_pa=escorc_n(inten_pa,inten_al,0,0)
   prob_id_pa=rtest(r_id_pa, dimsizes(inten_id&time),0)  
   print("r_id_pa="+r_id_pa)
   print("prob_id_pa="+prob_id_pa)




    id_inten(i,:)=tofloat(inten_id)
    pa_inten(i,:)=tofloat(inten_pa)
    al_inten(i,:)=tofloat(inten_al)

    ; id_inten(i,:)=(/runave(id_inten(i,:), 3, 0)/)
    ; pa_inten(i,:)=(/runave(pa_inten(i,:), 3, 0)/)
    ; al_inten(i,:)=(/runave(al_inten(i,:), 3, 0)/)
    delete([/inten_id,inten_pa,inten_al/])
  end do

    
    id_inten_avg=dim_avg_n_Wrap(id_inten, 0)
    pa_inten_avg=dim_avg_n_Wrap(pa_inten, 0)
    al_inten_avg=dim_avg_n_Wrap(al_inten, 0)
    
    id_inten_avg!0="time"
    pa_inten_avg!0="time"
    al_inten_avg!0="time"
    id_inten_avg&time=ispan(1, 34, 1)
    pa_inten_avg&time=ispan(1, 34, 1)
    al_inten_avg&time=ispan(1, 34, 1)
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
  ssta=(/rmMonAnnCycTLL(sst(:,:,:))/)
  copy_VarCoords(sst(:,:,:),ssta)

  ; ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
  ; copy_VarCoords(ssta, ssta_de)
  ssta_de=ssta
  ssta_djf=month_to_season(ssta_de, "DJF")

  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
  copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)

  nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))  

;;;
r_id_nino3_4=escorc_n(id_inten, nino3_4, 1, 0)
r_pa_nino3_4=escorc_n(pa_inten, nino3_4, 1, 0)

prob_id_nino3_4=rtest(r_id_nino3_4, dimsizes(id_inten_avg&time), 0)
prob_pa_nino3_4=rtest(r_pa_nino3_4, dimsizes(pa_inten_avg&time), 0)
print("r_id_nino3_4="+r_id_nino3_4)
print("prob_id_nino3_4="+prob_id_nino3_4)
print("r_pa_nino3_4"+r_pa_nino3_4)
print("prob_pa_nino3_4"+prob_pa_nino3_4)


  r_id_sst=escorc_n(id_inten_avg,ssta_s,0,0)  
  r_pa_sst=escorc_n(pa_inten_avg,ssta_s,0,0)
  r_al_sst=escorc_n(al_inten_avg,ssta_s,0,0)
  copy_VarCoords(ssta_s(0,:,:), r_id_sst)
  copy_VarCoords(ssta_s(0,:,:), r_pa_sst)
  copy_VarCoords(ssta_s(0,:,:), r_al_sst)
   


  prob_id_sst=rtest(r_id_sst,dimsizes(id_inten_avg&time),0)
  copy_VarCoords(r_id_sst, prob_id_sst)
  prob_id_sst=where(prob_id_sst.eq.0,999 , prob_id_sst)

  prob_pa_sst=rtest(r_pa_sst,dimsizes(pa_inten_avg&time),0)
  copy_VarCoords(r_pa_sst, prob_pa_sst)
  prob_pa_sst=where(prob_pa_sst.eq.0,999 , prob_pa_sst)
  

  prob_al_sst=rtest(r_al_sst,dimsizes(al_inten_avg&time),0)
  copy_VarCoords(r_al_sst, prob_al_sst)
  prob_al_sst=where(prob_al_sst.eq.0,999 , prob_al_sst)
  printVarSummary(prob_al_sst)
  
;;;


;;;计算线性一致性
  re_id_sst=regCoef_n(id_inten_avg,ssta_s,0,0)  
  re_pa_sst=regCoef_n(pa_inten_avg,ssta_s,0,0)
  re_al_sst=regCoef_n(al_inten_avg,ssta_s,0,0)
  copy_VarCoords(ssta_s(0,:,:), re_id_sst)
  copy_VarCoords(ssta_s(0,:,:), re_pa_sst)
  copy_VarCoords(ssta_s(0,:,:), re_al_sst)
  printVarSummary(re_id_sst)

  rc_id=regCoef_n(ispan(1, 34, 1), id_inten_avg,0,0)
  rc_pa=regCoef_n(ispan(1, 34, 1), pa_inten_avg,0,0)
  rc_al=regCoef_n(ispan(1, 34, 1), al_inten_avg,0,0)
  
  re_id_sst=(/re_id_sst*rc_id/)
  re_pa_sst=(/re_pa_sst*rc_pa/)
  re_al_sst=(/re_al_sst*rc_al/)
  

  rc_sst=regCoef_n(ispan(1, 34, 1), ssta_s,0,0)
  copy_VarCoords(ssta_s(0,:,:), rc_sst)

  rc_id_per=rc_sst
  rc_pa_per=rc_sst
  rc_al_per=rc_sst

  rc_id_per=(/re_id_sst/rc_sst/)
  rc_pa_per=(/re_pa_sst/rc_sst/) 
  rc_al_per=(/re_al_sst/rc_sst/)
  




  rc_id_per=where(rc_sst.eq.rc_sst@_FillValue, rc_id_per@_FillValue, rc_id_per)
  rc_pa_per=where(rc_sst.eq.rc_sst@_FillValue, rc_pa_per@_FillValue, rc_pa_per)
  rc_al_per=where(rc_sst.eq.rc_sst@_FillValue, rc_al_per@_FillValue, rc_al_per)



;;;画图

  wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/集合平均数据强度指数与海温的相关-全球")
  ; gsn_define_colormap(wks1, "BlRe")
    ; gsn_define_colormap(wks, "BlRe ")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  ; res@cnFillPalette="sunshine_9lev"
  res@cnFillPalette="BlueWhiteOrangeRed"
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
  ; res@lbOrientation = "Horizontal"

  res@mpLimitMode           = "LatLon"

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
  plot_con=new(3,graphic)
  plot_zero=new(3,graphic)
  plot(0)= gsn_csm_contour_map(wks1,rc_id_per,res)

  plot(1)= gsn_csm_contour_map(wks1,rc_pa_per ,res)
  plot(2)=gsn_csm_contour_map(wks1, rc_al_per, res)
  

  res@cnFillOn=False
  res@cnLinesOn=True
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   0.02
  res@cnMaxLevelValF       =   0.02
  res@cnLevelSpacingF      =  0.005
  res@gsnLeftString=""
  res@gsnContourNegLineDashPattern=1
  res@cnLineLabelsOn=True

  plot_con(0)= gsn_csm_contour(wks1,re_id_sst,res)
  plot_con(1)= gsn_csm_contour(wks1,re_pa_sst ,res)
  plot_con(2)=gsn_csm_contour(wks1, re_al_sst, res)
  

  overlay(plot(0), plot_con(0))
  overlay(plot(1), plot_con(1))
  overlay(plot(2), plot_con(2))





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

  plot_zero(0)=gsn_csm_contour(wks1,re_id_sst,res_z)
  plot_zero(1) =gsn_csm_contour(wks1,re_pa_sst,res_z)
  plot_zero(2) =gsn_csm_contour(wks1,re_al_sst,res_z)
  overlay(plot(0), plot_zero(0))
  overlay(plot(1), plot_zero(1))
  overlay(plot(2), plot_zero(2))
  




 

; ;;;打点
;   res_pro=True
;   res_pro@gsnAddCyclic=False
;   res_pro@gsnDraw=False
;   res_pro@gsnFrame=False
;   res_pro@cnLinesOn=False
;   res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
;   res_pro@cnMinLevelValF      = 0.00        ; set min contour level
;   res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
;   res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing
;   ; res_pro@cnLevels = ispan(-6,6,1)
;   res_pro@cnLineLabelsOn=False
;   res_pro@tmYRMode = "Automatic"
;   res_pro@cnInfoLabelOn=False

;   iplot=new(3,graphic)

;   iplot(0)=gsn_csm_contour(wks1,prob_id_sst,res_pro)
;   iplot(1)=gsn_csm_contour(wks1,prob_pa_sst,res_pro)
;   iplot(2)=gsn_csm_contour(wks1,prob_al_sst,res_pro)



;   opt                      = True
;   opt@gsnShadeFillType     = "pattern"
;   opt@gsnShadeMid     = 17         ; stipple
;   ; opt@gsnShadeMid   = "green"
;   opt@gsnShadeFillScaleF   = 1.0        ; add extra density
;   opt@gsnShadeFillDotSizeF = 0.0    ; make dots larger
  
;   iplot(0) = gsn_contour_shade(iplot(0), 0,0.05, opt)
;   iplot(1)= gsn_contour_shade(iplot(1), 0,0.05, opt)
;   iplot(2)= gsn_contour_shade(iplot(2), 0,0.05, opt)

; ;;;  
; do i=0,2
; overlay(plot(i), iplot(i))
; end do

pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks1, plot, (/3,1/), pres)

end





















;