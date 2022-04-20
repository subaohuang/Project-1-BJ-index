
function one_test(pv_ano,pv_clm,east_time)
begin
    pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east- pv_clm)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)

 return confi_pv_east
 end

;;NCEP1 U,V 
begin
    fils1="/mnt/f/share/MRE/ud_vd_mre_global.nc"
    
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    level=f1->level
    lat=f1->lat
    lon=f1->lon
    time=f1->time(timeind)
    uwnd_necp1_r=f1->ud_new({850},{30:-30},{90:300},timeind)   ;;;插值前先逆着取
    vwnd_necp1_r=f1->vd_new({850},{30:-30},{90:300},timeind)   ;;;
    vwnd_necp1_r=(/vwnd_necp1_r*-1/)
    delete(uwnd_necp1_r&time)
    delete(vwnd_necp1_r&time)
    uwnd_necp1_r&time=date(timeind)
    vwnd_necp1_r&time=date(timeind)
    
    uwnd_necp1=uwnd_necp1_r(time|:,lat|:,lon|:)
    vwnd_necp1=vwnd_necp1_r(time|:,lat|:,lon|:)
   
    delete(uwnd_necp1_r)
    delete(vwnd_necp1_r)

    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
;;;;

;;;readin potensial velocity
    f3=addfile("/mnt/f/share/MRE/velocity_potensial_mre.nc", "r")
    date=cd_calendar(f3->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
  
    pv=f3->pv(timeind,{850},{30:-30},:)
    pv=(/pv/10^6/)
    delete(pv&time)
    pv&time=date(timeind)
    delete(date)
    delete(timeind)


;;;
 ;;读入ERSST海温数据
    m=""+ispan(1979,2013,1)+""
    diri="/mnt/f/share/ersst/"
    fils=new((/dimsizes(m),12/),string)
    do i=0,dimsizes(m)-1
    fils2=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
    fils(i,:)= fils2
    end do
    fils_new=onedtond(fils, dimsizes(m)*12)
    fileh=new(dimsizes(m)*12,string)
    fileh(0:dimsizes(m)*12-1)=fils_new

    f2=addfiles(fileh, "join")
    printVarSummary(f2)
    date=cd_calendar(f2[:]->time, 1)
    timeind=ind(date.le.201312.and.date.ge.197901)
    sst=f2[:]->sst(timeind,0,{20:-20},:)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)

    ssta_djf=month_to_season(ssta_de, "DJF")
    ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    copy_VarCoords(ssta_djf(0:dimsizes(m)-2,:,:), ssta_s)



    emi=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))-0.5*dim_avg_n(ssta_s(:,{20:-10},{125:145}), (/1,2/))-0.5*dim_avg_n(ssta_s(:,{5:-15},{250:290}),(/1,2/))
    
    emi_center=dim_avg_n(ssta_s(:,{10:-10},{165:220}),(/1,2/))

    emi_std=dim_standardize_n(emi, 1, 0)
    std_emi=dim_stddev_n(emi, 0)

    nino3=dim_avg_n(ssta_s(:,{5:-5},{210:270}), (/1,2/))
    nino3_std=dim_standardize_n(nino3, 1, 0)
    std_nino3=dim_stddev_n(nino3, 0)
     
    nino3_4=dim_avg_n(ssta_s(:,{5:-5},{190:240}), (/1,2/))    
;;;

;;去趋势
   u_de=dtrend_msg_n(ispan(1,dimsizes(uwnd_necp1&time),1),uwnd_necp1,False,True,0)
   copy_VarCoords(uwnd_necp1, u_de)
   v_de=dtrend_msg_n(ispan(1,dimsizes(vwnd_necp1&time),1),vwnd_necp1,False,True,0)
   copy_VarCoords(vwnd_necp1, v_de)
   pv_de=dtrend_msg_n(ispan(1,dimsizes(pv&time),1),pv,False,True,0)
   copy_VarCoords(pv, pv_de)
;;;   
  printVarSummary(u_de)

;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
 ; u_son=month_to_season(u_de(time|:,lat|:,lon|:), "SON")
 u_djf=month_to_season(u_de(time|:,lat|:,lon|:), "DJF")

 ; v_son=month_to_season(v_de(time|:,lat|:,lon|:), "SON")
 v_djf=month_to_season(v_de(time|:,lat|:,lon|:), "DJF")

 ; pv_son=month_to_season(pv_de(time|:,lat|:,lon|:), "SON")
 pv_djf=month_to_season(pv_de(time|:,lat|:,lon|:), "DJF")

 ; u_son&time=ispan(1958,2012,1)
 u_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 ; v_son&time=ispan(1958,2012,1)
 v_djf&time=ispan(1979,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 
 ;;注意点
    u_djf(0:dimsizes(u_djf&time)-2,:,:)=(/u_djf(1:,:,:)/)
    u_djf(dimsizes(u_djf&time)-1,:,:)=(/u_de({201212},:,:)/)

    v_djf(0:dimsizes(v_djf&time)-2,:,:)=(/v_djf(1:,:,:)/)
    v_djf(dimsizes(v_djf&time)-1,:,:)=(/v_de({201212},:,:)/)
    
    pv_djf(0:dimsizes(pv_djf&time)-2,:,:)=(/pv_djf(1:,:,:)/)
    pv_djf(dimsizes(pv_djf&time)-1,:,:)=(/pv_de({201212},:,:)/)

    u_cold=u_djf
    u_cold=(/(u_djf(:,:,:))*1.0/)

    v_cold=v_djf
    v_cold=(/(v_djf(:,:,:))*1.0/)
    
    pv_cold=pv_djf
    pv_cold=(/(pv_djf(:,:,:))*1.0/)

    printVarSummary(u_cold)

 ;  msf_clm=dim_avg_n((/msf_son(:,:,:)+msf_djf(:,:,:)/)*0.5,0)
 ; copy_VarCoords(msf_son(1,:,:), msf_clm)


;;;挑选年份 modoki year
 modokiind=ind(emi_std.ge.std_emi)
 u_modoki=dim_avg_n(u_cold(modokiind,:,:),0)
 copy_VarCoords(u_djf(1,:,:), u_modoki)
 v_modoki=dim_avg_n(v_cold(modokiind,:,:),0)
 copy_VarCoords(v_djf(1,:,:), v_modoki)
;;;
printVarSummary(u_cold)
;;;挑选年份 canino  enino year
 caind=ind(nino3_std.ge.std_nino3)
 u_canonical=dim_avg_n(u_cold(caind,:,:),0)
 copy_VarCoords(u_djf(1,:,:), u_canonical)
 v_canonical=dim_avg_n(v_cold(caind,:,:),0)
 copy_VarCoords(v_djf(1,:,:), v_canonical)
;
;;;挑选年份 canino  lanina year
 laind=ind(nino3_std.le.(/-1*std_nino3/))
 u_lanina=dim_avg_n(u_cold(laind,:,:),0)
 copy_VarCoords(u_djf(1,:,:), u_lanina)
  v_lanina=dim_avg_n(v_cold(laind,:,:),0)
 copy_VarCoords(v_djf(1,:,:), v_lanina)
;;; normal time
 noind=ind(nino3_std.ge.(/-0.5*std_nino3/).and.nino3_std.le.(/0.5*std_nino3/))
 u_noind=dim_avg_n(u_cold(noind,:,:),0)
 copy_VarCoords(u_djf(1,:,:), u_noind)
 v_noind=dim_avg_n(v_cold(noind,:,:),0)
 copy_VarCoords(v_djf(1,:,:), v_noind)
;;;

;;;
;;;挑选年份 modoki year
 pv_modoki=dim_avg_n_Wrap(pv_cold(modokiind,:,:),0)
;;;
;;;挑选年份 canino  enino year
 pv_canonical=dim_avg_n_Wrap(pv_cold(caind,:,:),0)
;
;;;挑选年份 canino  lanina year
 pv_lanina=dim_avg_n_Wrap(pv_cold(laind,:,:),0)
;;; normal time
 pv_noind=dim_avg_n_Wrap(pv_cold(noind,:,:),0)
;;;


pv_clm=dim_avg_n_Wrap(pv_cold, 0)
u_clm=dim_avg_n_Wrap(u_cold, 0)
v_clm=dim_avg_n_Wrap(v_cold, 0)

;;;显著性检验

pv_cano_config=one_test(pv_cold,pv_clm,caind)
u_cano_config=one_test(u_cold,u_clm,caind)
v_cano_config=one_test(v_cold,v_clm,caind)

pv_modo_config=one_test(pv_cold,pv_clm,modokiind)
u_modo_config=one_test(u_cold,u_clm,modokiind)
v_modo_config=one_test(v_cold,v_clm,modokiind)


pv_lani_config=one_test(pv_cold,pv_clm,laind)
u_lani_config=one_test(u_cold,u_clm,laind)
v_lani_config=one_test(v_cold,v_clm,laind)




; u_canonical=where(u_cano_config.le.0.05,u_canonical , u_canonical@_FillValue)
; u_modoki  =where(u_modo_config.le.0.05,u_modoki , u_modoki@_FillValue)
; u_lanina=where(u_lanina.le.0.05,u_lanina , u_lanina@_FillValue)



; v_canonical=where(v_cano_config.le.0.05,v_canonical , v_canonical@_FillValue)
; v_modoki  =where(v_modo_config.le.0.05,v_modoki , v_modoki@_FillValue)
; v_lanina=where(v_lanina.le.0.05,v_lanina , v_lanina@_FillValue)


uw_canonical_t=where(u_cano_config.gt.v_cano_config,u_cano_config,v_cano_config)
uw_modoki_t  =where(u_modo_config.gt.v_cano_config,u_cano_config, v_cano_config)
uw_lanina_t=where(u_lanina.gt.v_cano_config,u_cano_config, v_cano_config)






 wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/1979-2012各ENSO事件的辐散速度势图")
 gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  
  res                  = True                      ; plot mods desired
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  ; res@gsnScalarContour        = True
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False
  res@lbLabelBarOn=False

  res@gsnLeftString="canonical enino 850hpa wind"
  res@mpLimitMode      = "LatLon"
  res@mpMinLatF        = min(u_modoki&lat)
  res@mpMaxLatF        = max(u_modoki&lat)
  res@mpMinLonF        = min(u_modoki&lon)
  res@mpMaxLonF        = max(u_modoki&lon)
  res@mpCenterLonF       =avg(u_modoki&lon)
  res@mpFillDrawOrder="PostDraw"
  res@cnLineLabelsOn=False
 
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =  -7
  res@cnMaxLevelValF       =  7
  res@cnLevelSpacingF      =  1
  res@mpShapeMode ="FreeAspect"
  ; res@vpHeightF=0.38
  ; res@vpWidthF=0.62

  ; plot=new(4,graphic)
  ; plot(0)= gsn_csm_vector_scalar_map(wks,u_canonical,v_canonical,pv_canonical,res)
  ; res@gsnLeftString="canonical lanina 850hpa wind"
  ; plot(1)= gsn_csm_vector_scalar_map(wks,u_lanina,v_lanina,pv_lanina, res)
  ; res@gsnLeftString="enino modoki 850hpa wind"
  ; plot(2)= gsn_csm_vector_scalar_map(wks,u_modoki,v_modoki,pv_modoki ,res)
  ; res@gsnLeftString="normal 850hpa wind"
  ; plot(3)= gsn_csm_vector_scalar_map(wks,u_noind,v_noind,pv_noind, res)


  plot=new(4,graphic)
  plot(0)= gsn_csm_contour_map(wks,pv_canonical,res)
  res@gsnLeftString="canonical lanina 850hpa wind"
  plot(1)= gsn_csm_contour_map(wks,pv_lanina, res)
  res@gsnLeftString="enino modoki 850hpa wind"
  plot(2)= gsn_csm_contour_map(wks,pv_modoki ,res)
  res@gsnLeftString="normal 850hpa wind"
  plot(3)= gsn_csm_contour_map(wks,pv_noind, res)
  
  ;;;打点
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@cnLinesOn=False
        res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
        res_pro@cnMinLevelValF      = 0.00        ; set min contour level
        res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
        res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

        res_pro@cnLineLabelsOn=False
        res_pro@tmYRMode = "Automatic"
        res_pro@cnInfoLabelOn=False

        iplot=new(3,graphic)

        iplot(0)=gsn_csm_contour(wks,pv_cano_config,res_pro)
        iplot(1)=gsn_csm_contour(wks,pv_lani_config,res_pro)
        iplot(2)=gsn_csm_contour(wks,pv_modo_config ,res_pro)
     

        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 0.2         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
        
        do i=0,2
        iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
      end do
        ;;;  
        do i=0,2
        overlay(plot(i), iplot(i))
      end do
        
  ;resWind
resWind=True
resWind@gsnFrame=False
resWind@gsnDraw=False
resWind@gsnAddCyclic=False
resWind@gsnSpreadColors      = True
resWind@lbLabelBarOn=False

resWind@vcRefMagnitudeF          = 2.5           ; define vector ref mag
resWind@vcRefLengthF             = 0.045          ; define length of vec ref
resWind@vcGlyphStyle             = "CurlyVector"      ; turn on curly vectors
resWind@vcMinDistanceF           = 0.015          ; thin vectors
resWind@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角


resWind@vcMonoLineArrowColor    = False
resWind@vcLevelSelectionMode= "ExplicitLevels"
resWind@vcLevels= (/0.05/)
resWind@vcLevelColors= (/ "black","gray"/)


plot_wind_ttest=new(4,graphic)
plot_wind_ttest(0)=gsn_csm_vector_scalar(wks,u_canonical,v_canonical ,uw_canonical_t,resWind)
plot_wind_ttest(1)=gsn_csm_vector_scalar(wks,u_lanina ,v_lanina ,uw_lanina_t ,resWind)
plot_wind_ttest(2)=gsn_csm_vector_scalar(wks,u_modoki,v_modoki ,uw_modoki_t ,resWind)
resWind@vcMonoLineArrowColor    = True
plot_wind_ttest(3)=gsn_csm_vector(wks,u_noind,v_noind,resWind)







do i=0,3
overlay(plot(i),plot_wind_ttest(i))
end do
pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks, plot, (/2,2/), pres)


end
