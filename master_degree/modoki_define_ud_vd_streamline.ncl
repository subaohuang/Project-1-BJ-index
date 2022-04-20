
;;NCEP1 U,V 
begin
    fils1="/mnt/f/share/MRE/ud_vd_mre_global.nc"
    
    year1=1958
    year2=2012
    dims=year2- year1+1
    
    yearstart=195801
    yearend=201212

    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.yearstart.and.date.le.201212)
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
    timeind=ind(date.ge.yearstart.and.date.le.201312)
    
    pv=f3->pv(timeind,{850},{30:-30},:)
    delete(pv&time)
    pv&time=date(timeind)
    delete(date)
    delete(timeind)


;;;

    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.yearstart.and.date.le.201312)
    sst=lonFlip(f2->sst(timeind,{20:-20},:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta_de=dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)
    copy_VarCoords(ssta, ssta_de)
    ssta_son=month_to_season(ssta_de, "SON")
    ssta_djf=month_to_season(ssta_de, "DJF")
 
    
    ssta_s=(/ssta_son(0:dims-1,:,:)+ssta_djf(1:,:,:)/)*0.5
    
    copy_VarCoords(ssta_son(0:dims-1,:,:), ssta_s)
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
 u_son=month_to_season(u_de(time|:,lat|:,lon|:), "SON")
 u_djf=month_to_season(u_de(time|:,lat|:,lon|:), "DJF")

 v_son=month_to_season(v_de(time|:,lat|:,lon|:), "SON")
 v_djf=month_to_season(v_de(time|:,lat|:,lon|:), "DJF")

 pv_son=month_to_season(pv_de(time|:,lat|:,lon|:), "SON")
 pv_djf=month_to_season(pv_de(time|:,lat|:,lon|:), "DJF")

 u_son&time=ispan(year1,2012,1)
 u_djf&time=ispan(year1,2012,1)  ;;;因为这里第一个是year1 1月和2月的平均，第二个是year112,195901,195902的平均
 v_son&time=ispan(year1,2012,1)
 v_djf&time=ispan(year1,2012,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
 
 ;;注意点
    u_djf(0:dimsizes(u_djf&time)-2,:,:)=(/u_djf(1:,:,:)/)
    u_djf(dimsizes(u_djf&time)-1,:,:)=(/u_de({201212},:,:)/)

    v_djf(0:dimsizes(v_djf&time)-2,:,:)=(/v_djf(1:,:,:)/)
    v_djf(dimsizes(v_djf&time)-1,:,:)=(/v_de({201212},:,:)/)
    
    pv_djf(0:dimsizes(pv_djf&time)-2,:,:)=(/pv_djf(1:,:,:)/)
    pv_djf(dimsizes(pv_djf&time)-1,:,:)=(/pv_de({201212},:,:)/)

    u_cold=u_son
    u_cold=(/(u_djf(:,:,:))*1.0/)

    v_cold=v_son
    v_cold=(/(v_djf(:,:,:))*1.0/)
    
    pv_cold=pv_son
    pv_cold=(/(pv_djf(:,:,:))*1.0/)

    printVarSummary(u_cold)

 ;  msf_clm=dim_avg_n((/msf_son(:,:,:)+msf_djf(:,:,:)/)*0.5,0)
 ; copy_VarCoords(msf_son(1,:,:), msf_clm)


;;;挑选年份 modoki year
 modokiind=ind(emi_std.ge.std_emi)
 u_modoki=dim_avg_n(u_cold(modokiind,:,:),0)
 copy_VarCoords(u_son(1,:,:), u_modoki)
 v_modoki=dim_avg_n(v_cold(modokiind,:,:),0)
 copy_VarCoords(v_son(1,:,:), v_modoki)
;;;
printVarSummary(u_cold)
;;;挑选年份 canino  enino year
 caind=ind(nino3_std.ge.std_nino3)
 u_canonical=dim_avg_n(u_cold(caind,:,:),0)
 copy_VarCoords(u_son(1,:,:), u_canonical)
 v_canonical=dim_avg_n(v_cold(caind,:,:),0)
 copy_VarCoords(v_son(1,:,:), v_canonical)
;

;;;挑选年份 canino  lanina year
 laind=ind(nino3_std.le.(/-1*std_nino3/))
 u_lanina=dim_avg_n(u_cold(laind,:,:),0)
 copy_VarCoords(u_son(1,:,:), u_lanina)
  v_lanina=dim_avg_n(v_cold(laind,:,:),0)
 copy_VarCoords(v_son(1,:,:), v_lanina)
;;; normal time
 noind=ind(nino3_std.ge.(/-0.5*std_nino3/).and.nino3_std.le.(/0.5*std_nino3/))
 u_noind=dim_avg_n(u_cold(noind,:,:),0)
 copy_VarCoords(u_son(1,:,:), u_noind)
 v_noind=dim_avg_n(v_cold(noind,:,:),0)
 copy_VarCoords(v_son(1,:,:), v_noind)
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

 wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/各ENSO事件的辐散速度势图")
 gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  
  res                  = True                      ; plot mods desired
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@gsnScalarContour        = True
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
  res@cnLineLabelsOn=False
 

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =  -6*10^6
  res@cnMaxLevelValF       =  6*10^6
  res@cnLevelSpacingF      =  10^6

  ; res@vcRefMagnitudeF          = 3.0           ; define vector ref mag
  ; res@vcRefLengthF             = 0.04          ; define length of vec ref
  ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcMinDistanceF           = 0.012           ; thin vectors
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角


  res@vcRefMagnitudeF          = 2.5           ; define vector ref mag
  res@vcRefLengthF             = 0.045          ; define length of vec ref
  res@vcGlyphStyle             = "CurlyVector"      ; turn on curly vectors
  res@vcMinDistanceF           = 0.012           ; thin vectors
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角



  plot=new(4,graphic)
  ;;;画流线
  ; plot(0)= gsn_csm_vector_map(wks,u_canonical,v_canonical,res)
  ; res@gsnLeftString="canonical lanina 850hpa wind"
  ; plot(1)= gsn_csm_vector_map(wks,u_lanina,v_lanina,res)
  ; res@gsnLeftString="enino modoki 850hpa wind"
  ; plot(2)= gsn_csm_vector_map(wks,u_modoki,v_modoki,res)
  ; res@gsnLeftString="normal 850hpa wind"
  ; plot(3)= gsn_csm_vector_map(wks,u_noind,v_noind,res)
 
  plot(0)= gsn_csm_vector_scalar_map(wks,u_canonical,v_canonical,pv_canonical,res)
  res@gsnLeftString="canonical lanina 850hpa wind"
  plot(1)= gsn_csm_vector_scalar_map(wks,u_lanina,v_lanina,pv_lanina, res)
  res@gsnLeftString="enino modoki 850hpa wind"
  plot(2)= gsn_csm_vector_scalar_map(wks,u_modoki,v_modoki,pv_modoki ,res)
  res@gsnLeftString="normal 850hpa wind"
  plot(3)= gsn_csm_vector_scalar_map(wks,u_noind,v_noind,pv_noind, res)


pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks, plot, (/2,2/), pres)


end
