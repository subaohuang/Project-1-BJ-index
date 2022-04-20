;;********************************************************
;;计算ERA的u,v,slp的场到P,  Q-term以及其动力项，还有观测的U,V,SLP的趋势
;;
;;lnshuheng 2020.1.4
;**********************************************************
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


function linear_congruent(diri,variable,startmon,endmon,p_t)
begin
;==========================================================
;;;回归
;;==========================================================
   ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
    ;;;
    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;
    f1 = addfile(diri,"r")
    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    if variable.eq."sst" then 
      sst_r=lonFlip(f1->$variable$(timeind,:,:))
      sst = sst_r(:,{lat1:lat2},{lon1:lon2})
      delete(sst_r)
    else if variable.eq."msl" then 
      sst=short2flt(f1->$variable$(timeind,{lat1:lat2},{lon1:lon2}))
    else if variable.eq."u".or.variable.eq."v" then 
      sst=short2flt(f1->$variable$(timeind,{850},{lat1:lat2},{lon1:lon2}))
    else 
      sst=f1->$variable$(timeind,{lat1:lat2},{lon1:lon2})
    end if 
    end if
    end if


   
   ;;;处理成年异常数据
     sst = (/rmMonAnnCycTLL(sst)/)
     sst_annual = month_to_annual(sst, 1)
     delete(sst)
     p_t = (/rmMonAnnCycTLL(p_t)/)
     p = dim_avg_n_Wrap(p_t(:,{0:15},{125:160}), (/1,2/))
     
     p_annual = month_to_annual(p, 1)
     p_annual = (/dim_standardize(p_annual, 1)/)
   

     ;;;计算趋势,计算线性一致性之前先把趋势计算了
     
     sst_rc  = regCoef_n(ispan(1, yearnum, 1), sst_annual, 0, 0)
     ; sst_rc = sst_rc*10

     sst_rc = sst_rc*yearnum

     copy_VarCoords(sst_annual(1,:,:), sst_rc)

     p_rc  = regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
     p_rc = p_rc*yearnum
  
     sst_annual = (/dtrend_msg_n(ispan(1,dim(0),1),sst_annual,False,True,0)/)
     p_annual = (/dtrend_msg_n(ispan(1,dim(0),1),p_annual,False,True,0)/)

   ;;回归
    ;;1: 回归到P
    sst_p_rc = regCoef_n(p_annual, sst_annual, 0, 0)
    copy_VarCoords(sst_annual(1,:,:), sst_p_rc)
   ;;

    ;;;显著性检验
    prob_sst = rc_test(sst_rc)

    prob_sst_p = rc_test(sst_p_rc)
    ;;;;

  ;;;趋势一致性
    sst_p_rc = (/sst_p_rc*p_rc/)
  ;;;;
  
  if variable.eq."m_zonal".or.variable.eq."m_meri" then
    sst_rc = (/sst_rc*10^4/)
    sst_p_rc = (/sst_p_rc*10^4/)
  end if
  
  return  [/sst_rc , sst_p_rc , prob_sst , prob_sst_p/]
end 


begin
  ;;;time span
    startmon =197901
    endmon = 201612
  ;;;

  ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
   ;;;

;;;读入GPCP降水场
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p_t=f2->precip(timeind,:,:)
;;;

;;;SST 
   diri1 = "/home/ys17-19/lsh/data/sst/HadISST_sst.nc"
   diri2 =  "/home/yangsong3/data-observation/ERA-interim/monthly/uwnd.mon.mean.nc"
   diri3 =  "/home/yangsong3/data-observation/ERA-interim/monthly/vwnd.mon.mean.nc"
   diri4 =  "/home/yangsong3/data-observation/ERA-interim/monthly/slp.mon.mean.nc"
   diri5 = "/home/ys17-19/lsh/data/era/era_m_zonal_meri_ERAI_500hpa.nc"
;;;
   sst_info = linear_congruent(diri1,"sst",startmon,endmon,p_t) 
   uwnd_info = linear_congruent(diri2,"u",startmon,endmon,p_t)   
   vwnd_info = linear_congruent(diri3,"v",startmon,endmon,p_t)
   slp_info = linear_congruent(diri4,"msl",startmon,endmon,p_t)
   m_zonal_info = linear_congruent(diri5,"m_zonal",startmon,endmon,p_t)
   m_meri_info = linear_congruent(diri5,"m_meri",startmon,endmon,p_t)
;;;;

;;;;显著性检验
    
    slp_rc = slp_info[0]
    prob_slp = slp_info[2]

    uwnd_rc = uwnd_info[0]
    prob_uwnd = uwnd_info[2]
    
    vwnd_rc = vwnd_info[0]
    prob_vwnd = vwnd_info[2]

    slp_p_rc = slp_info[1]
    prob_slp_p = slp_info[3]

    uwnd_p_rc = uwnd_info[1]
    prob_uwnd_p = uwnd_info[3]

    vwnd_p_rc = vwnd_info[1]
    prob_vwnd_p = vwnd_info[3]

  
    ; slp_rc = where(prob_slp.le.0.1, slp_rc, slp_rc@_FillValue)
    ; slp_p_rc = where(prob_slp_p.le.0.1, slp_p_rc, slp_p_rc@_FillValue)
 

    uwnd_rc = where(prob_uwnd.le.0.1.or.prob_vwnd.le.0.1, uwnd_rc, uwnd_rc@_FillValue)
    uwnd_p_rc = where(prob_uwnd_p.le.0.1.or.prob_vwnd_p.le.0.1, uwnd_p_rc, uwnd_p_rc@_FillValue)

    vwnd_rc = where(prob_vwnd.le.0.1.or.prob_vwnd.le.0.1, vwnd_rc, vwnd_rc@_FillValue)
    vwnd_p_rc = where(prob_vwnd_p.le.0.1.or.prob_vwnd_p.le.0.1, vwnd_p_rc, vwnd_p_rc@_FillValue)
;;;;;



;;;;draw 

;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "SST,U,V,SLP,Local&Meri WC congruent_on_P_2016")
  plot=new(8,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level

  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.4

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2 

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
 

 

  ;;;;;SST 
    res@gsnLeftString = "(a)observed SST trend" 
    plot(0) = gsn_csm_contour_map(wks,sst_info[0],res)

    res@gsnLeftString = "(b)linearly congruent with Pr" 
    res@cnMinLevelValF       =   -0.5
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.1

    plot(1) = gsn_csm_contour_map(wks,sst_info[1],res)
  ;;;;;;



  ;;;;WC

    res@cnMinLevelValF       =   -6.0
    res@cnMaxLevelValF       =   6.0
    res@cnLevelSpacingF      =  1.0

    res@gsnLeftString = "(e)observed lcoal WC linear trend" 

    plot(4) = gsn_csm_contour_map(wks,m_zonal_info[0] ,res)

    res@gsnLeftString = "(f)linearly congruent with Pr"

    plot(5) = gsn_csm_contour_map(wks,m_zonal_info[1] ,res)
    
    res@gsnLeftString = "(g)observed lcoal HC linear trend" 

    plot(6) = gsn_csm_contour_map(wks,m_meri_info[0] ,res)

    res@gsnLeftString = "(h)linearly congruent with Pr"

    plot(7) = gsn_csm_contour_map(wks,m_meri_info[1] ,res)
  ;;;;;


  ;;;;;U,V.SLP

    res@cnMinLevelValF       =   -160
    res@cnMaxLevelValF       =   160
    res@cnLevelSpacingF      =  40


    res@gsnScalarContour        = True
    ;;;;;;;LineArrow
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.04
    res@vcRefMagnitudeF          = 4.0
    ; res@vcRefMagnitudeF          = 0.1
    res@vcRefLengthF             = 0.08     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.008
    res@vcLineArrowHeadMinSizeF=0.008
    res@vcVectorDrawOrder="PostDraw"
    res@mpFillOn =False

    res@gsnLeftString = "(c)observed U,V,SLP linear trend" 
    plot(2) = gsn_csm_vector_scalar_map(wks,uwnd_rc,vwnd_rc,slp_rc ,res)
    res@vcRefMagnitudeF          = 2.0
    res@cnMinLevelValF       =   -60
    res@cnMaxLevelValF       =   60
    res@cnLevelSpacingF      =  10

    res@gsnLeftString = "(d)linearly congruent with Pr"

    plot(3) = gsn_csm_vector_scalar_map(wks,uwnd_p_rc,vwnd_p_rc,slp_p_rc ,res)
  ;;;;;



  ;;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = ispan(1,8,1) + ""
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   do i = 0,7
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do
;;;;


;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(8,graphic)

    iplot(0)=gsn_csm_contour(wks,sst_info[2],res_pro)
    iplot(1)=gsn_csm_contour(wks,sst_info[3],res_pro)
    iplot(2)=gsn_csm_contour(wks,slp_info[2],res_pro)
    iplot(3)=gsn_csm_contour(wks,slp_info[3],res_pro)

    iplot(4)=gsn_csm_contour(wks,m_zonal_info[2],res_pro)
    iplot(5)=gsn_csm_contour(wks,m_zonal_info[3],res_pro)
    
    iplot(6)=gsn_csm_contour(wks,m_meri_info[2],res_pro)
    iplot(7)=gsn_csm_contour(wks,m_meri_info[3],res_pro)
    

    do i = 0,7
    overlay(plot(i), iplot(i))
    end do

    pres = True
    pres@gsnPanelLabelBar = False
    gsn_panel(wks, plot, (/4,2/), pres)
end
    







    
   



  