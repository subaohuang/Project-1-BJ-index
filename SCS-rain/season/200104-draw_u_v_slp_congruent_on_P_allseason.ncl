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

begin
;;;读入GPCP降水场
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.201612)
    p_t=f2->precip(timeind,:,:)
;;;;



  
;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"uwnd","vwnd","slp","hgt","shum","2_metre_temperature"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
        ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
    ;;;

    ;;;time span
    startmon =197901
    endmon = 201612

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;


    f1 = addfile(filein(0),"r")

    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)

    u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    
    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    slp = short2flt(f3->msl(timeind,{lat1:lat2},{lon1:lon2}))        

    print("==========finish read data =============")
    
   ;;;处理成年异常数据
     
     u = (/rmMonAnnCycTLL(u)/)
     v = (/rmMonAnnCycTLL(v)/)
     slp = (/rmMonAnnCycTLL(slp)/)
     p_t = (/rmMonAnnCycTLL(p_t)/)

     p = dim_avg_n_Wrap(p_t(:,{0:15},{110:155}), (/1,2/))

     season = (/"MAM","JJA","SON","DJF"/)

     wks   = gsn_open_wks("pdf", "U,V,SLP_congruent_on_P-allseason")

     plot=new(8,graphic)
     str = ispan(1, 8, 1) +""

  do i=0,3
     p_annual = month_to_season(p, season(i))
     u_annual = month_to_season(u, season(i))
     v_annual = month_to_season(v, season(i))
     slp_annual = month_to_season(slp, season(i))

     p_annual = (/dim_standardize(p_annual, 1)/)


     ;;;计算趋势,计算线性一致性之前先把趋势计算了
     
     u_rc  = regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
     u_rc = (/u_rc*yearnum/)
     copy_VarCoords(u_annual(1,:,:), u_rc)
     
     v_rc  = regCoef_n(ispan(1, yearnum, 1), v_annual, 0, 0)
     v_rc = (/v_rc*yearnum/)
     copy_VarCoords(v_annual(1,:,:), v_rc)
     
     slp_rc  = regCoef_n(ispan(1, yearnum, 1), slp_annual, 0, 0)
     slp_rc = (/slp_rc*yearnum/)
     copy_VarCoords(slp_annual(1,:,:), slp_rc)

     p_rc  = regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
     p_rc = (/p_rc*yearnum/)

    
     u_annual = (/dtrend_msg_n(ispan(1,dim(0),1),u_annual,False,True,0)/)
     v_annual = (/dtrend_msg_n(ispan(1,dim(0),1),v_annual,False,True,0)/)
     slp_annual = (/dtrend_msg_n(ispan(1,dim(0),1),slp_annual,False,True,0)/)
     p_annual = (/dtrend_msg_n(ispan(1,dim(0),1),p_annual,False,True,0)/)


   
   ;;回归
    ;;1: 回归到P
    u_p_rc = regCoef_n(p_annual, u_annual, 0, 0)
    v_p_rc = regCoef_n(p_annual, v_annual, 0, 0)
    slp_p_rc = regCoef_n(p_annual, slp_annual, 0, 0)

    copy_VarCoords(u_annual(1,:,:), u_p_rc)
    copy_VarCoords(v_annual(1,:,:), v_p_rc)
    copy_VarCoords(slp_annual(1,:,:), slp_p_rc)
  
   ;;;

   ;;;显著性检验
    
    prob_slp = rc_test(slp_rc)
    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)


    prob_slp_p = rc_test(slp_p_rc)
    prob_u_p = rc_test(u_p_rc)
    prob_v_p = rc_test(v_p_rc)
   ;;;;
   
   
    ;;;
    slp_rc = where(prob_slp.le.0.1, slp_rc, slp_rc@_FillValue)
    slp_p_rc = where(prob_slp_p.le.0.1, slp_p_rc, slp_p_rc@_FillValue)
   
    

    u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)
    u_p_rc = where(prob_u_p.le.0.1.or.prob_v_p.le.0.1, u_p_rc, u_p_rc@_FillValue)
  
;;趋势一致性

    slp_p_rc = (/slp_p_rc*p_rc/)
    u_p_rc = (/u_p_rc*p_rc/)
    v_p_rc = (/v_p_rc*p_rc/)
;;;;
    
print("==========finish preprocessing data =============")

;==============================================================
; Mission 4. draw plot
; ==============================================================
     ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
    ;;;

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -160
  res@cnMaxLevelValF       =   160
  res@cnLevelSpacingF      =  40.0
  ;   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -4
  ; res@cnMaxLevelValF       =   4
  ; res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False


  if i.eq.3
   res@lbLabelBarOn=True
   res@lbOrientation = "Horizontal"
   res@pmLabelBarWidthF = 0.5
   res@pmLabelBarHeightF = 0.1
   res@pmLabelBarOrthogonalPosF = 0.3
  end if
  
  res@lbLabelStride=2
  ; res@lbOrientation="Vertical"

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

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
  res@gsnLeftString = "(a)"+season(i)+"-observed"+" U,V,SLP linear trend" 
  plot(i*2) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,slp_rc ,res)
  res@vcRefMagnitudeF          = 2.0
  res@cnMinLevelValF       =   -80
  res@cnMaxLevelValF       =   80
  res@cnLevelSpacingF      =  20

  res@gsnLeftString = "(b)"+season(i)+"-U,V,SLP linearly congruent with GPCP Pr"
    
  plot(i*2+1) = gsn_csm_vector_scalar_map(wks,u_p_rc,v_p_rc,slp_p_rc ,res)

  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   
   plot@$str(i*2)$ = gsn_add_polyline(wks, plot(i*2), pgx_pool, pgy_pool, pgres)
   plot@$str(i*2+1)$ = gsn_add_polyline(wks, plot(i*2+1), pgx_pool, pgy_pool, pgres)

   delete([/u_rc,v_rc,slp_rc,u_p_rc,v_p_rc,slp_p_rc/])
   delete([/lat1,lat2,lon1,lon2/])
end do

pres = True
pres@gsnPanelLabelBar = False
; pres@lbOrientation  = "Vertical"
gsn_panel(wks, plot, (/4,2/), pres)
  


end
    
  







    
    