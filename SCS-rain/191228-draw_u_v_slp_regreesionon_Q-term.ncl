;;********************************************************
;;回归ERA的u,v,slp,z到水汽收支方程的不同的变化上
;;
;;lnshuheng 2019.12.24
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
;;;读入Q-term
    f3 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Q-term-5.nc", "r")
    series_Q = f3->series_Q(:)
    series_thermaldyna_q = f3->series_thermaldyna_q
    series_dyna_q = f3->series_dyna_q
    series_quad_q = f3->series_quad_q
    series_eddy_q = f3->series_eddy_q
;;;

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
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    
    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    slp = short2flt(f3->msl(timeind,{lat1:lat2},{lon1:lon2}))        

    f4 = addfile(filein(3),"r")
    date :=cd_calendar(f4->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    z = short2flt(f4->z(timeind,{500},{lat1:lat2},{lon1:lon2}))   


    f5 = addfile(filein(4),"r")
    date :=cd_calendar(f5->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    q = short2flt(f5->q(timeind,:,{lat1:lat2},{lon1:lon2}))   


    ;;;积分q
    pver = q&level
    era_wgt_q = q(:,1,:,:)
    dp=dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(q),dp,1)          ; confirm dp to caculate sum(var*dp)

    era_wgt_q=(/dim_sum_n(dp_cfm*q,1)/)               ; sum(var*dp)

    ; slp = (/era_wgt_q/)
    ;;;;

    f6 = addfile(filein(5),"r")
    date :=cd_calendar(f6->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    t2m = short2flt(f6->t2m(timeind,{lat1:lat2},{lon1:lon2}))   

    ; slp = (/t2m/)
    ;;;;

    ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)

    sst=lonFlip(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    sst_annual = month_to_annual(ssta, 1)
    sst_annual = (/dim_standardize(sst_annual, 1)/)
    print("=====finish read data==========")
   ;;;
 
    
   ;;;处理成年异常数据
     u = (/rmMonAnnCycTLL(u)/)
     v = (/rmMonAnnCycTLL(v)/)
     z = (/rmMonAnnCycTLL(z)/)
     slp = (/rmMonAnnCycTLL(slp)/)

     u_annual = month_to_annual(u, 1)
     v_annual = month_to_annual(v, 1)
     z_annual = month_to_annual(z, 1)
     slp_annual = month_to_annual(slp, 1)

     u_annual = (/dtrend_msg_n(ispan(1,dim(0),1),u_annual,False,True,0)/)
     v_annual = (/dtrend_msg_n(ispan(1,dim(0),1),v_annual,False,True,0)/)
     z_annual = (/dtrend_msg_n(ispan(1,dim(0),1),z_annual,False,True,0)/)
     slp_annual = (/dtrend_msg_n(ispan(1,dim(0),1),slp_annual,False,True,0)/)


   ;;;回归前处理

    series_dyna_q = (/dtrend_msg_n(ispan(1,dim(0),1),series_dyna_q,False,True,0)/)
    series_thermaldyna_q = (/dtrend_msg_n(ispan(1,dim(0),1),series_thermaldyna_q,False,True,0)/)
    series_quad_q = (/dtrend_msg_n(ispan(1,dim(0),1),series_quad_q,False,True,0)/)
    series_eddy_q = (/dtrend_msg_n(ispan(1,dim(0),1),series_eddy_q,False,True,0)/)
   
    series_dyna_q = (/dim_standardize(series_dyna_q, 1)/)
    series_thermaldyna_q = (/dim_standardize(series_thermaldyna_q, 1)/)
    series_quad_q = (/dim_standardize(series_quad_q, 1)/)
    series_eddy_q = (/dim_standardize(series_eddy_q, 1)/)



  
   ;;回归
    u_dyna_rc = regCoef_n(series_dyna_q, u_annual, 0, 0)
    v_dyna_rc = regCoef_n(series_dyna_q, v_annual, 0, 0)
    slp_dyna_rc = regCoef_n(series_dyna_q, slp_annual, 0, 0)
    z_dyna_rc = regCoef_n(series_dyna_q, z_annual, 0, 0)
    
    copy_VarCoords(u_annual(1,:,:), u_dyna_rc)
    copy_VarCoords(v_annual(1,:,:), v_dyna_rc)
    copy_VarCoords(slp_annual(1,:,:), slp_dyna_rc)
    copy_VarCoords(z_annual(1,:,:), z_dyna_rc)

    u_quad_rc = regCoef_n(series_quad_q, u_annual, 0, 0)
    v_quad_rc = regCoef_n(series_quad_q, v_annual, 0, 0)
    slp_quad_rc = regCoef_n(series_quad_q, slp_annual, 0, 0)
    z_quad_rc = regCoef_n(series_quad_q, z_annual, 0, 0)
    
    copy_VarCoords(u_annual(1,:,:), u_quad_rc)
    copy_VarCoords(v_annual(1,:,:), v_quad_rc)
    copy_VarCoords(slp_annual(1,:,:), slp_quad_rc)
    copy_VarCoords(z_annual(1,:,:), z_quad_rc)

    u_thermaldyna_rc = regCoef_n(series_thermaldyna_q, u_annual, 0, 0)
    v_thermaldyna_rc = regCoef_n(series_thermaldyna_q, v_annual, 0, 0)
    slp_thermaldyna_rc = regCoef_n(series_thermaldyna_q, slp_annual, 0, 0)
    z_thermaldyna_rc = regCoef_n(series_thermaldyna_q, z_annual, 0, 0)
   
   copy_VarCoords(u_annual(1,:,:), u_thermaldyna_rc)
   copy_VarCoords(v_annual(1,:,:), v_thermaldyna_rc)
   copy_VarCoords(slp_annual(1,:,:), slp_thermaldyna_rc)
   copy_VarCoords(z_annual(1,:,:), z_thermaldyna_rc)
   
    u_eddy_rc = regCoef_n(series_eddy_q, u_annual, 0, 0)
    v_eddy_rc = regCoef_n(series_eddy_q, v_annual, 0, 0)
    slp_eddy_rc = regCoef_n(series_eddy_q, slp_annual, 0, 0)
    z_eddy_rc = regCoef_n(series_eddy_q, z_annual, 0, 0)

   copy_VarCoords(u_annual(1,:,:), u_eddy_rc)
   copy_VarCoords(v_annual(1,:,:), v_eddy_rc)
   copy_VarCoords(slp_annual(1,:,:), slp_eddy_rc)
   copy_VarCoords(z_annual(1,:,:), z_eddy_rc)
   ;;;


  
   ;;;显著性检验
    prob_slp_dyna = rc_test(slp_dyna_rc)
    prob_u_dyna = rc_test(u_dyna_rc)
    prob_v_dyna = rc_test(v_dyna_rc)
    prob_z_dyna = rc_test(z_dyna_rc)

    prob_slp_thermaldyna = rc_test(slp_thermaldyna_rc)
    prob_u_thermaldyna = rc_test(u_thermaldyna_rc)
    prob_v_thermaldyna = rc_test(v_thermaldyna_rc)
    prob_z_thermaldyna = rc_test(z_thermaldyna_rc)

    prob_slp_quad = rc_test(slp_quad_rc)
    prob_u_quad = rc_test(u_quad_rc)
    prob_v_quad = rc_test(v_quad_rc)
    prob_z_quad = rc_test(z_quad_rc)

    prob_slp_eddy = rc_test(slp_eddy_rc)
    prob_u_eddy = rc_test(u_eddy_rc)
    prob_v_eddy = rc_test(v_eddy_rc)
    prob_z_eddy = rc_test(z_eddy_rc)
   ;;;
    ; prob_slp = rtest(slp_rc,35,0)
    ; prob_u = rtest(u_rc,35,0)
    ; prob_v = rtest(v_rc,35,0)
    ; prob_z = rtest(z_rc,35,0)

   

    slp_dyna_rc = where(prob_slp_dyna.le.0.1, slp_dyna_rc, slp_dyna_rc@_FillValue)
    slp_thermaldyna_rc = where(prob_slp_thermaldyna.le.0.1, slp_thermaldyna_rc, slp_thermaldyna_rc@_FillValue)
    slp_quad_rc = where(prob_slp_quad.le.0.1, slp_quad_rc, slp_quad_rc@_FillValue)
    slp_eddy_rc = where(prob_slp_eddy.le.0.1, slp_eddy_rc, slp_eddy_rc@_FillValue)

   ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

   ; ;;;趋势一致性
   ; ;;;
   ; slp_rc = (/slp_rc*p_rc/)
   ; u_rc = (/u_rc*p_rc/)
   ; v_rc = (/v_rc*p_rc/)
   ; z_rc = (/z_rc*p_rc/)
    


;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "t2m_regression_on_Q-term")
  plot=new(4,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -50
  ; res@cnMaxLevelValF       =   50
  ; res@cnLevelSpacingF      =  10
  ;   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -4
  ; res@cnMaxLevelValF       =   4
  ; res@cnLevelSpacingF      =  1.0
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
  
  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.02           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  res@vcRefMagnitudeF          = 1.0
  ; res@vcRefMagnitudeF          = 0.1
  res@vcRefLengthF             = 0.08     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.0
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  res@vcVectorDrawOrder="PostDraw"
  res@mpFillOn =False

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
 
  plot(0) = gsn_csm_vector_scalar_map(wks,u_dyna_rc,v_dyna_rc,slp_dyna_rc ,res)

  plot(1) = gsn_csm_vector_scalar_map(wks,u_thermaldyna_rc,v_thermaldyna_rc,slp_thermaldyna_rc ,res)

  plot(2) = gsn_csm_vector_scalar_map(wks,u_quad_rc,v_quad_rc,slp_quad_rc ,res)
  
  plot(3) = gsn_csm_vector_scalar_map(wks,u_eddy_rc,v_eddy_rc,slp_eddy_rc ,res)

    str1 = unique_string("dum1")  ; "unique_string" will return a unique    

    ;;;画z的等直线
    res_z=True
    res_z@gsnDraw=False
    res_z@gsnFrame=False

    res_z@cnLinesOn=True
    res_z@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_z@cnMinLevelValF      = -10.0        ; set min contour level
    res_z@cnMaxLevelValF      = 10.0      ; set max contour level
    res_z@cnLevelSpacingF     = 1      ; set contour spacin

     res_z@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_z@cnMinLevelValF      = -50.0        ; set min contour level
    res_z@cnMaxLevelValF      = 50.0      ; set max contour level
    res_z@cnLevelSpacingF     = 5      ; set contour spacin
    res_z@cnLineLabelsOn     =True
    ; res_z@cnLineLabelBackgroundColor="Transparent"

    res_z@gsnContourNegLineDashPattern=2
    res_z@gsnContourPosLineDashPattern = 0

    res_z@cnLineLabelInterval        = 2             ; default = 2
    res_z@cnInfoLabelOn=False
    res_z@cnLineLabelFont=3.0
    res_z@cnLineLabelFontHeightF=0.008
    res_z@cnInfoLabelOn=False
    res_z@cnLineThicknessF=1.0
    res_z@gsnContourZeroLineThicknessF=0.001

    ; plot_z = gsn_csm_contour(wks , z_rc ,res_z)

    ; overlay(plot, plot_z)


  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3","4"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   do i = 0,3
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do
   ; ;;;打点
   ;  res_pro=True
   ;  res_pro@gsnAddCyclic = False
   ;  res_pro@gsnDraw=False
   ;  res_pro@gsnFrame=False
   ;  res_pro@cnFillOn = True

   ;  res_pro@cnLevelSelectionMode = "ExplicitLevels"
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  res_pro@lbLabelBarOn = False
   ;  res_pro@cnLinesOn=False
   ;  ; res_pro@cnMonoLineColor =False
   ;  ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

   ;  res_pro@cnMonoFillPattern = True
   ;  res_pro@cnFillPattern = 17
   ;  res_pro@cnFillDotSizeF = 0.0045
   ;  res_pro@cnMonoFillColor = False
   ;  res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
   ;  ; res_pro@cnLevels = ispan(-6,6,1)
   ;  res_pro@cnLineLabelsOn=False
   ;  res_pro@tmYRMode = "Automatic"
   ;  res_pro@cnInfoLabelOn=False
   ;  ; res_pro@cnLineLabelBackgroundColor="Transparent"
   ;  ; res_pro@cnLineThicknessF = 3.
   ;  ; res_pro@cnFillOn        = True
    ; iplot = new(4,graphic)

   ;  iplot(0)=gsn_csm_contour(wks,prob_slp,res_pro)
    ; do i = 0,3
    ; overlay(plot(i), iplot(i))
    ; end do

    pres = True
    pres@gsnPanelLabelBar = True
    gsn_panel(wks, plot, (/2,2/), pres)
end
    
  







    
   