;;********************************************************
;;计算水汽通量算度，并且分解水汽通量散度
;; ERA-interim:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
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

function cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
begin
  ;;;;Q component
    Q_rc  = regCoef_n(ispan(1, yearnum, 1), era_divQ_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), Q_rc)
    Q_rc = (/Q_rc*yearnum/)
    
    era_wgt_qu_delta_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_qu_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), era_wgt_qu_delta_rc)
    era_wgt_qu_delta_rc = (/era_wgt_qu_delta_rc*yearnum/)

    era_wgt_qv_delta_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_qv_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), era_wgt_qv_delta_rc)
    era_wgt_qv_delta_rc = (/era_wgt_qv_delta_rc*yearnum/)
 return [/Q_rc,era_wgt_qu_delta_rc,era_wgt_qv_delta_rc/]
 end


begin
  ;;;read era data
		path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
		var_name = (/"uwnd","vwnd","shum"/)
		filein = path_wind + var_name +".mon.mean.nc"
        
        ;;;region 
          lat1 = -10
          lat2 = 50
          lon1 = 60
          lon2 = 160
        ;;;

        ; ;;;global region 
        ; lat1 = -50
        ; lat2 = 50
        ; lon1 = 0
        ; lon2 = 360
        ; ;;;

    

        ;;;time span
        startmon =197901
        endmon = 201312

        startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
        endyear_r = str_split_by_length(endmon, 4 )

        startyear = stringtoint(startyear_r(0))
        endyear = stringtoint(endyear_r(0))

        yearnum = endyear - startyear +1
        ;;


		f1 = addfile(filein(0),"r")

		date=cd_calendar(f1->time, 1)
		timeind=ind(date.ge.startmon.and.date.le.endmon)

		u = short2flt(f1->u(timeind,:,{lat1:lat2},{lon1:lon2}))

		f2 = addfile(filein(1),"r")
		v = short2flt(f2->v(timeind,:,{lat1:lat2},{lon1:lon2}))

		
		f3 = addfile(filein(2),"r")
		date :=cd_calendar(f3->time, 1)
		timeind :=ind(date.ge.startmon.and.date.le.endmon)
		q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
   ;;;
   print("=====finish read data==========")
;==============================================================
; Mission 1. calculate Q
; ==============================================================
    era_qu  = q
    era_qv  = q
    era_qu  = u*q
    printVarSummary(era_qu)
    era_qv  = v*q
    era_wgt_qu = era_qu(:,0,:,:)
    era_wgt_qv = era_qv(:,0,:,:)


   ;;;------Cal DivQ-------
    ;-----------Follow is the vertical integral-----------
    pver=era_qu&level
    pver=pver*100
    pver@units="Pa"

    dp=dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(era_qu),dp,1)          ; confirm dp to caculate sum(var*dp)

    era_wgt_qu=(/dim_sum_n(dp_cfm*era_qu,1)/)               ; sum(var*dp)
    era_wgt_qv=(/dim_sum_n(dp_cfm*era_qv,1)/)               ; sum(var*dp)

    era_wgt_qu = (/era_wgt_qu/9.8/)
    era_wgt_qv = (/era_wgt_qv/9.8/)

    g_lon   =   era_wgt_qu&longitude
    g_lat   =   era_wgt_qu&latitude
    g_nlat  =   dimsizes(g_lat)

    DivQ = uv2dv_cfd(era_wgt_qu,era_wgt_qv,g_lat,g_lon,2)
    ; DivQ = uv2dvF_Wrap(era_wgt_qu(:,::-1,:),era_wgt_qv(:,::-1,:)) 
    DivQ   =   (/-DivQ*86400/) 
    copy_VarCoords(era_wgt_qu, DivQ)
    printVarSummary(DivQ)
    ; print(DivQ(0,:,:))
    print("=====finish DivQ calculate==========")
    
    era_divQ_bar = month_to_annual(DivQ, 1)
    era_divQ_delta = dim_rmvmean_n_Wrap(era_divQ_bar, 0)
    
    printVarSummary(era_divQ_delta)

    ;;;得到divQ_bar的两个箭头
      era_wgt_qu_annual = month_to_annual(era_wgt_qu, 1)
      era_wgt_qu_delta =  dim_rmvmean_n_Wrap(era_wgt_qu_annual, 0)
      
      era_wgt_qv_annual = month_to_annual(era_wgt_qv, 1)
      era_wgt_qv_delta =  dim_rmvmean_n_Wrap(era_wgt_qv_annual, 0)

      era_wgt_qu_delta = (/era_wgt_qu_delta*9.8/)
      era_wgt_qv_delta = (/era_wgt_qv_delta*9.8/)

    ;;;


;==============================================================
; Mission 2. calculate dynamic and thermaldynamic term of Q
; ==============================================================
   
   ;;;prepare variable

   
   q_bar = month_to_annual(q, 1)
   u_bar = month_to_annual(u, 1)
   v_bar = month_to_annual(v, 1)

   q_delta = dim_rmvmean_n_Wrap(q_bar, 0)
   u_delta = dim_rmvmean_n_Wrap(u_bar, 0)
   v_delta = dim_rmvmean_n_Wrap(v_bar, 0)
   
   q_clm = dim_avg_n_Wrap(q_bar, 0)
   u_clm = dim_avg_n_Wrap(u_bar, 0)   
   v_clm = dim_avg_n_Wrap(v_bar, 0)
  
   q_clm_cfm = conform_dims(dimsizes(q_bar),q_clm,(/1,2,3/))
   u_clm_cfm = conform_dims(dimsizes(u_bar),u_clm,(/1,2,3/))          
   v_clm_cfm = conform_dims(dimsizes(v_bar),v_clm,(/1,2,3/))

   
   dim = dimsizes(q_bar(1,:,:,:))
   q_mon = q
   u_mon = u
   v_mon = v
   do i = 0,yearnum-1,1
   q_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), q_bar(i,:,:,:), (/1,2,3/))
   u_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), u_bar(i,:,:,:), (/1,2,3/))
   v_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), v_bar(i,:,:,:), (/1,2,3/))
   end do 
   
    q_eddy = q 
    q_eddy = (/q - q_mon/)     ;;;;eddy term
    u_eddy = u 
    u_eddy = (/u - u_mon/)     ;;;;eddy term
    v_eddy = v 
    v_eddy = (/v - v_mon/)     ;;;;eddy term

   
    
    uq_eddy = q_eddy*u_eddy
    copy_VarCoords(q_eddy, uq_eddy)

   
    
    vq_eddy = q_eddy*v_eddy
    copy_VarCoords(q_eddy, vq_eddy)
    
    

    uq_eddy_bar = month_to_annual(uq_eddy, 1)
    vq_eddy_bar = month_to_annual(vq_eddy, 1)
    


    uq_eddy_delta = dim_rmvmean_n_Wrap(uq_eddy_bar, 0)
    vq_eddy_delta = dim_rmvmean_n_Wrap(vq_eddy_bar, 0)

    print("=====end variable set==========")
   


  
  ;;;;;
    print("=====now start componte calculate==========")
   ;;; Calculate the MFC dynamic term
   ;; MFC_dynamic = d(δv*q)  

    era_dyna_qu = u_delta*q_clm_cfm
    era_dyna_qv = v_delta*q_clm_cfm
   

    era_dynamic_q = uv2dv_cfd(era_dyna_qu,era_dyna_qv,g_lat,g_lon,2)
    dp_cfm_new = conform_dims(dimsizes(u_bar),dp,1)          ; confirm dp to caculate sum(var*dp)

    era_wgt_dyna_q = (/dim_sum_n(dp_cfm_new*era_dynamic_q,1)/)               ; sum(var*dp)
    era_wgt_dyna_q = (/era_wgt_dyna_q*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_dyna_q)
    
    ;;;得到他们的箭头项
    era_wgt_dyna_qu = (/dim_sum_n(dp_cfm_new*era_dyna_qu,1)/)               ; sum(var*dp)
    era_wgt_dyna_qv = (/dim_sum_n(dp_cfm_new*era_dyna_qv,1)/)               ; sum(var*dp)   

    copy_VarCoords(q_bar(:,1,:,:), era_wgt_dyna_qu)
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_dyna_qv)
 
  

    ;;;;

    


    print("=====end dynamic componte calculate==========")
  
  ; Calculate the MFC thermaldynamic term
  ;   MFC_thermaldynamic = d(v*δq)
    
    era_thermaldyna_qu = u_clm_cfm*q_delta
    era_thermaldyna_qv = v_clm_cfm*q_delta
   
    era_thermaldyna_q = uv2dv_cfd(era_thermaldyna_qu,era_thermaldyna_qv,g_lat,g_lon,2)
    
    era_wgt_thermaldyna_q = (/dim_sum_n(dp_cfm_new*era_thermaldyna_q,1)/)               ; sum(var*dp)
    era_wgt_thermaldyna_q = (/era_wgt_thermaldyna_q*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_thermaldyna_q)

    ;;;得到他们的箭头项
    era_wgt_thermaldyna_qu = (/dim_sum_n(dp_cfm_new*era_thermaldyna_qu,1)/)               ; sum(var*dp)
    era_wgt_thermaldyna_qv = (/dim_sum_n(dp_cfm_new*era_thermaldyna_qv,1)/)               ; sum(var*dp)   

    ; era_wgt_thermaldyna_qu = (/era_wgt_thermaldyna_qu*86400/9.8/)
    ; era_wgt_thermaldyna_qv = (/era_wgt_thermaldyna_qv*86400/9.8/)

    copy_VarCoords(q_bar(:,1,:,:), era_wgt_thermaldyna_qu)
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_thermaldyna_qv)
    ;;;;
     print("=====end thermaldynamic componte calculate==========")


  ;;;quadratic term d(δv*δq)

    era_quad_qu = u_delta*q_delta
    era_quad_qv = v_delta*q_delta

    era_quad_q = uv2dv_cfd(era_quad_qu,era_quad_qv,g_lat,g_lon,2)
    
    era_wgt_quad_q = (/dim_sum_n(dp_cfm_new*era_quad_q,1)/)               ; sum(var*dp)
    era_wgt_quad_q = (/era_wgt_quad_q*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_quad_q)

     ;;;得到他们的箭头项
    era_wgt_quad_qu = (/dim_sum_n(dp_cfm_new*era_quad_qu,1)/)               ; sum(var*dp)
    era_wgt_quad_qv = (/dim_sum_n(dp_cfm_new*era_quad_qv,1)/)               ; sum(var*dp)   


    copy_VarCoords(q_bar(:,1,:,:), era_wgt_quad_qu)
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_quad_qv)
    ;;;;
     print("=====end quadratic componte calculate==========")
  ;;;



  ;;;quandratic climate term d((δv*δq)clm)

     era_quad_qu_clm = dim_avg_n_Wrap(era_quad_qu, 0)
     era_quad_qv_clm = dim_avg_n_Wrap(era_quad_qv, 0)

     era_quad_clm_q = uv2dv_cfd(era_quad_qu_clm ,era_quad_qv_clm,g_lat,g_lon,2)

      dp_cfm_new2 = conform_dims(dimsizes(u_bar(1,:,:,:)),dp,0)          ; confirm dp to caculate sum(var*dp)
      era_wgt_quad_q_clm = (/dim_sum_n(dp_cfm_new2*era_quad_clm_q,0)/)               ; sum(var*dp)
      era_wgt_quad_q_clm = (/era_wgt_quad_q_clm*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(1,1,:,:), era_wgt_quad_q_clm)
      printVarSummary(era_wgt_quad_q_clm)
     
       print("=====end quandratic clmate componte calculate==========")

  ;;;eddy term δd((bar(v'q'))

    era_eddy_q = uv2dv_cfd(uq_eddy_delta,vq_eddy_delta,g_lat,g_lon,2)

    era_wgt_eddy_q = (/dim_sum_n(dp_cfm_new*era_eddy_q,1)/)               ; sum(var*dp)
    era_wgt_eddy_q = (/era_wgt_eddy_q*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_eddy_q)


    ;;;得到他们的箭头项
    era_wgt_eddy_qu = (/dim_sum_n(dp_cfm_new*uq_eddy_delta,1)/)               ; sum(var*dp)
    era_wgt_eddy_qv = (/dim_sum_n(dp_cfm_new*vq_eddy_delta,1)/)               ; sum(var*dp)   


    copy_VarCoords(q_bar(:,1,:,:), era_wgt_eddy_qu)
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_eddy_qv)
    ;;;;
     print("=====end eddy componte calculate==========")



  ;;;


; ;==============================================================
; ; Mission 3.  trend
; ; ==============================================================
    ;;;;Q component
    Q_rc =  cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
    
    q_rc = Q_rc[0]
    qu_rc  =  Q_rc[1]
    qv_rc  =  Q_rc[2]

    prob_q = rc_test(q_rc)
    prob_qu = rc_test(qu_rc)
    prob_qv = rc_test(qv_rc)
    ;;;;
    Q_dyna_rc  =  cal_trend(yearnum,era_wgt_dyna_q,era_wgt_dyna_qu,era_wgt_dyna_qv)    
    q_dyna_rc = Q_dyna_rc[0]
    qu_dyna_rc  =  Q_dyna_rc[1]
    qv_dyna_rc  =  Q_dyna_rc[2]

    prob_q_dyna = rc_test(q_dyna_rc)
    prob_qu_dyna = rc_test(qu_dyna_rc)
    prob_qv_dyna = rc_test(qv_dyna_rc)
    ;;;;

    Q_thermaldyna_rc  =  cal_trend(yearnum,era_wgt_thermaldyna_q,era_wgt_thermaldyna_qu,era_wgt_thermaldyna_qv)    
    q_thermaldyna_rc = Q_thermaldyna_rc[0]
    qu_thermaldyna_rc  =  Q_thermaldyna_rc[1]
    qv_thermaldyna_rc  =  Q_thermaldyna_rc[2]

    prob_q_thermaldyna = rc_test(q_thermaldyna_rc)
    prob_qu_thermaldyna = rc_test(qu_thermaldyna_rc)
    prob_qv_thermaldyna = rc_test(qv_thermaldyna_rc)
    ;;;

    Q_quad_rc  =  cal_trend(yearnum,era_wgt_quad_q,era_wgt_quad_qu,era_wgt_quad_qv)    
    q_quad_rc = Q_quad_rc[0]
    qu_quad_rc  =  Q_quad_rc[1]
    qv_quad_rc  =  Q_quad_rc[2]

    prob_q_quad = rc_test(q_quad_rc)
    prob_qu_quad = rc_test(qu_quad_rc)
    prob_qv_quad = rc_test(qv_quad_rc)
   
    ;;;;;;;;
    Q_eddy_rc  =  cal_trend(yearnum,era_wgt_eddy_q,era_wgt_eddy_qu,era_wgt_eddy_qv)    
    q_eddy_rc = Q_eddy_rc[0]
    qu_eddy_rc  =  Q_eddy_rc[1]
    qv_eddy_rc  =  Q_eddy_rc[2]

    prob_q_eddy = rc_test(q_eddy_rc)
    prob_qu_eddy = rc_test(qu_eddy_rc)
    prob_qv_eddy = rc_test(qv_eddy_rc)
    ;;;
    
    ; qu_rc = where(prob_qu.le.0.2, qu_rc,qu_rc@_FillValue)
    ; qv_rc = where(prob_qv.le.0.2, qv_rc,qv_rc@_FillValue)


;==============================================================
; Mission 6. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "q_delta以及各个项的趋势以及箭头")
  plot=new(5,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="MPL_RdBu"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -4.0
  res@cnMaxLevelValF       =   4.0
  res@cnLevelSpacingF      =  0.5
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -5
  res@mpMaxLatF             = 45
  res@mpMinLonF             = 90
  res@mpMaxLonF             = 160
  
  ; ;;;global region
  ;  res@mpMinLatF             = -30
  ;  res@mpMaxLatF             = 30
  ;  res@mpMinLonF             = 0
  ;  res@mpMaxLonF             = 360
  ; ;;;

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

  res@mpFillOn = False
  res@gsnScalarContour        = True
    ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.02           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  ; res@vcRefMagnitudeF          = 1.0
 
  res@vcRefLengthF             = 0.08     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.0
  res@vcVectorDrawOrder="PostDraw"
  res@vcRefMagnitudeF          = 1000.0
  res@gsnLeftString = "Moisture Convergence Flux"
  plot(0)=gsn_csm_vector_scalar_map(wks,qu_rc ,qv_rc, q_rc, res)
  res@gsnLeftString = "Dynamic component"
    res@vcRefMagnitudeF          = 1000.0
  plot(1)=gsn_csm_vector_scalar_map(wks,qu_dyna_rc ,qv_dyna_rc, q_dyna_rc, res)
  res@gsnLeftString = "Thermaldynamic component"
    res@vcRefMagnitudeF          = 200.0
  plot(2)=gsn_csm_vector_scalar_map(wks,qu_thermaldyna_rc ,qv_thermaldyna_rc, q_thermaldyna_rc, res)
  res@gsnLeftString = "Quadratic component"
    res@vcRefMagnitudeF          = 30.0
  plot(3)=gsn_csm_vector_scalar_map(wks,qu_quad_rc ,qv_quad_rc, q_quad_rc, res)
  res@gsnLeftString = "Transient component"
    res@vcRefMagnitudeF          = 250.0
  plot(4)=gsn_csm_vector_scalar_map(wks,qu_eddy_rc ,qv_eddy_rc, q_eddy_rc, res)
 

  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = "" + ispan(1,5,1) + ""
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   do i = 0,4
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False

    ; res_pro@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res_pro@cnMinLevelValF       =   0
    ; res_pro@cnMaxLevelValF       =   0.05
    ; res_pro@cnLevelSpacingF      =  0.05

    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(5,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_q,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob_q_dyna,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob_q_thermaldyna,res_pro)
    iplot(3)=gsn_csm_contour(wks,prob_q_quad,res_pro)
    iplot(4)=gsn_csm_contour(wks,prob_q_eddy,res_pro)
  

    do i=0,4
    overlay(plot(i), iplot(i))
    end do

 pres = True
 pres@gsnPanelLabelBar = True
 pres@lbOrientation = "vertical"
 pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)

 gsn_panel(wks, plot, (/3,2/), pres)
 

end 
    
   
   

  





    
   ;;;

   






