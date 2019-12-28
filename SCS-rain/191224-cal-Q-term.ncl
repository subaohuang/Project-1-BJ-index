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

    era_dQdx_flt  = era_wgt_qu
    era_dQdx   = flt2dble(era_dQdx_flt)
    era_dQdy_flt  = era_wgt_qv
    era_dQdy   = flt2dble(era_dQdy_flt)


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
    ; era_divQ_bar = month_to_season(DivQ, "MAM")

;==============================================================
; Mission 2. calculate individual term of Q
; ==============================================================
   
   ;;; Calculate the MFC_advection term
   ;; MFC_advect = -(u*(dq/dx)+v*(dq/dy) ) 

   q_bar = month_to_annual(q, 1)
   u_bar = month_to_annual(u, 1)
   v_bar = month_to_annual(v, 1)

   ; q_bar = month_to_season(q, "MAM")
   ; u_bar = month_to_season(u, "MAM")
   ; v_bar = month_to_season(v, "MAM") 



   printVarSummary(q_bar)
   
   dq = grad_latlon_cfd(q_bar, g_lat, g_lon, False, False)

   era_adv_q_bar = q_bar                         
   era_adv_q_bar = u_bar*dq[1] + v_bar*dq[0]    ;;dq[o]:dq/dy，dq[1]:dq/dx

   dp_cfm_new = conform_dims(dimsizes(u_bar),dp,1)          ; confirm dp to caculate sum(var*dp)
   era_wgt_adv_q_bar=(/dim_sum_n(dp_cfm_new*era_adv_q_bar,1)/)               ; sum(var*dp)
   era_wgt_adv_q_bar = (/era_wgt_adv_q_bar*86400/(-9.8)/)       ;;convert to day
   copy_VarCoords(q_bar(:,1,:,:), era_wgt_adv_q_bar)

  ; Calculate the MFC_convergence term
  ;   MFC_conv   = -q*((du/dx)+(dv/dy) )      ; con(div)-vergence

    duv  = uv2dv_cfd(u_bar,v_bar,g_lat,g_lon,2)       ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)

    era_con_q_bar  = q_bar*duv
    era_wgt_con_q_bar=(/dim_sum_n(dp_cfm_new*era_con_q_bar,1)/)               ; sum(var*dp)
    
    era_wgt_con_q_bar=(/dim_sum_n(dp_cfm_new*era_con_q_bar,1)/)               ; sum(var*dp)
    era_wgt_con_q_bar = (/era_wgt_con_q_bar*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar)
   ;;;

    ; ;;;eddy term

    ; ;;;test 

    ; series_era_divQ_bar = dim_avg_n_Wrap(era_divQ_bar(:,{0:15},{110:155}), (/1,2/))
    ; series_adv_q_bar = dim_avg_n_Wrap(era_wgt_adv_q_bar(:,{0:15},{110:155}), (/1,2/))
    ; series_con_q_bar = dim_avg_n_Wrap(era_wgt_con_q_bar(:,{0:15},{110:155}), (/1,2/))
    
    ; wks   = gsn_open_wks("pdf", "Q_bar以及adve和convc的时间序列")
    ; res  = True
    ; res@xyLineColors       =  (/"black","blue","red","green","gray"/)          ; change line color
    ; res@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines"/)
    ; plot = gsn_csm_xy(wks, ispan(1,yearnum,1), (/series_era_divQ_bar,series_adv_q_bar,series_con_q_bar/), res)
    ; exit()

;==============================================================
; Mission 3. trend 
; ==============================================================

   
    divQ_rc  = regCoef_n(ispan(1, yearnum, 1), era_divQ_bar, 0, 0)
    copy_VarCoords(era_divQ_bar(1,:,:), divQ_rc)
    prob1 = rc_test(divQ_rc)
    divQ_rc = (/divQ_rc*yearnum/)

   
    era_wgt_adv_q_bar_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_adv_q_bar, 0, 0)
    copy_VarCoords(era_wgt_adv_q_bar(1,:,:), era_wgt_adv_q_bar_rc )


    prob2 = rc_test(era_wgt_adv_q_bar_rc)
    era_wgt_adv_q_bar_rc = (/era_wgt_adv_q_bar_rc*yearnum/)

    era_wgt_con_q_bar_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_con_q_bar, 0, 0)
    copy_VarCoords(era_wgt_con_q_bar(1,:,:), era_wgt_con_q_bar_rc )
    prob3 = rc_test(era_wgt_con_q_bar_rc)
    era_wgt_con_q_bar_rc = (/era_wgt_con_q_bar_rc*yearnum/)


    ;;;各成分的箭头
     era_wgt_qu_annual = month_to_annual(era_wgt_qu, 1)
     era_wgt_qv_annual = month_to_annual(era_wgt_qv, 1)

     era_wgt_qu_rc = regCoef_n(ispan(1, yearnum, 1), era_wgt_qu_annual, 0, 0)
     copy_VarCoords(era_divQ_bar(1,:,:), era_wgt_qu_rc)
     prob_qu = rc_test(era_wgt_qu_rc)
     era_wgt_qu_rc = (/era_wgt_qu_rc*yearnum/)


     era_wgt_qv_rc = regCoef_n(ispan(1, yearnum, 1), era_wgt_qv_annual, 0, 0)
     copy_VarCoords(era_divQ_bar(1,:,:),era_wgt_qv_rc )
     prob_qv = rc_test(era_wgt_qv_rc)
     era_wgt_qv_rc = (/era_wgt_qv_rc*yearnum/)

     ;;;



;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "Q_bar趋势以及adve和convc趋势-箭头")
  plot=new(3,graphic)
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
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -5
  res@mpMaxLatF             = 45
  res@mpMinLonF             = 90
  res@mpMaxLonF             = 160
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



  ; res@mpFillOn = False
  ; res@gsnScalarContour        = True
  ;   ;;;;;;;LineArrow
  ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; ; res@vcGlyphStyle="CurlyVector"

  ; res@vcMinDistanceF           = 0.02           ; thin vectors
  ; ; res@vcMinMagnitudeF          =0.04
  ; ; res@vcRefMagnitudeF          = 1.0
  ; ; res@vcRefMagnitudeF          = 0.1
  ; res@vcRefLengthF             = 0.08     ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  ; res@vcLineArrowThicknessF=1.0
  ; res@vcVectorDrawOrder="PostDraw"

  ; plot(0)=gsn_csm_vector_scalar_map(wks, era_wgt_qu_rc, era_wgt_qv_rc, divQ_rc, res)
   plot(0)=gsn_csm_contour_map(wks, divQ_rc, res)
  plot(1)=gsn_csm_contour_map(wks, era_wgt_adv_q_bar_rc, res)
  plot(2)=gsn_csm_contour_map(wks, era_wgt_con_q_bar_rc, res)

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
  
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   do i = 0,2
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
    iplot = new(3,graphic)

    iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob2,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob3,res_pro)
  

    ; do i=0,2
    ; overlay(plot(i), iplot(i))
    ; end do

 pres = True
 gsn_panel(wks, plot, (/3,1/), pres)
 

end 
    
  





    
   ;;;

   






