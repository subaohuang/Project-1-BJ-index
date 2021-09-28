;;********************************************************
;;计算水汽通量算度，并且分解水汽通量散度
;; ERA-interim:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
;**********************************************************
load  "/home/ys17-19/lsh/Project/other/guanshun/add_TP.ncl"
  
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

;==============================================================
; Mission 0. plot(0) 
; ==============================================================
begin
    plot = new(4,graphic)
    iplot = plot

    str = ispan(1,4,1)+""

   
    ;;;time span
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;


    f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    ; f1=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201512)
    p=f1->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
  
    p_pi = p

    delete(p)

    p_annual = month_to_season(p_pi, "MAM")
    p_clm = dim_avg_n_Wrap(p_annual, 0)

    
    rc_mam  = regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
    copy_VarCoords(p_annual(1,:,:), rc_mam)
    rc_mam = (/rc_mam*yearnum/)

    prob_mam=rc_test(rc_mam)
    copy_VarCoords(p_pi(1,:,:), prob_mam)

    delete(p_pi)
    ; print(avg(rc_mam({0:15},{125:160})))
 
    ; print(avg(rc_mam({0:15},{125:160}))/avg(p_clm({0:15},{125:160})))




    ;;;
;;;;;
;==============================================================
; Mission read  wind data
; ==============================================================
  ;;;read era data
  path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
  var_name = (/"uwnd","vwnd","shum/shum.total"/)
  filein = path_wind + var_name +".mon.mean.nc"

  ;;;region 
    lat1 = -10
    lat2 = 50
    lon1 = 60
    lon2 = 180
  ;;;


    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; timeind :=ind(date.ge.197901.and.date.le.201512)

    yearnum = 37

    u = short2flt(f1->u(timeind,:,{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(1),"r")
    v = short2flt(f1->v(timeind,:,{lat1:lat2},{lon1:lon2}))


    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    timeind :=ind(date.ge.197901.and.date.le.201512)
    q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
    ;;;
print("=====finish read data==========")
;;;;

;;;calculate wind trend 
    u_annual = month_to_season(u(:,{850},:,:), "MAM")
    v_annual = month_to_season(v(:,{850},:,:), "MAM")

    u_rc  = regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
    copy_VarCoords(u_annual(1,:,:), u_rc)
    u_rc = (/u_rc*yearnum/)
    
    v_rc  = regCoef_n(ispan(1, yearnum, 1),v_annual, 0, 0)
    copy_VarCoords(u_annual(1,:,:), v_rc)
    v_rc = (/v_rc*yearnum/)

    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)

    u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)
    ; v_rc = where(prob_u.le.01.or.prob_v.le.0.1, v_rc, v_rc@_FillValue)
    u_rc_test = u_rc
    v_rc_test = v_rc
    u_rc_test = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

    printVarSummary(u_rc)

    delete(u_annual)
    delete(v_annual)

;;;   




;;;;
    wks   = gsn_open_wks("pdf", "fig2_review.pdf")
       
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic  = False

    res@cnFillPalette="NCV_blu_red"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=1
    res@pmLabelBarOrthogonalPosF = 0.17
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 30
    res@mpMinLonF             =100
    res@mpMaxLonF             = 170
    res@lbOrientation="horizontal"
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpFillOn = False
 

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025
    res@gsnLeftStringFontHeightF = 0.03
    res@tmYLMinorOn = True

    res@tmXBValues = (/100,120,140,160/) 
    res@tmXBLabels = (/100,120,140,160/) +""+"E" 


    res@tmYLValues = (/0,10,20,30/) 
    res@tmYLLabels = (/"0","10N","20N","30N"/) 

    res@tmYLMinorValues = ispan(-5,30,5)

    ; res@tmBorderThicknessF=1.0
    res@cnFillPalette = "precip_diff_12lev"
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ; res@cnLevels = (/-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/4,5,6,7,8,9,10,11,12,13,14/) - 2

    res@cnLevels = (/-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    res@cnFillColors= (/5,6,7,8,10,11,12,13,14/) - 2


    res@lbLabelFontHeightF = 0.02
    res@lbAutoManage = False

    res@gsnLeftString="(a) MAM Precip& 850UV trend " 
    res@gsnLeftStringOrthogonalPosF  = -0.12

    res@gsnScalarContour        = True

    ;;;;;;;LineArrow
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcMinDistanceF           = 0.025           ; thin vectors
    
    ; res@vcRefAnnoJust = "TopLeft"
    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False
    res@vcRefAnnoOrthogonalPosF = -1.2
    res@vcRefMagnitudeF          = 2.0
    res@vcMinMagnitudeF          =0.2
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006
    res@vcLineArrowHeadMaxSizeF = 0.018           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.012
    


    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "2 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013

    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False

    plot(0) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,rc_mam ,res)

    ; plot(0) = gsn_csm_contour_map(wks, rc_mam, res)
    

    ;;叠加矢量图
     ;;;;;;;LineArrow
    resv = True
    resv@gsnDraw = False
    resv@gsnFrame = False
    resv@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    resv@vcMinDistanceF           = 0.025           ; thin vectors
    ; res@vcMinMagnitudeF          =0.04
    resv@vcRefMagnitudeF          = 4.0
    resv@vcRefLengthF             = 0.1     ; define length of v
    resv@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    resv@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    resv@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    resv@vcLineArrowHeadMinSizeF = 0.006
    resv@vcVectorDrawOrder="PostDraw"

    resv@vcMonoLineArrowColor = True
    ; resv@vcLineArrowColor = "green"
    resv@vcLineArrowColor = "blue"
    resv@vcLineArrowThicknessF=1.5

    resv@vcRefAnnoOn    = False
    resv@gsnLeftString         = " "
    resv@gsnRightString        = " "


    plot_vect := gsn_csm_vector(wks,u_rc_test({-5:30},{100:170}),v_rc_test({-5:30},{100:170}),resv)
    ; overlay (plot(0),plot_vect)


    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=1
    ;;;画box
    lat1 = 0
    lat2 = 20
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    plot@$str(0)$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)


    ;;;打点
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False

    ; iplot(0) = gsn_csm_contour(wks,prob_mam,res_pro)    

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; ; opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeMid     = 17         ; stipple
    ; ; opt@gsnShadeFillDotSizeF =  0.0025
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數
    ; iplot(0) = gsn_contour_shade(iplot(0) , 0,0.1,opt)
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False


    res_pro@cnLinesOn=False


    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False

  
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)

    iplot=gsn_csm_contour(wks,prob_mam,res_pro)

    ; overlay(plot(i*3) , iplot(i*3) )
    overlay(plot(0), iplot(0))



;==============================================================
; Mission 1. calculate Q
; ==============================================================
    era_qu  = q
    era_qv  = q
    era_qu  = u*q
  
    era_qv  = v*q
    era_wgt_qu = era_qu(:,0,:,:)
    era_wgt_qv = era_qv(:,0,:,:)


  
   ;;;------Cal DivQ-------
    ;-----------Follow is the horizontal integral-----------
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
;;;;
;;;;
; ;;;choose season
   season = (/"MAM","JJA","SON","DJF"/)
   i = 0   
      era_divQ_bar = month_to_season(DivQ, season(i))
      era_divQ_delta = dim_rmvmean_n_Wrap(era_divQ_bar, 0)

      printVarSummary(era_divQ_delta)

      ;;;得到divQ_bar的两个箭头
      era_wgt_qu_annual = month_to_season(era_wgt_qu, season(i))
      era_wgt_qu_delta =  dim_rmvmean_n_Wrap(era_wgt_qu_annual, 0)


      era_wgt_qv_annual = month_to_season(era_wgt_qv, season(i))
      era_wgt_qv_delta =  dim_rmvmean_n_Wrap(era_wgt_qv_annual, 0)

      era_wgt_qu_delta = (/era_wgt_qu_delta*9.8/)
      era_wgt_qv_delta = (/era_wgt_qv_delta*9.8/)
      ;;;
  ;;;;

   
  ;;;;计算动力项和热力项

  ;;;;change δq ̅∇∙V ̅_clm+q ̅_clm ∇∙δV ̅+δV ̅∙∇q ̅_clm+V ̅_clm∙∇δq ̅

     q_bar = month_to_season(q, "MAM")
     u_bar = month_to_season(u, "MAM")
     v_bar = month_to_season(v, "MAM") 

     q_clm = dim_avg_n_Wrap(q_bar, 0)
     u_clm = dim_avg_n_Wrap(u_bar, 0)
     v_clm = dim_avg_n_Wrap(v_bar, 0)

     
     printVarSummary(q_bar)

     yearnum = 37


     q_change =  regCoef_n(ispan(1, yearnum, 1), q_bar, 0, 0)
     copy_VarCoords(q_bar(1,:,:,:), q_change)
     q_change = (/q_change*yearnum/)

     u_change =  regCoef_n(ispan(1, yearnum, 1), u_bar, 0, 0)
     copy_VarCoords(q_bar(1,:,:,:), u_change)
     u_change = (/u_change*yearnum/)

     v_change =  regCoef_n(ispan(1, yearnum, 1), v_bar, 0, 0)
     copy_VarCoords(q_bar(1,:,:,:), v_change)
     v_change = (/v_change*yearnum/)

   ;;;;
   

   ;;;;;
   ;;;σq▽(Vclm)

    duv_clm = u_clm
    duv_clm  = uv2dv_cfd(u_clm,v_clm,g_lat,g_lon,2)       ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)
    
    ;;; all level 
    era_con_q_change  = q_change*duv_clm

    dp_cfm_new2 = conform_dims(dimsizes(u_clm),dp,0)          ; confirm dp to caculate sum(var*dp)
    copy_VarCoords(u_clm, dp_cfm_new2)

    era_wgt_con_q_change = (/dim_sum_n(dp_cfm_new2*era_con_q_change,0)/)               ; sum(var*dp)
    era_wgt_con_q_change = (/era_wgt_con_q_change*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_clm(1,:,:), era_wgt_con_q_change)
    ;;;

  ;;;qclm▽(σV)

    duv_change = u_change
    duv_change  = uv2dv_cfd(u_change,v_change,g_lat,g_lon,2)       ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)
    
    ;;; all level 
    era_con_u_change = q_clm
    era_con_u_change  = duv_change*q_clm

    era_wgt_con_u_change = (/dim_sum_n(dp_cfm_new2*era_con_u_change,0)/)               ; sum(var*dp)
    era_wgt_con_u_change = (/era_wgt_con_u_change*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_clm(1,:,:), era_wgt_con_u_change)
    printVarSummary(dp_cfm_new2)
  ;;;;


  ;;;;σV▽(qclm)
   
    dq_clm = grad_latlon_cfd(q_clm, g_lat, g_lon, False, False)

    era_adv_u_change = q_clm                         
    era_adv_u_change = u_change*dq_clm[1] + v_change*dq_clm[0]    ;;dq[o]:dq/dy，dq[1]:dq/dx

    era_wgt_adv_u_change = (/dim_sum_n(dp_cfm_new2*era_adv_u_change,0)/)               ; sum(var*dp)
    era_wgt_adv_u_change = (/era_wgt_adv_u_change*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_clm(1,:,:), era_wgt_adv_u_change)
    

   ;;;;Vclm▽(σq)
   
    dq_change = grad_latlon_cfd(q_change, g_lat, g_lon, False, False)

    era_adv_q_change = u_clm                         
    era_adv_q_change = u_clm*dq_change[1] + v_clm*dq_change[0]    ;;dq[o]:dq/dy，dq[1]:dq/dx

    era_wgt_adv_q_change = (/dim_sum_n(dp_cfm_new2*era_adv_q_change,0)/)               ; sum(var*dp)
    era_wgt_adv_q_change = (/era_wgt_adv_q_change*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_clm(1,:,:), era_wgt_adv_q_change)


    delete(q)
    delete(u)
    delete(v)


; ;==============================================================
; ; Mission 3.  trend
; ; ==============================================================
    ;;;;Q component
    Q_rc =  cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
    yearnum = 37
    
    q_rc = Q_rc[0]
    qu_rc  =  Q_rc[1]
    qv_rc  =  Q_rc[2]

    prob_q = rc_test(q_rc)
    prob_qu = rc_test(qu_rc)
    prob_qv = rc_test(qv_rc)

    ; prob_q = 1- q_rc@prob
    ; prob_qu = 1- qu_rc@prob
    ; prob_qv = 1- qv_rc@prob

    copy_VarCoords(q_rc, prob_q)
    copy_VarCoords(q_rc, prob_qv)
    copy_VarCoords(q_rc, prob_qu)
    

    qu_rc_test = qu_rc
    qv_rc_test = qv_rc
    qu_rc_test = where(prob_qu.le.0.1.or.prob_qv.le.0.1, qu_rc,qu_rc@_FillValue)
    
    qu_rc = where(prob_qu.le.0.2.or.prob_qv.le.0.2, qu_rc,qu_rc@_FillValue)
;;;;;


;;;;计算相对变化量

    ; q_rc_avg = avg(dim_avg_n(q_rc({0:20},{125:160}),0))

    ; dyn_Q_trend_avg = avg(dim_avg_n(dyn_Q_trend({0:20},{125:160}),0))

    ; therdyn_Q_trend_avg = avg(dim_avg_n(therdyn_Q_trend({0:20},{125:160}),0))

    ; print("q_rc_avg="+q_rc_avg)
    ; print("dyn_Q_trend_avg="+dyn_Q_trend_avg)
    ; print("therdyn_Q_trend_avg="+therdyn_Q_trend_avg)
    percent = era_wgt_con_u_change
    percent = era_wgt_con_u_change / q_rc 
    printVarSummary(percent)
    print(avg(era_wgt_con_u_change({0:20},{125:160})))
    print(avg(q_rc({0:20},{125:160})))

    print(avg(era_wgt_con_u_change({0:20},{125:160}))/avg(q_rc({0:20},{125:160})))
    ; exit()

    print(avg(percent({0:20},{125:160})))

    ; print(percent({0:20},{125:160}))
    
    ; percent = dyn_Q_trend / q_rc 


;==============================================================
; Mission 6. draw plot
; ==============================================================
  colordata = read_colormap_file("/home/yangsong3/data-observation/linshh/data/colordata/bluewhiteorangered_11color.rgb")
  
  ; gsn_define_colormap(wks, colordata)
  res2=True
  res2@gsnFrame=False
  res2@gsnDraw=False
  res2@gsnAddCyclic = False

  res2@cnFillPalette=colordata(::-1,:)

  res2@cnFillOn=True
  res2@cnLinesOn=False

  res2@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res2@cnMinLevelValF       =   -3.0
  res2@cnMaxLevelValF       =   3.0
  res2@cnLevelSpacingF      =  1.0
  res2@cnLineLabelsOn=False
  res2@lbLabelBarOn=True
  res2@lbLabelStride=1
  res2@pmLabelBarOrthogonalPosF = 0.17

  res2@lbOrientation="horizontal"
  res2@mpLimitMode           = "LatLon"

  ;;;global region
  res2@mpMinLatF             = -5
  res2@mpMaxLatF             = 30
  res2@mpMinLonF             = 100
  res2@mpMaxLonF             = 170
  
  res2@mpCenterLonF          = (res2@mpMinLonF+res2@mpMaxLonF)/2.0
  res2@tmXTOn=False
  res2@tmYROn=False
  res2@tmXBOn=True
  res2@tmYROn=False

  res2@tmXBLabelFontHeightF=0.025
  res2@tmYLLabelFontHeightF=0.025
  res2@tmXBMajorThicknessF=2.0
  res2@tmYLMajorThicknessF=2.0
  ; res2@tmBorderThicknessF=1.0

  res2@gsnLeftStringFontHeightF = 0.03
  res2@tmXBValues = (/100,120,140,160/) 
  res2@tmXBLabels = (/100,120,140,160/) +""+"E" 

  res2@tmYLValues = (/0,10,20,30/) 
  res2@tmYLLabels = (/"0","10N","20N","30N"/) 

  res2@tmYLMinorValues = ispan(-5,30,5)



  res2@mpFillOn = False
  ; res2@gsnScalarContour        = True
  ; ;;;;;;;LineArrow
  ; res2@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; ; res2@vcGlyphStyle="CurlyVector"
  ; res2@vcMonoLineArrowColor = True
  ; res2@vcLineArrowColor = "black"
  ; res2@vcMinDistanceF           = 0.025           ; thin vectors

  ; res2@vcRefLengthF             = 0.1     ; define length of v
  ; res2@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  ; res2@vcVectorDrawOrder="PostDraw"

  ; res2@vcLineArrowThicknessF=1.5
  ; ; res2@vcLineArrowHeadMaxSizeF=0.01
  ; ; res2@vcLineArrowHeadMinSizeF=0.01
  ; ; res2@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  ; ; res2@vcLineArrowHeadMinSizeF = 0.006
  ; res2@vcLineArrowHeadMaxSizeF = 0.018           ;change the vector line arrows
  ; res2@vcLineArrowHeadMinSizeF = 0.012

  ; res2@vcRefMagnitudeF          = 1000.0
  ; res2@vcMinMagnitudeF          =100.0
  res2@gsnLeftString = "(d) MAM CVIMF Trend"
  ; res2@gsnLeftStringOrthogonalPosF  = -0.12


  ; res2@vcRefAnnoString1On = True
  ; res2@vcRefAnnoString1 = "1000 kg*m/s"
  ; res2@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
  ; res2@vcRefAnnoParallelPosF = 1.02
  ; res2@vcRefAnnoFontHeightF = 0.013

  ; res2@vcRefAnnoBackgroundColor = "Transparent"
  ; res2@vcRefAnnoPerimOn = False
  
  ; plot(3)=gsn_csm_vector_scalar_map(wks,qu_rc,qv_rc, q_rc, res2)
  plot(3)=gsn_csm_contour_map(wks, q_rc, res2)


  plot@$str(3)$ = gsn_add_polyline(wks, plot(3), pgx_pool, pgy_pool, pgres)
   

   ;;;打点
    res_pro2=True
    res_pro2@gsnAddCyclic = False
    res_pro2@gsnDraw=False
    res_pro2@gsnFrame=False
    res_pro2@cnFillOn = True

    res_pro2@cnLevelSelectionMode = "ExplicitLevels"
    res_pro2@cnLevels = (/0,0.1/)
    res_pro2@lbLabelBarOn = False


    res_pro2@cnLinesOn=False


    res_pro2@cnMonoFillPattern = True
    res_pro2@cnFillPattern = 17
    ; res_pro2@cnFillDotSizeF = 0.0045
    res_pro2@cnFillDotSizeF = 0.002
    res_pro2@cnMonoFillColor = False
    res_pro2@cnFillColors = (/"Transparent","White","Transparent"/)
  
    res_pro2@cnLineLabelsOn=False
    res_pro2@tmYRMode = "Automatic"
    res_pro2@cnInfoLabelOn=False

    iplot(3)=gsn_csm_contour(wks,prob_q,res_pro2)    
    overlay(plot(3), iplot(3))

    res_pro2@cnLevels = (/0.7,1.0/)
    res_pro2@cnFillPattern = 6
    res_pro2@cnFillColors = (/"Transparent","gray25","gray25"/)

    iplot_percent = gsn_csm_contour(wks,percent,res_pro2)
    overlay(plot(3), iplot_percent)
    
    res_pro2@cnLevels = (/0,0.1/)
    res_pro2@cnFillPattern = 17
    res_pro2@cnFillColors = (/"Transparent","White","Transparent"/)

 ;==============================================================
; Mission 3. calculate OLR
; ==============================================================
    f6 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
    ;;;region
    lat1 = -5
    lat2 = 45
    lon1 = 90
    lon2 = 170

    date:=cd_calendar(f6->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    olr = short2flt(f6->olr(timeind,{lat1:lat2},{lon1:lon2}))


    olr_annual = month_to_season(olr, season(i))
    printVarSummary(olr_annual)
    
    olr_rc = regCoef_n(ispan(1, yearnum, 1), olr_annual, 0, 0)
    copy_VarCoords(olr_annual(1,:,:),olr_rc)
    olr_rc = (/olr_rc*yearnum/)
    
    print(avg(olr_rc({0:15},{125:160})))
  


    prob_olr_rc = rc_test(olr_rc)
    ;;;; draw olr 
    res3=True
    res3@gsnFrame=False
    res3@gsnDraw=False
    res3@gsnAddCyclic = False

    res3@cnFillOn=True
    res3@cnLinesOn= False


    res3@cnLineLabelsOn=False
    res3@lbLabelBarOn=True
    res3@lbLabelStride=1
    res3@lbOrientation="horizontal"
    res3@pmLabelBarOrthogonalPosF = 0.17

   
    res3@mpLimitMode           = "LatLon"
    res3@mpMinLatF             = -5
    res3@mpMaxLatF             = 30
    res3@mpMinLonF             = 100
    res3@mpMaxLonF             = 170

    res3@tmXBValues = (/100,120,140,160/) 
    res3@tmXBLabels = (/100,120,140,160/) +""+"E" 

    res3@tmYLValues = (/0,10,20,30/) 
    res3@tmYLLabels = (/"0","10N","20N","30N"/) 

    res3@tmYLMinorValues = ispan(-5,30,5)

    res3@mpCenterLonF          = (res3@mpMinLonF+res3@mpMaxLonF)/2.0

    res3@tmXTOn=False
    res3@tmYROn=False
    res3@tmXBOn=True
    res3@tmYROn=False

    res3@tmXBLabelFontHeightF=0.025
    res3@tmYLLabelFontHeightF=0.025
    res3@tmXBMajorThicknessF=2.0
    res3@tmYLMajorThicknessF=2.0
    ; res3@tmBorderThicknessF=1.0

    res3@gsnScalarContour        = True
    res3@mpFillOn =False

    res3@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res3@cnMinLevelValF        := -20.0        ; set min contour level
    res3@cnMaxLevelValF        := 20.0          ; set max contour level
    res3@cnLevelSpacingF       := 5.0          ; set contour interval
    res3@cnFillPalette = "BlueWhiteOrangeRed"
    

    res3@cnLevelSelectionMode = "ExplicitLevels"

    ; gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    res3@cnLevels := (/-20.0,-15.0,-10.0,-5.0,-2.5,2.5,5.0,10.0,15.0,20.0/)
    res3@lbLabelStrings = res3@cnLevels
    res3@lbLabelFontHeightF = 0.02
    res3@lbAutoManage = False


    ; res3@lbLabelStrings = (/-20.0,-10.0,10.0,20.0/)


    ; res3@cnFillColors := (/2,5,30,40,3,146,177,200,251/) -4

    
    res3@gsnLeftString  ="(b) MAM OLR trend"
    res3@gsnLeftStringFontHeightF = 0.03
    res3@lbAutoManage = False


    plot(1) = gsn_csm_contour_map(wks,olr_rc,res3)


    plot@$str(1)$ = gsn_add_polyline(wks, plot(1), pgx_pool, pgy_pool, pgres)

  ;;;打点
    res_pro2@cnFillColors = (/"Transparent","White","Transparent"/)

    iplot(1)=gsn_csm_contour(wks,prob_olr_rc,res_pro2)

    overlay(plot(1), iplot(1))
  ;;;;
;==============================================================
; Mission 4. plot(3) 
; ==============================================================

    ;;;region 
    lat1 = -5
    lat2 = 45
    lon1 = 90
    lon2 = 170
    ;;;
    fils3 := systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
    f3 := addfiles(fils3, "r")

    date:=cd_calendar(f3[:]->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude
    lon = f3[1]->longitude
    level_era = f3[1]->level

    indwant = ind(level_era.le.1000.and.level_era.ge.850)
    levdim = dimsizes(indwant)
    
    div_r = new((/dimsizes(time),levdim,dimsizes(lat),dimsizes(lon)/),float)

    

    do i = 0,dimsizes(timeind)/12-1

       div_r(i*12:(i+1)*12-1,:,:,:) = short2flt(f3[i]->d(:,{1000:850},:,:))

    end do
    ;
    div = dim_avg_n_Wrap(div_r(:,:,{-10:30},{100:180}),1)

    ; div = div_r

    div = (/div*10^6/)
    div = (/-1*div/)
   
    w = div


    w_annual = month_to_season(w,"MAM")

    ; yearnum = 37
    wgt_w_rc := regCoef_n(ispan(1,yearnum,1), w_annual, 0, 0)
    copy_VarCoords(w_annual(1,:,:), wgt_w_rc)

    prob_w = rc_test(wgt_w_rc)
    wgt_w_rc = (/wgt_w_rc*yearnum/) 

    res4=True
    res4@gsnFrame=False
    res4@gsnDraw=False
    res4@gsnAddCyclic = False
  
    res4@cnFillOn=True
    res4@cnLinesOn=False

    res4@cnLineLabelsOn=False
    res4@lbLabelBarOn=True
    res4@lbLabelStride=1
    res4@lbOrientation="horizontal"
    res4@pmLabelBarOrthogonalPosF = 0.17


    res4@mpLimitMode           = "LatLon"
    res4@mpMinLatF             = -5
    res4@mpMaxLatF             = 30
    res4@mpMinLonF             = 100
    res4@mpMaxLonF             = 170 
    res4@mpCenterLonF          = (res4@mpMinLonF+res4@mpMaxLonF)/2.0
    res4@tmXTOn=False
    res4@tmYROn=False
    res4@tmXBOn=True
    res4@tmYROn=False
    res4@tmXBMode="Explicit"
    res4@tmYLMode="Explicit"


    res4@tmXBLabelFontHeightF=0.025
    res4@tmYLLabelFontHeightF=0.025
    res4@tmXBMajorThicknessF=2.0
    res4@tmYLMajorThicknessF=2.0
    ; res4@tmBorderThicknessF=1.0
    res4@tmXBValues = (/100,120,140,160/) 
    res4@tmXBLabels = (/100,120,140,160/) +""+"E" 
    res4@tmYLValues = (/0,10,20,30/) 
    res4@tmYLLabels = (/"0","10N","20N","30N"/) 

    res4@tmYLMinorValues = ispan(-5,30,5)

    res4@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
    res4@cnFillPalette = "MPL_PuOr"
    ; gsn_define_colormap(wks2, "MPL_PuOr")
    
    res4@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
    res4@lbLabelStrings = res4@cnLevels
    res4@lbLabelFontHeightF = 0.018
    ; res4@cnFillColors = (/22,32,45,60,-1,80,97,100,110/) -4

    
    res4@gsnLeftString = "(c) MAM low-level wind convergence trend "
    res4@gsnLeftStringFontHeightF = 0.03
    plot(2) = gsn_csm_contour_map(wks, wgt_w_rc, res4)


    plot@$str(2)$ = gsn_add_polyline(wks, plot(2), pgx_pool, pgy_pool, pgres)


  ;;;打点 ;;;打点
    res_pro4=True
    res_pro4@gsnAddCyclic = False
    res_pro4@gsnDraw=False
    res_pro4@gsnFrame=False
    res_pro4@cnFillOn = True

    res_pro4@cnLevelSelectionMode = "ExplicitLevels"
    res_pro4@cnLevels = (/0.0,0.1/)
    res_pro4@lbLabelBarOn = False

    res_pro4@cnLinesOn=False
    ; res_pro4@cnMonoLineColor =False
    ; res_pro4@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro4@cnMonoFillPattern = True
    res_pro4@cnFillPattern = 17
    res_pro4@cnFillDotSizeF = 0.001
    res_pro4@cnMonoFillColor = False
    res_pro4@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro4@cnLevels = ispan(-6,6,1)
    res_pro4@cnLineLabelsOn=False
    res_pro4@tmYRMode = "Automatic"
    res_pro4@cnInfoLabelOn=False



    iplot(2)=gsn_csm_contour(wks,prob_w,res_pro4)

    overlay(plot(2), iplot(2))
  ;;;;

;==============================================================
; Mission  panel  
; ==============================================================
pres = True

gsn_panel(wks, plot, (/2,2/), pres)


end