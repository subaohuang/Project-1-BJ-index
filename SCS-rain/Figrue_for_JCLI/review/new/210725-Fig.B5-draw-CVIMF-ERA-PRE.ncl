;==============================================================
; Mission read  wind data
; ==============================================================
function rc(data)
begin 
    dim = dimsizes(data)
    ; yearnum = 37
    wgt_Q_rc := regCoef_n(ispan(1,dim(0),1), data, 0, 0)
    copy_VarCoords(data(1,:,:), wgt_Q_rc)
    wgt_Q_rc = (/wgt_Q_rc*dim(0)/) 
return wgt_Q_rc
end 

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

function vertical_avg(data)
begin
    psfc =1013.0
    plev = data&level
    
    pver=pver*100
    pver@units="Pa"

    dp = pver              ; caculate dp
    dp(0:dimsizes(dp)-2) = abs(pver(1:) - pver(0:dimsizes(pver)-2))
    dp(dimsizes(dp)-1) = dp(dimsizes(dp)-2)

    dp_4d = conform_dims(dimsizes(data), dp, 1)
    
    q1dp = data(:,1,:,:)
    q1dp =  dim_sum_n_Wrap(data*dp_4d, 1)/dim_sum_n(dp_4d, 1)  ;;;积分
return q1dp
end 
;;;;;;;;;;;;
begin
 startmon = 197901
 endmon = 201512

  ;;;read era data
  path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
  var_name = (/"uwnd","vwnd","shum/shum.total"/)
  filein = path_wind + var_name +".mon.mean.nc"

  
  ;;;region 
    lat1 = -10
    lat2 = 30
    lon1 = 90
    lon2 = 180
  ;;;


    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; timeind :=ind(date.ge.197901.and.date.le.201512)

    yearnum = 37

    u = short2flt(f1->u(timeind,{1000:0},{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(1),"r")
    v = short2flt(f1->v(timeind,{1000:0},{lat1:lat2},{lon1:lon2}))

    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    timeind :=ind(date.ge.197901.and.date.le.201512)
    q = short2flt(f3->q(timeind,{1000:0},{lat1:lat2},{lon1:lon2}))        
    ;;;
print("=====finish read data==========")
;;;;


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

    dp = pver              ; caculate dp
    dp(0:dimsizes(dp)-2) = abs(pver(1:) - pver(0:dimsizes(pver)-2))
    dp(dimsizes(dp)-1) = dp(dimsizes(dp)-2)

    
    dp_cfm = conform_dims(dimsizes(era_qu),dp,1)          ; confirm dp to caculate sum(var*dp)
    copy_VarCoords(era_qu, dp_cfm)
    
    g_lon   =   era_qu&longitude
    g_lat   =   era_qu&latitude
    g_nlat  =   dimsizes(g_lat)

    Divqv = uv2dv_cfd(era_qu,era_qv,g_lat,g_lon,2)
    ; era_wgt_qu = (/dim_sum_n(dp_cfm*era_qu,1)/)               ; sum(var*dp)
    ; era_wgt_qv = (/dim_sum_n(dp_cfm*era_qv,1)/)               ; sum(var*dp)

    ; era_wgt_qu = (/era_wgt_qu/9.8/)
    ; era_wgt_qv = (/era_wgt_qv/9.8/)
    

    ;;;垂直积分
    DivQ = (/dim_sum_n(dp_cfm*Divqv,1)/)  
    DivQ   =   (/-DivQ*86400/9.8/) 
    copy_VarCoords(era_qu(:,1,:,:), DivQ)
    printVarSummary(DivQ)

    DivQ_annual = month_to_season(DivQ, "MAM")
    DivQ_delta = dim_rmvmean_n_Wrap(DivQ_annual,0)
    wgt_Q_rc = rc(DivQ_annual)
    wgt_Q_delta_rc = rc(DivQ_delta)
 

    prob_Q = rc_test(wgt_Q_rc)

    
    print("=====finish DivQ calculate==========")

    ;;;;;

    ;;;;计算动力项和热力项

     q_bar = month_to_season(q, "MAM")
     u_bar = month_to_season(u, "MAM")
     v_bar = month_to_season(v, "MAM") 

     printVarSummary(q_bar)
     
     dq = grad_latlon_cfd(q_bar, g_lat, g_lon, False, False)

     era_adv_q_bar = q_bar                         
     era_adv_q_bar = u_bar*dq[1] + v_bar*dq[0]    ;;dq[o]:dq/dy，dq[1]:dq/dx

     dp_cfm_new = conform_dims(dimsizes(u_bar),dp,1)          ; confirm dp to caculate sum(var*dp)
     copy_VarCoords(u_bar, dp_cfm_new)
     era_wgt_adv_q_bar=(/dim_sum_n(dp_cfm_new*era_adv_q_bar,1)/)               ; sum(var*dp)
     era_wgt_adv_q_bar = (/era_wgt_adv_q_bar*86400/(-9.8)/)       ;;convert to day
     copy_VarCoords(q_bar(:,1,:,:), era_wgt_adv_q_bar)

  ; Calculate the MFC_convergence term
  ;   MFC_conv   = -q*((du/dx)+(dv/dy) )      ; con(div)-vergence
    duv = u_bar
    duv  = uv2dv_cfd(u_bar,v_bar,g_lat,g_lon,2)       ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)
    
    ;;; all level 
    era_con_q_bar  = q_bar*duv
    era_wgt_con_q_bar=(/dim_sum_n(dp_cfm_new*era_con_q_bar,1)/)               ; sum(var*dp)
    era_wgt_con_q_bar = (/era_wgt_con_q_bar*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar)
   ;;;

    ;;;;1000-850
    era_con_q_bar_0085  = q_bar(:,{1000:850},:,:)*duv(:,{1000:850},:,:)
    era_wgt_con_q_bar_0085=(/dim_sum_n(dp_cfm_new(:,{1000:850},:,:)*era_con_q_bar_0085,1)/)               ; sum(var*dp)
    era_wgt_con_q_bar_0085 = (/era_wgt_con_q_bar_0085*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar_0085)
    ;;;


    ;;;850-100积分
    era_con_q_bar_8500  = q_bar(:,{825:100},:,:)*duv(:,{825:100},:,:)
    era_wgt_con_q_bar_8500=(/dim_sum_n(dp_cfm_new(:,{825:100},:,:)*era_con_q_bar_8500,1)/)               ; sum(var*dp)
    era_wgt_con_q_bar_8500 = (/era_wgt_con_q_bar_8500*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar_8500)
    ;;;

    ;;;求变化
    dyn_Q_trend  = rc(era_wgt_con_q_bar)
    prob_dyn = rc_test(dyn_Q_trend)

    dyn_Q_trend_0085  = rc(era_wgt_con_q_bar_0085)
    prob_dyn_0085 = rc_test(dyn_Q_trend_0085)

    dyn_Q_trend_8500  = rc(era_wgt_con_q_bar_8500)
    prob_dyn_8500 = rc_test(dyn_Q_trend_8500)
    
    ;;;;
   


    ;;;;平流项
    therdyn_Q_trend  = regCoef_n(ispan(1, yearnum, 1), era_wgt_adv_q_bar, 0, 0)
    copy_VarCoords(era_wgt_adv_q_bar(1,:,:), therdyn_Q_trend)
    therdyn_Q_trend = (/therdyn_Q_trend*yearnum/)

    prob_therdyn = rc_test(therdyn_Q_trend)
    
    therdyn_Q_trend = wgt_Q_rc - dyn_Q_trend

    delete(q)
    delete(u)
    delete(v)
;;;   
    wgt_Q_rc_areamean = avg(DivQ_annual(:,{0:20},{125:160}))

    ; print(wgt_Q_rc_areamean)
    ; exit()
;;;;;;;

;;;;;;;;;;;;ERA-interim datasets 

  f2=addfile("/home/yangsong3/data-observation/ERA-interim/monthly/precipitation/ERAI-pre-1degree.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201512)
    lsp=f2->lsp(timeind,:,:)  ;;;units are m
    cp=f2->cp(timeind,:,:)   
    
    ;;convert units to mm/day
    p_pi = cp
    p_pi = (/lsp+cp/)
    p_pi = p_pi*10^3
    
    p_pi_annual = month_to_season(p_pi, "MAM")

    rc_pre =  rc(p_pi_annual)
    prob_pre = rc_test(rc_pre)
;;;;;;



;;;print percent 

   percent1 =avg(dyn_Q_trend({0:20},{125:160}))/avg(wgt_Q_rc({0:20},{125:160}))
   percent2 =avg(therdyn_Q_trend({0:20},{125:160}))/avg(wgt_Q_rc({0:20},{125:160}))

   percent3 =avg(dyn_Q_trend_0085({0:20},{125:160}))/avg(dyn_Q_trend({0:20},{125:160}))
   percent4 =avg(dyn_Q_trend_8500({0:20},{125:160}))/avg(dyn_Q_trend({0:20},{125:160}))

   ; print(percent3)
   ; print(percent4)
   print(avg(wgt_Q_rc({0:20},{125:160})))
   print(avg(dyn_Q_trend({0:20},{125:160})))



    wks = gsn_open_wks("pdf", "Fig.R5.CMVIF&ERA-interim PRE trend")

    plot = new(6,graphic)
    iplot = plot
;;;;;
    res4=True
    res4@gsnFrame=False
    res4@gsnDraw=False
    res4@gsnAddCyclic = False
  
    res4@cnFillOn=True
    res4@cnLinesOn=False

    res4@cnLineLabelsOn=False
    res4@lbLabelBarOn=False
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

    
    ; res4@cnFillPalette = "MPL_PuOr"
    ; colordata = read_colormap_file("/home/yangsong3/data-observation/linshh/data/colordata/bluewhiteorangered_11color.rgb")
    ; res4@cnFillPalette=colordata(::-1,:)
    
    
    ; res4@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res4@cnMinLevelValF       =   -3.0
    ; res4@cnMaxLevelValF       =   3.0
    ; res4@cnLevelSpacingF      =  1.0
    
    res4@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res4@cnLevels = (/-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    res4@cnFillPalette = "precip_diff_12lev"
    ; res4@cnFillPalette = "CBR_coldhot"
    res4@cnFillColors= (/2,4,6,8,10,11,12,13,14/) - 2
      
    res4@gsnLeftString = "(a) CVIMF trend "
    res4@gsnLeftStringFontHeightF = 0.03
    ; plot = gsn_csm_contour_map(wks, wgt_w_rc, res4)
    plot(0) = gsn_csm_contour_map(wks, wgt_Q_rc, res4)
    res4@gsnLeftString = "(c) Moisture convergence (MC) term trend "
    plot(2) = gsn_csm_contour_map(wks, dyn_Q_trend, res4)
    res4@gsnLeftString = "(d) Moisture advection (MA) term trend"
    plot(3) = gsn_csm_contour_map(wks, therdyn_Q_trend, res4)

   
    res5 = res4

    ; res5@cnFillPalette = "precip_diff_12lev"
    ; res5@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ; res5@cnLevels = (/-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    ; res5@cnFillColors= (/5,6,7,8,10,11,12,13,14/) - 2

    res5@gsnLeftString = "(b) ERA-interim precipitation trend"
    plot(1) = gsn_csm_contour_map(wks, rc_pre, res5)
    

    res4@gsnLeftString = "(e) MC term 1000-850hPa "
    plot(4) = gsn_csm_contour_map(wks, dyn_Q_trend_0085, res4)
 
    res4@gsnLeftString = "(f) MC term above 850hPa"
    plot(5) = gsn_csm_contour_map(wks, dyn_Q_trend_8500, res4)



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
    str = ispan(1, 6, 1) + ""

    do i = 0,5
      plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
     end do 



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



    iplot(0) = gsn_csm_contour(wks,prob_Q,res_pro4)

    iplot(1) = gsn_csm_contour(wks,prob_pre,res_pro4)

    iplot(2) = gsn_csm_contour(wks,prob_dyn,res_pro4)

    iplot(3) = gsn_csm_contour(wks,prob_therdyn,res_pro4)

    iplot(4) = gsn_csm_contour(wks,prob_dyn_0085,res_pro4)

    iplot(5) = gsn_csm_contour(wks,prob_dyn_8500,res_pro4)

   
    do i=0,5
    overlay(plot(i), iplot(i))
    end do 
  ;;;;

pres =True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.01
gsn_panel(wks, plot, (/3,2/), pres)
end



