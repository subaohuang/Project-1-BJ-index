;==============================================================
; Mission read  wind data
; ==============================================================
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

    u = short2flt(f1->u(timeind,{1000:100},{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(1),"r")
    v = short2flt(f1->v(timeind,{1000:100},{lat1:lat2},{lon1:lon2}))

    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    timeind :=ind(date.ge.197901.and.date.le.201512)
    q = short2flt(f3->q(timeind,{1000:100},{lat1:lat2},{lon1:lon2}))        
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

    era_wgt_qu = (/dim_sum_n(dp_cfm*era_qu,1)/)               ; sum(var*dp)
    era_wgt_qv = (/dim_sum_n(dp_cfm*era_qv,1)/)               ; sum(var*dp)

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

    DivQ_annual = month_to_season(DivQ,"MAM")
    ; yearnum = 37
    wgt_Q_rc := regCoef_n(ispan(1,yearnum,1), DivQ_annual, 0, 0)
    copy_VarCoords(DivQ_annual(1,:,:), wgt_Q_rc)

    prob_Q = rc_test(wgt_Q_rc)
    wgt_Q_rc = (/wgt_Q_rc*yearnum/) 


  
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
    
    dyn_Q_trend  = regCoef_n(ispan(1, yearnum, 1), era_wgt_con_q_bar, 0, 0)
    copy_VarCoords(era_wgt_con_q_bar(1,:,:), dyn_Q_trend)
    dyn_Q_trend = (/dyn_Q_trend*yearnum/)

    therdyn_Q_trend  = regCoef_n(ispan(1, yearnum, 1), era_wgt_adv_q_bar, 0, 0)
    copy_VarCoords(era_wgt_adv_q_bar(1,:,:), therdyn_Q_trend)
    therdyn_Q_trend = (/therdyn_Q_trend*yearnum/)


    prob_dyn = rc_test(dyn_Q_trend)
    prob_therdyn = rc_test(therdyn_Q_trend)
    
    ; therdyn_Q_trend = wgt_Q_rc - dyn_Q_trend -therdyn_Q_trend

    delete(q)
    delete(u)
    delete(v)
;;;
 
    
 wgt_Q_rc_areamean = avg(DivQ_annual(:,{0:20},{125:160}))

 ; print(wgt_Q_rc_areamean)
 ; exit()


wks = gsn_open_wks("pdf", "CVIMF-trend-1000-100")

plot = new(3,graphic)
iplot = plot
;;;;;
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

    
    res4@cnFillPalette = "MPL_PuOr"
    
    
    res4@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res4@cnMinLevelValF       =   -2.0
    res4@cnMaxLevelValF       =   2.0
    res4@cnLevelSpacingF      =  0.5

    
    res4@gsnLeftString = "CVIMF"
    res4@gsnLeftStringFontHeightF = 0.03
    ; plot = gsn_csm_contour_map(wks, wgt_w_rc, res4)
    plot(0) = gsn_csm_contour_map(wks, wgt_Q_rc, res4)
    res4@gsnLeftString = "DYNAMIC"
    plot(1) = gsn_csm_contour_map(wks, dyn_Q_trend, res4)
    res4@gsnLeftString = "Thermaldynamic"
    plot(2) = gsn_csm_contour_map(wks, therdyn_Q_trend, res4)

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
    str = ispan(1, 3, 1) + ""

    do i = 0,2 
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
    iplot(1) = gsn_csm_contour(wks,prob_dyn,res_pro4)
    iplot(2) = gsn_csm_contour(wks,prob_therdyn,res_pro4)


    do i=0,2
    overlay(plot(i), iplot(i))
    end do 
  ;;;;

pres =True

gsn_panel(wks, plot, (/3,1/), pres)
end



