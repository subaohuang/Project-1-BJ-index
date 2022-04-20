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
    lat1 = -10.0
    lat2 = 30.0
    lon1 = 100.0
    lon2 = 180.0
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

    ;;;垂直积分
    DivQ = (/dim_sum_n(dp_cfm*Divqv,1)/)  
    DivQ   =   (/-DivQ*86400/9.8/) 
    copy_VarCoords(era_qu(:,1,:,:), DivQ)
    printVarSummary(DivQ)

    DivQ_annual = month_to_season(DivQ, "MAM")
    wgt_Q_rc = rc(DivQ_annual)

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

;   ; Calculate the MFC_convergence term
;   ;   MFC_conv   = -q*((du/dx)+(dv/dy) )      ; con(div)-vergence
;     duv = u_bar
;     duv  = uv2dv_cfd(u_bar,v_bar,g_lat,g_lon,2)       ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)
    
;     ;;; all level 
;     era_con_q_bar  = q_bar*duv
;     era_wgt_con_q_bar=(/dim_sum_n(dp_cfm_new*era_con_q_bar,1)/)               ; sum(var*dp)
;     era_wgt_con_q_bar = (/era_wgt_con_q_bar*86400/(-9.8)/)       ;;convert to day
;     copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar)
;    ;;;

;     ;;;;1000-850
;     era_con_q_bar_0085  = q_bar(:,{850:},:,:)*duv(:,{850:},:,:)
;     era_wgt_con_q_bar_0085=(/dim_sum_n(dp_cfm_new(:,{850:},:,:)*era_con_q_bar_0085,1)/)               ; sum(var*dp)
;     era_wgt_con_q_bar_0085 = (/era_wgt_con_q_bar_0085*86400/(-9.8)/)       ;;convert to day
;     copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar_0085)
;     ;;;


;     ;;;850-100积分
;     era_con_q_bar_8500  = q_bar(:,{850:100},:,:)*duv(:,{850:100},:,:)
;     era_wgt_con_q_bar_8500=(/dim_sum_n(dp_cfm_new(:,{850:100},:,:)*era_con_q_bar_8500,1)/)               ; sum(var*dp)
;     era_wgt_con_q_bar_8500 = (/era_wgt_con_q_bar_8500*86400/(-9.8)/)       ;;convert to day
;     copy_VarCoords(q_bar(:,1,:,:), era_wgt_con_q_bar_8500)
;     ;;;

;     ;;求变化
;     dyn_Q_trend  = rc(era_wgt_con_q_bar)
;     prob_dyn = rc_test(dyn_Q_trend)

;     dyn_Q_trend_0085  = rc(era_wgt_con_q_bar_0085)
;     prob_dyn_0085 = rc_test(dyn_Q_trend_0085)

;     dyn_Q_trend_8500  = rc(era_wgt_con_q_bar_8500)
;     prob_dyn_8500 = rc_test(dyn_Q_trend_8500)
    
;     ;;;
   


;     ;;;;平流项
;     therdyn_Q_trend  = regCoef_n(ispan(1, yearnum, 1), era_wgt_adv_q_bar, 0, 0)
;     copy_VarCoords(era_wgt_adv_q_bar(1,:,:), therdyn_Q_trend)
;     therdyn_Q_trend = (/therdyn_Q_trend*yearnum/)

;     prob_therdyn = rc_test(therdyn_Q_trend)
    
;     therdyn_Q_trend = wgt_Q_rc - dyn_Q_trend

;     delete(q)
;     delete(u)
;     delete(v)
; ;;;   
;     wgt_Q_rc_areamean = avg(DivQ_annual(:,{0:20},{125:160}))

;     ; print(wgt_Q_rc_areamean)
;     ; exit()
;;;;;;;

  ;;;;change δq ̅∇∙V ̅_clm+q ̅_clm ∇∙δV ̅+δV ̅∙∇q ̅_clm+V ̅_clm∙∇δq ̅
     q_clm = dim_avg_n_Wrap(q_bar, 0)
     u_clm = dim_avg_n_Wrap(u_bar, 0)
     v_clm = dim_avg_n_Wrap(v_bar, 0)

     q_change =  regCoef_n(ispan(1, yearnum, 1), q_bar, 0, 0)
     copy_VarCoords(q_bar(1,:,:,:), q_change)
     q_change = (/q_change*yearnum/)

     u_change =  regCoef_n(ispan(1, yearnum, 1), u_bar, 0, 0)
     copy_VarCoords(q_bar(1,:,:,:), u_change)
     u_change = (/u_change*yearnum/)

     v_change =  regCoef_n(ispan(1, yearnum, 1), v_bar, 0, 0)
     copy_VarCoords(q_bar(1,:,:,:), v_change)
     v_change = (/v_change*yearnum/)
  ;;;
    
  
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


  ;;;;

   
  ;;;qclm▽(σV)

    duv_change = u_change
    duv_change  = uv2dv_cfd(u_change,v_change,g_lat,g_lon,2)       ; (1/m)(m/s) => (1/s) ; (time,level,lat,lon)
    
    ;;; all level 
    era_con_u_change = q_clm
    era_con_u_change  = duv_change*q_clm

    era_wgt_con_u_change = (/dim_sum_n(dp_cfm_new2*era_con_u_change,0)/)               ; sum(var*dp)
    era_wgt_con_u_change = (/era_wgt_con_u_change*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_clm(1,:,:), era_wgt_con_u_change)
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
    
   
    printVarSummary(dp_cfm_new2)
   
;     ;;;;1000-850
;     era_con_u_change_0085  = era_con_u_change({1000:850},:,:)
;     era_wgt_con_u_change_0085=(/dim_sum_n(dp_cfm_new2({1000:850},:,:)*era_con_u_change_0085,0)/)               ; sum(var*dp)
;     era_wgt_con_u_change_0085 = (/era_wgt_con_u_change_0085*86400/(-9.8)/)       ;;convert to day
;     copy_VarCoords(q_clm(1,:,:), era_wgt_con_u_change_0085)
;     ;;;

;     ;;;850-100积分
;     era_con_u_change_8500  = era_con_u_change({825:0},:,:)
;     era_wgt_con_u_change_8500=(/dim_sum_n(dp_cfm_new2({825:0},:,:)*era_con_u_change_8500,0)/)               ; sum(var*dp)
;     era_wgt_con_u_change_8500 = (/era_wgt_con_u_change_8500*86400/(-9.8)/)       ;;convert to day
;     copy_VarCoords(q_clm(1,:,:), era_wgt_con_u_change_8500)
;     ;;;

;    percent1 =avg(era_wgt_con_u_change_0085({0:20},{125:160}))/avg(era_wgt_con_u_change({0:20},{125:160}))

;    percent2 = avg(era_wgt_con_u_change_8500({0:20},{125:160}))/avg(era_wgt_con_u_change({0:20},{125:160})) +percent1

; print("percent1="+percent1)
; print("percent2 = "+percent2)
; print(avg(era_wgt_con_u_change_0085({0:20},{125:160})))
; print(avg(era_wgt_con_u_change_8500({0:20},{125:160})))
; print(avg(era_wgt_con_u_change({0:20},{125:160})))

   ; print(avg(era_wgt_con_u_change_0085({0:20},{125:160})))

   ; print(avg(era_wgt_adv_u_change_8500({0:20},{125:160}))) 

;;;print percent 

;    percent1 =avg(dyn_Q_trend({0:20},{125:160}))/avg(wgt_Q_rc({0:20},{125:160}))
;    percent2 =avg(therdyn_Q_trend({0:20},{125:160}))/avg(wgt_Q_rc({0:20},{125:160}))

;    percent3 =avg(dyn_Q_trend_0085({0:20},{125:160}))/avg(dyn_Q_trend({0:20},{125:160}))
;    percent4 =avg(dyn_Q_trend_8500({0:20},{125:160}))/avg(dyn_Q_trend({0:20},{125:160}))

;    ; print(percent3)
;    ; print(percent4)
;    print(avg(wgt_Q_rc({0:20},{125:160})))
;    print(avg(dyn_Q_trend({0:20},{125:160})))
; exit()


lat1 = 0.0
lat2 = 20.0
lon1 = 125.0
lon2 = 160.0
era_adv_u_change_mean = avg(era_wgt_adv_u_change({lat1:lat2},{lon1:lon2}))
era_adv_q_change_mean = avg(era_wgt_adv_q_change({lat1:lat2},{lon1:lon2}))

era_con_u_change_mean = avg(era_wgt_con_u_change({lat1:lat2},{lon1:lon2}))
era_con_q_change_mean = avg(era_wgt_con_q_change({lat1:lat2},{lon1:lon2}))

wgt_Q_rc_areamean  = avg(wgt_Q_rc({lat1:lat2},{lon1:lon2}))


;==============================================================
; Mission 5. bar plot
; ==============================================================
  
  ; x = (/1,2,3,4,5,6,7,8/)

  ; y = (/q_delta_rc,q_dyna_adv_rc,q_dyna_con_rc,q_thermaldyna_adv_rc,\
  ;      q_thermaldyna_con_rc,q_quad_adv_rc,q_quad_con_rc,q_eddy_rc/)
  ; labels = (/"TOTAL","ADV(DY)","CON(DY)","ADV(TH)","CON(TH)",\
  ; "ADV(QD)","CON(QD)","EDDY"/)

  x = (/1,2,3,4,5/)
printVarSummary(wgt_Q_rc_areamean)
printVarSummary(era_wgt_con_u_change)
printVarSummary(era_wgt_adv_u_change)
printVarSummary(era_wgt_adv_q_change)
printVarSummary(era_wgt_con_q_change)

  y = (/wgt_Q_rc_areamean,era_con_u_change_mean,era_adv_u_change_mean,era_con_q_change_mean,era_adv_q_change_mean/)

   labels = (/"DIVQ","DY(CON)","DY(ADV)","TH(CON)","TH(ADV)"/)
;********************************************************
; create the plot
;********************************************************   
  wks1 = gsn_open_wks("pdf","bar-q-component-change-trendmethod")                  ; send graphics to PNG file
  res                       = True                 ; plot mods desired         
  res@gsnFrame              = False                ; don't advance frame yet
  res@gsnXYBarChart         = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  ; res@gsnXYBarChartColors = (/"firebrick","red","orange","green", \
  ; "navy","blue","SkyBlue","SlateBlue"/)            ; choose colors
  
  res@tmXBOn                = True                ; turn off tickmarks at bot
  res@tmXTOn                = False
  res@trYMinF               = -0.4                    ; bring bars down to zero
  res@trYMaxF               = 1.2                   ; bring bars down to zero

  res@trXMinF               = 0                    ; adds space on either end
  res@trXMaxF               = 7                    ; of the 1st and last bars

  res@gsnYRefLine = 0.      ; set the reference line equal to 0
  res@gsnAboveYRefLineColor = "Gray" 
  res@gsnBelowYRefLineColor = "white"
  
  ; res@tiMainString          = "Estimated Cancer Deaths for 2002"
;---Put city labels on Y axis
  res@tmXBMode                 = "Explicit"
  res@tmXBValues               = x
  res@tmXBLabels               = labels
  res@tmXBLabelFontHeightF = 0.01
  
  res@vpHeightF = 0.4
  res@vpWidthF = 0.6
  plot = gsn_csm_xy (wks1,x,y,res)                  ; create plot
  
 draw(plot)
 frame(wks1)

end



