;;********************************************************
;;计算蒸发的趋势以及查看P-E的大小
;;oaflux的蒸发数据，文件的单位为cm/day，短型
;;以及GCCP的降水数据，
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
;****************************************
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
   ;;read  eva data 
    years = ispan(1979, 2016, 1)
    fils =  "/home/yangsong3/data-observation/OAflux/monthly/evaporation/evapr_oaflux_"+years+".nc"

    f1 = addfiles(fils, "r")
    date = cd_calendar(f1[:]->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    evp_r =short2flt(f1[:]->evapr(timeind,:,:))
    evp_r@missing_value = -9.96921e+36
    evp_r@_FillValue =   -9.96921e+36
    evp_r&lat@units = "degrees_north"
    evp_r&lon@units = "degrees_east"
    
    evp_r = where(evp_r.eq.3276.6, -9.96921e+36, evp_r)
    evp =  evp_r*10 / 365.0   ;;convert t0 day
    copy_VarMeta(evp_r, evp)
    delete(evp_r)
   
    printVarSummary(evp)
    
    
    ;;;处理数据
        ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

        yo = fspan(-88.75, 88.75, 72)
        xo = fspan(1.25, 358.75, 144)

        evp_new = linint2_Wrap(evp&lon, evp&lat, evp, True, xo, yo, 0)

        delete(evp_new&LAT)
        delete(evp_new&LON)
        evp_new!1 = "lat"
        evp_new!2 = "lon"
        evp_new&lat = yo
        evp_new&lon = xo
        printVarSummary(evp_new)
        ;;;

        ; ;;116.25存在插值后缺侧的现象，现在插值下
        ; m =  linmsg(u_av_1({100:120}),0)
        ; copy_VarCoords(u_av_1({100:120}), m)
        ; u_av_1({116.25}) = m({116.25})
        ; delete(m)
        ; ;;;
;;;;

;;;read pre data

    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.201612)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)  
;;;

;;;时间序列
    evp_new = (/rmMonAnnCycTLL(evp_new)/)
    p = (/rmMonAnnCycTLL(p)/)

    evp_annual = month_to_annual(evp_new, 1)
    p_annual  = month_to_annual(p, 1)
    p_e = p_annual
    p_e = (/p_annual - evp_annual/)
    copy_VarCoords(p_annual, p_e)

    dim = dimsizes(p_e)
    yearnum = dim(0)
    p_e_t = dim_avg_n_Wrap(p_e(:,{0:15},{125:160}), (/1,2/))
    p_annual_t = dim_avg_n_Wrap(p_annual(:,{0:15},{125:150}), (/1,2/))
    evp_annual_t = dim_avg_n_Wrap(evp_annual(:,{0:15},{125:160}), (/1,2/))
;;;

;;;读入Q-term
    f3 = addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/Q-term-8-philipian.nc", "r")
    series_Q = f3->series_Q
    
    series_dyna_adv_q  = f3->series_dyna_adv_q 
    series_dyna_con_q  = f3->series_dyna_con_q 
    series_thermaldyna_adv_q = f3->series_thermaldyna_adv_q
    series_thermaldyna_con_q = f3->series_thermaldyna_con_q
 
    series_quad_adv_q  = f3->series_quad_adv_q 
    series_quad_con_q  = f3->series_quad_con_q 

    series_eddy_q = f3->series_eddy_q
;;;
 
;;;; 计算趋势trend 
    p_e_rc = regline(ispan(1,yearnum,1), p_e_t)
    p_e_rc = (/p_e_rc*yearnum/)

    q_delta_rc = regline(ispan(1,yearnum,1), series_Q)
    q_delta_rc = (/q_delta_rc*yearnum/)

    print("q_delta_rc="+q_delta_rc)

    q_dyna_adv_rc = regline(ispan(1,yearnum,1), series_dyna_adv_q)
    q_dyna_adv_rc = (/q_dyna_adv_rc*yearnum/)
    
    print("q_dyna_adv_rc="+q_dyna_adv_rc)

    q_dyna_con_rc = regline(ispan(1,yearnum,1), series_dyna_con_q)
    q_dyna_con_rc = (/q_dyna_con_rc*yearnum/)
    
    print("q_dyna_con_rc="+q_dyna_con_rc)

    q_thermaldyna_adv_rc = regline(ispan(1,yearnum,1), series_thermaldyna_adv_q)
    q_thermaldyna_adv_rc = (/q_thermaldyna_adv_rc*yearnum/)
    
    q_thermaldyna_con_rc = regline(ispan(1,yearnum,1), series_thermaldyna_con_q)
    q_thermaldyna_con_rc = (/q_thermaldyna_con_rc*yearnum/)

   
    q_quad_adv_rc = regline(ispan(1,yearnum,1), series_quad_adv_q)
    q_quad_adv_rc = (/q_quad_adv_rc*yearnum/)
    
    q_quad_con_rc = regline(ispan(1,yearnum,1), series_quad_con_q)
    q_quad_con_rc = (/q_quad_con_rc*yearnum/)

    
    q_eddy_rc = regline(ispan(1,yearnum,1), series_eddy_q)
    q_eddy_rc = (/q_eddy_rc*yearnum/)
;;;;



;;;
;==============================================================
; Mission 5. bar plot
; ==============================================================
  
  x = (/1,2,3,4,5,6,7,8,9/)
  y = (/p_e_rc, q_delta_rc,q_dyna_adv_rc,q_dyna_con_rc,q_thermaldyna_adv_rc,\
       q_thermaldyna_con_rc,q_quad_adv_rc,q_quad_con_rc,q_eddy_rc/)

  labels = (/"P-E","DIVQ","ADV(DY)","CON(DY)","ADV(TH)","CON(TH)",\
  "ADV(QD)","CON(QD)","EDDY"/)
;********************************************************
; create the plot
;********************************************************   
  wks1 = gsn_open_wks("pdf","bar-q-component-trend-2016")                  ; send graphics to PNG file
  res                       = True                 ; plot mods desired         
  res@gsnFrame              = False                ; don't advance frame yet
  res@gsnXYBarChart         = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  ; res@gsnXYBarChartColors = (/"firebrick","red","orange","green", \
  ; "navy","blue","SkyBlue","SlateBlue"/)            ; choose colors
  
  res@tmXBOn                = True                ; turn off tickmarks at bot
  res@trYMinF               = -0.4                    ; bring bars down to zero
  res@trYMaxF               = 1.0                   ; bring bars down to zero

  res@trXMinF               = 0                    ; adds space on either end
  res@trXMaxF               = 10                    ; of the 1st and last bars

  res@gsnYRefLine = 0.      ; set the reference line equal to 0
  res@gsnAboveYRefLineColor = "Gray" 
  res@gsnBelowYRefLineColor = "white"
  
  ; res@tiMainString          = "Estimated Cancer Deaths for 2002"
;---Put city labels on Y axis
  res@tmXBMode                 = "Explicit"
  res@tmXBValues               = x
  res@tmXBLabels               = labels
  res@tmXBLabelFontHeightF = 0.01

  plot = gsn_csm_xy (wks1,x,y,res)                  ; create plot
  
 draw(plot)
 frame(wks1)



  end