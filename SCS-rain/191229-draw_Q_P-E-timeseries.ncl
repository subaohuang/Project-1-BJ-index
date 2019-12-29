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
    fils = systemfunc("ls /home/yangsong3/data-observation/OAflux/monthly/evaporation/evapr_oaflux_"+"*"+".nc")

    f1 = addfiles(fils, "r")
    date = cd_calendar(f1[:]->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
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
    timeind:=ind(date.ge.197901.and.date.le.201312)
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

    p_e_t = dim_avg_n_Wrap(p_e(:,{0:15},{110:155}), (/1,2/))
    p_annual_t = dim_avg_n_Wrap(p_annual(:,{0:15},{110:155}), (/1,2/))
    evp_annual_t = dim_avg_n_Wrap(evp_annual(:,{0:15},{110:155}), (/1,2/))
;;;

;;;读入Q-term
    f3 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Q-term-5.nc", "r")
    series_Q = f3->series_Q(:)
;;;

;;;

;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "Q以及P-E的时间序列")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 2013
  resxy@trXMinF = 1979
  resxy@tmXBValues = ispan(1979,2013,5)
  resxy@tmXBLabels = ispan(1979,2013,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(1979,2013,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
  resxy@xyMonoLineThickness = True
  resxy@xyLineColors       =  (/"black","blue","red","green","gray"/)          ; change line color
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines","lines"/)
  resxy@xyExplicitLegendLabels  =  (/"main","dyn","thermaldynamic","quandratic", "eddy"/)     ;加线的标注
  

  time = ispan(1979, 2013, 1)
  plot = gsn_csm_xy(wks,time,(/p_e_t,series_Q,p_annual_t,evp_annual_t/), resxy)

  draw(plot)
  frame(wks)
  end