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
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 
begin
   startmon = 198001
   endmon = 202012
   ;;read  eva data 
    year = ispan(1980,2020,1) 
   ;;;;
    ; f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    ; f2=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")
    ; date:=cd_calendar(f2->time, 1)
    ; timeind:=ind(date.ge.startmon.and.date.le.endmon)
    ; p=f2->precip(timeind,:,:)

    ; p = (/rmMonAnnCycTLL(p)/)

    f2=addfile("/home/yangsong3/data-observation/NCEP_NCAR/monthly/prate.sfc.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p_r=f2->prate(timeind,{-90:90},{0:360})
    p_r = (/p_r*86400/)
    p_r = (/rmMonAnnCycTLL(p_r)/)

    ;;;处理数据
    ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

    yo = fspan(-88.75, 88.75, 72)
    xo = fspan(1.25, 358.75, 144)
    p = linint2_Wrap(p_r&lon, p_r&lat, p_r, True, xo, yo, 0)

    delete(p&LAT)
    delete(p&LON)
    p!1 = "lat"
    p!2 = "lon"
    p&lat = yo
    p&lon = xo

   
;;;;
    f3=addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/evapr_oaflux_chazhi.nc", "r")
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    eva=f3->evapr(timeind,:,:)
    eva = (/eva*10/365/) ;;convert to mm/day
    eva = (/rmMonAnnCycTLL(eva)/)
;;;;


    ;;;read era data

    f4 = addfile("/home/yangsong3/data-observation/NCEP_NCAR/monthly/shum.mon.mean.nc","r")
    date :=cd_calendar(f4->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    q = short2flt(f4->shum(timeind,{1000:300},:,:))
    q = q/1000.0
    ;;density of water is 10^3 kg/m3, intergrated qv units: mm/s
  
    pver=q&level
    pver=pver*100.0
    pver@units="Pa"

    dp=dpres_plevel(pver,100000.0,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(q),dp,1)          ; confirm dp to caculate sum(var*dp)
    q_col = (/dim_sum_n(dp_cfm*q,1)/)               ; sum(var*dp)
    copy_VarCoords(q(:,1,:,:), q_col)

    q_col = q_col/(9.8)  ;;; divide pg, convert to mm

    q_col = (/center_finite_diff_n(q_col, 1, False, 0, 0)/)  ;;mm/month

    q_col = q_col/30.0 ;;;convert to mm/day
    
;;;;    



;;;时间序列
    evp_annual = month_to_annual(eva, 1)
    p_annual  = month_to_annual(p, 1)
    q_annual  = month_to_annual(q_col, 1)
    
    evp_12 = split_12(eva)
    p_12 = split_12(p)
    q_12 = split_12(q_col)


    evp_annual = (/dim_avg_n_Wrap(evp_12(0:4,:,:,:), 0)/)
    p_annual = (/dim_avg_n_Wrap(p_12(0:4,:,:,:), 0)/)
    q_annual = (/dim_avg_n_Wrap(q_12(0:4,:,:,:), 0)/)


    p_e = p_annual
    p_e = (/p_annual - evp_annual/)
    copy_VarCoords(p_annual, p_e)
    
    ;;;
    lat1 = 0
    lat2 = 20
    lon1 = 250
    lon2 = 280

    p_e_t = dim_avg_n_Wrap(p_e(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    p_annual_t = dim_avg_n_Wrap(p_annual(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    evp_annual_t = dim_avg_n_Wrap(evp_annual(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    q_annual_t = dim_avg_n_Wrap(q_annual(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
;;;

 ;=================================================================
; ; Mission 2.  trend
; ; ==============================================================
;;;读入Q-term
  f3 = addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/JFMAM-Q-term-8-NEP-2020.nc", "r")
  series_Q = f3->series_Q 
  series_dyna_adv_q = f3->series_dyna_adv_q
  series_dyna_con_q = f3->series_dyna_con_q
  series_thermaldyna_adv_q = f3->series_thermaldyna_adv_q
  series_thermaldyna_con_q = f3->series_thermaldyna_con_q
  series_quad_adv_q = f3->series_quad_adv_q
  series_quad_con_q = f3->series_quad_con_q
  series_eddy_q = f3->series_eddy_q

;;;
startyear = 1980
endyear  = 2020
;;;;
yearnum  = endyear - startyear +1
;==============================================================
; ; Mission 4.  trend
; ; ==============================================================

  rc = regCoef_n(ispan(1,yearnum,1), q_annual_t, 0, 0)
  rc = rc*yearnum


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


    p_e_rc = regline(ispan(1,yearnum,1), p_e_t)
    p_e_rc = (/p_e_rc*yearnum/)

    e_rc = regline(ispan(1,yearnum,1), evp_annual_t)
    e_rc = (/e_rc*yearnum/)

    p_rc = regline(ispan(1,yearnum,1), p_annual_t)
    p_rc = (/p_rc*yearnum/)


    q_rc = regline(ispan(1,yearnum,1), q_annual_t)
    q_rc = (/q_rc*yearnum/)

    p_e_rc = (/p_rc - e_rc/)

;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "P-E以及Q的各项的时间序列-2020-NEP")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = endyear
  resxy@trXMinF = startyear
  resxy@tmXBValues = ispan(startyear,endyear,5)
  resxy@tmXBLabels = ispan(startyear,endyear,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(startyear,endyear,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknesses =  (/2,2,2,1,1,1,1,1,1,1/)          ; make second line thicker

  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0,0,1,1,1,1,1,1,1/)
  
  resxy@xyLineColors       =  (/"black","blue","red", "green","orange","pink","purple","yellow","gold4","gray"/)          ; change line color
  resxy@xyMarkLineModes   =  (/"lines"/)
  resxy@xyMonoMarkLineMode = True
  
  
  time = ispan(startyear, endyear, 1)


  plot = gsn_csm_xy(wks,time,(/q_annual_t, p_e_t,series_Q,series_dyna_adv_q,\
                              series_dyna_con_q,series_thermaldyna_adv_q,series_thermaldyna_con_q,series_quad_adv_q,series_quad_con_q,\
                          series_eddy_q/), resxy)

  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       = resxy@xyLineColors 
  lineres@lgLineThicknesses =    resxy@xyLineThicknesses                 ; line thicknesses
  lineres@LineLengthPercent = 7.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = resxy@xyDashPatterns


  textres = True
  textres@lgLabels = (/"TPW","P-E","DIVQ","DY(ADV)","DY(CON)","TH(ADV)","TH(CON)", "QD(ADV)","QD(CON)","EDDY"/) 
  textres@lgLabelFontHeightF = 0.005


  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)
  

;==============================================================
; Mission 5. bar plot
; ==============================================================
  
  x = (/1,2,3,4,5,6,7,8,9/)
  y = (/p_e_rc,q_delta_rc,q_dyna_adv_rc,q_dyna_con_rc,q_thermaldyna_adv_rc,\
       q_thermaldyna_con_rc,q_quad_adv_rc,q_quad_con_rc,q_eddy_rc/)

  labels = (/"P-E","DIVQ","DY(ADV)","DY(CON)","TH(ADV)","TH(CON)",\
  "QD(ADV)","QD(CON)","EDDY"/)
;********************************************************
; create the plot
;********************************************************   
  wks1 = gsn_open_wks("pdf","bar-q-component-trend")                  ; send graphics to PNG file
  res                       = True                 ; plot mods desired         
  res@gsnFrame              = False                ; don't advance frame yet
  res@gsnXYBarChart         = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  ; res@gsnXYBarChartColors = (/"firebrick","red","orange","green", \
  ; "navy","blue","SkyBlue","SlateBlue"/)            ; choose colors
  
  res@tmXBOn                = True                ; turn off tickmarks at bot
  res@trYMinF               = -0.3                   ; bring bars down to zero
  res@trYMaxF               = 1.5                   ; bring bars down to zero

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

  plot1 = gsn_csm_xy (wks1,x,y,res)                  ; create plot
  
 draw(plot1)
 frame(wks1)

  end