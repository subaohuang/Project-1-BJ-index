;;********************************************************
;;draw individual diabatic heat term profile averaged over SCS and 
;;philippines sea
;;diabatic and convection data from CFSR
;;                               lnshuheng 2020.1.2 
;**********************************************************
function get_data(var1,month)
begin
    if month.ne."1" then
      var1_annual = month_to_season(var1, month)
      if season(i).eq."DJF"
        var1_annual(0:dimsizes(var1_annual&time)-2,:,:,:) = var1_annual(1:,:,:,:) 
        var1_annual(dimsizes(var1_annual&time)-1,:,:,:) = var1(dimsizes(var1&time)-1,:,:,:)
     end if 
    end if
    
    var1_annual  = month_to_annual(var1, 1)

    var1_annual!0="time" 
    dimz = dimsizes(var1_annual)
    yr_length = dimz(0)

    var1_annual_rc = regCoef_n(ispan(1,yr_length,1), var1_annual,0,0)
    copy_VarCoords(var1_annual(1,:,:,:), var1_annual_rc)
    var1_annual_rc = var1_annual_rc*yr_length

    var1_clim = dim_avg_n_Wrap(var1_annual,0)
    ;;;区域平均

    g_latS=0.
    g_latN=20.
    g_lonW=125.
    g_lonE=160.

    var1_clim_area = wgt_areaave_Wrap(var1_clim(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
    var1_annual_rc_area = wgt_areaave_Wrap(var1_annual_rc(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)

    y = var1_annual(:1,:,0,0)
    y(0,:)= (/var1_clim_area/)
    y(1,:)= y(0,:)+(/var1_annual_rc_area/)
  return y
end 
begin
; ==============================================================
; Mission 1. read data
; ==============================================================
;;read diabatic data
  f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/review/new/individual_heating_CFSR_ANN_clm_change_profile-WNP.nc", "r")

  lschr_heat_change = f1->lschr_q(1,:)   ;;;large scale condensate heating
  lschr_heat_clm = f1->lschr_q(0,:)   ;;;

  level = lschr_heat_clm&lev
  level = level/100.0

  dchr_heat_change = f1->dchr_q(1,:)   ;;;Deep convective heating rate 
  dchr_heat_clm = f1->dchr_q(0,:)   ;;;

  
  srhr_heat_change = f1->srhr_q(1,:)   ;;;solar radiation heating rate 
  srhr_heat_clm = f1->srhr_q(0,:)   ;;;
  
  
  lwrhr_heat_change = f1->lwrhr_q(1,:)   ;;;long-wave radiative heating rate 
  lwrhr_heat_clm = f1->lwrhr_q(0,:)   ;;;


  lschr_heat_clm = lschr_heat_clm*86400.0
  lschr_heat_change = lschr_heat_change*86400.0

  dchr_heat_clm = dchr_heat_clm*86400.0
  dchr_heat_change = dchr_heat_change*86400.0

  srhr_heat_clm = srhr_heat_clm*86400.0
  srhr_heat_change = srhr_heat_change*86400.0

  lwrhr_heat_clm = lwrhr_heat_clm*86400.0
  lwrhr_heat_change = lwrhr_heat_change*86400.0
;;;

;;;读入浅对流和扩散
  schr_heat_change = f1->schr_q(1,:)   ;;;shallow convection rate
  schr_heat_clm = f1->schr_q(0,:)   ;;;

  vdhr_heat_change = f1->vdhr_q(1,:)   ;;;l
  vdhr_heat_clm = f1->vdhr_q(0,:)   ;;;

  schr_heat_clm = schr_heat_clm*86400.0
  schr_heat_change = schr_heat_change*86400.0

  vdhr_heat_clm = vdhr_heat_clm*86400.0
  vdhr_heat_change = vdhr_heat_change*86400.0
;;;;

;;;read Q1 monthly data from linshuheng
  fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/q1/q1"+"*") ;;;daily数据计算的q1

  f2 = addfiles(fils, "r")
  date = cd_calendar(f2[:]->time, 1)

  timeind = ind(date.le.201512.and.date.ge.197901)

  q1 = f2[:]->q1(timeind,:,:,:)

  q1 = (/q1*86400.0/)
 

;;;read Q2 monthly data from linshuheng
  fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/q2/q2"+"*") ;;;daily数据计算的q1

  f3 = addfiles(fils, "r")
  date = cd_calendar(f3[:]->time, 1)

  timeind = ind(date.le.201512.and.date.ge.197901)

  q2 = f3[:]->q2(timeind,:,:,:)

  q2 = (/q2*-86400.0/)
  cp = 1004
  q2 = q2/1004.0
;;;


;;;change for q1 and q2 

  q1_clm_change = get_data(q1,1)
  q2_clm_change = get_data(q2,1)

  q1_clm = q1_clm_change(0,:)
  q1_change = q1_clm_change(1,:)

  q2_clm = q2_clm_change(0,:)
  q2_change = q2_clm_change(1,:)
;;;;


;;; clm for latent heating 
    lat_heating_clm = dchr_heat_clm + schr_heat_clm + lschr_heat_clm
;;;


;;; vertical divergence of the vertical eddy transport of moisture

   ver_clm = q2_clm - lat_heating_clm
;;;





;;;;

    wks = gsn_open_wks("pdf","Q2以及潜热释放以及eddy输出的潜热的气候态的profile-ANN-WNP")                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@trXMaxF = 5
    resxy@trXMinF = -5
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =(/"black","red","green","blue"/)           ; change line color
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode    =  "MarkLines"

    resxy@xyDashPatterns = (/0,0,0,0/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           = (/1,1,1,1/)
    resxy@xyMarkerSizeF = 0.03
    resxy@xyMonoMarkerSize = True
    resxy@xyMarkerColors = resxy@xyLineColors   
    resxy@trYReverse        = True                       ; reverse Y-axis
    resxy@gsnXRefLine = 0

    plot1 = gsn_csm_xy(wks, (/ver_clm,q2_clm,lat_heating_clm, dchr_heat_clm/),level ,resxy)
    draw(plot1)
    frame(wks)
    ;;;
    end



  