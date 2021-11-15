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
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 


begin
    ;;;read era data
    path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    var_name = (/"uwnd","vwnd","shum"/)
    filein = path_wind + var_name +".mon.mean.nc"

    ;;;region 
    lat1 = -30
    lat2 = 30
    lon1 = 200
    lon2 = 300
    ;;;

    ;;;time span
    startmon =198001
    endmon = 202012 

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))
    yearnum = endyear - startyear + 1
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    u = f1->uwnd(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    f2 = addfile(filein(1),"r")
    v = f2->vwnd(timeind,{1000:300},{lat1:lat2},{lon1:lon2})


    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    q = f3->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q = (/q/1000.0/) ;; convert to kg/kg

    level = q&level
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

    dp = dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm = era_qu
    dp_cfm = conform_dims(dimsizes(era_qu),dp,1)          ; confirm dp to caculate sum(var*dp)
   
    levelind = ind(level.le.1000.0.and.level.ge.300.0)
    
    era_wgt_qu  = (/dim_sum_n(dp_cfm(:,levelind,:,:)*era_qu(:,levelind,:,:),1)/)               ; sum(var*dp)
    era_wgt_qv = (/dim_sum_n(dp_cfm(:,levelind,:,:)*era_qv(:,levelind,:,:),1)/)               ; sum(var*dp)

    era_wgt_qu = (/era_wgt_qu/9.8/)
    era_wgt_qv = (/era_wgt_qv/9.8/)
    print(era_wgt_qv&lat)  

    ;;;;integration 
    q_col = q(:,0,:,:)
    q_col =  (/dim_sum_n(dp_cfm(:,levelind,:,:)*q(:,levelind,:,:),1)/)   
    
    q_col_t = q_col
    q_col_t = center_finite_diff_n(q_col, 1, False, 0, 0) 

;;;;



;;;parameter 
     pi = get_pi("float")
     dx = u&lat
     dx = 6378388.*cos(0.0174533*u&lat)
     dx!0 = "lat"
     copy_VarCoords(u(1,1,:,1), dx)

     dx = (/dx*2.5/) ;;;longitude range 2.5



     dy = (pi/180)*1*6378388.

     dy = (/dy*2.5/) ;;;;;;latitude range 2.5
     
     ; dz = dpres_plevel(depth, 45.58, 0, 0)
     ; print(dz)
     ; print(depth)
     ; dz= dz(::-1)
;;;;;


;;;局地项积分

    q_col_sum = q_col_t(:,1,1) 

    dx_2d = conform_dims(dimsizes(q_col_t(:,{0:20},{250:280})), dx({0:20}), (/1/))

    q_col_sum = (/dim_sum_n_Wrap(q_col_t(:,{0:20},{250:280})*dx_2d*dy,(/1,2/))/)

    q_col_sum = (/q_col_sum/(9.8*86400*30)/)  ;;; /month to /s 



    
    q_col_12 = split_12(q_col_sum)

    q_t_annual = month_to_annual(q_col_sum, 1)
    q_t_annual = (/dim_avg_n_Wrap(q_col_12(0:4,:), 0)/)
;;;;


   
;;;;
    q_t_rc  = regCoef_n(ispan(1, yearnum, 1), q_t_annual, 0, 0)
    q_t_rc = (/q_t_rc*yearnum/)
    print(q_t_rc)
  
;;;

; ;;;;积分边界
   

    boundary_south = era_wgt_qv(:,0,0)
    boundary_north = boundary_south

    boundary_west = era_wgt_qv(:,0,0)
    boundary_east = boundary_west

    boundary_south = (/dim_sum_n_Wrap(era_wgt_qv(:,{0},{250:280})*dx({0}),1)/)
    boundary_north = (/dim_sum_n_Wrap(era_wgt_qv(:,{20},{250:280})*dx({20}),1)/)
    boundary_west = (/dim_sum_n_Wrap(era_wgt_qu(:,{0:20},{250})*dy,1)/)
    boundary_east = (/dim_sum_n_Wrap(era_wgt_qu(:,{0:20},{280})*dy,1)/)
;;;;


;;;JFMAM MEAN 
   boundary_south_12 = split_12(boundary_south)
   boundary_north_12 = split_12(boundary_north)
   boundary_east_12 = split_12(boundary_east)
   boundary_west_12 = split_12(boundary_west)
;;;;




;;;;
   boundary_south_jfmam = dim_avg_n_Wrap(boundary_south_12(0:4,:),0)
   boundary_north_jfmam = dim_avg_n_Wrap(boundary_north_12(0:4,:),0)
   boundary_west_jfmam = dim_avg_n_Wrap(boundary_west_12(0:4,:),0)
   boundary_east_jfmam = dim_avg_n_Wrap(boundary_east_12(0:4,:),0)

;;;;;



;;;trend 
    boundary_south_rc  = regCoef_n(ispan(1, yearnum, 1), boundary_south_jfmam, 0, 0)
    boundary_south_rc = (/boundary_south_rc*yearnum/)

    boundary_north_rc  = regCoef_n(ispan(1, yearnum, 1), boundary_north_jfmam, 0, 0)
    boundary_north_rc = (/boundary_north_rc*yearnum/)

    boundary_east_rc  = regCoef_n(ispan(1, yearnum, 1), boundary_east_jfmam, 0, 0)
    boundary_east_rc = (/boundary_east_rc*yearnum/)
    
    boundary_west_rc  = regCoef_n(ispan(1, yearnum, 1), boundary_west_jfmam, 0, 0)
    boundary_west_rc = (/boundary_west_rc*yearnum/)
;;;;;


print(boundary_east_rc)
print(boundary_west_rc)
print(boundary_south_rc)
print(boundary_north_rc)


;;;;
    f3=addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/evapr_oaflux_chazhi.nc", "r")
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    eva = f3->evapr(timeind,:,:)
    ; eva = (/eva*10/365/) ;;convert to mm/day ;;untis : m/month 
    eva = (/(eva*0.01)/(365*86400)/) ;;convert cm/yr to m/s 
    eva = (/eva*1000.0/) ;;; water density 1000 kg/m3 convert to kg/m^2s
    ; eva = (/rmMonAnnCycTLL(eva)/)
  
;;;;

   dx_eva = eva&lat
   dx_eva = 6378388.*cos(0.0174533*eva&lat)
   dx_eva!0 = "lat"
   copy_VarCoords(eva(1,:,1), dx_eva)

   dx_eva = (/dx_eva*2.5/) ;;;longitude range 2.5

;;;
    eva_sum = eva(:,1,1) 

    dx_eva_2d = conform_dims(dimsizes(eva(:,{0:20},{250:280})), dx_eva({0:20}), (/1/))

    eva_sum = (/dim_sum_n_Wrap(eva(:,{0:20},{250:280})*dx_eva_2d*dy,(/1,2/))/)


    evp_12 = split_12(eva_sum)

    evp_annual = month_to_annual(eva_sum, 1)
    evp_annual = (/dim_avg_n_Wrap(evp_12(0:4,:), 0)/)

   
;;;;

    evp_rc  = regCoef_n(ispan(1, yearnum, 1), evp_annual, 0, 0)
    evp_rc = (/evp_rc*yearnum/)
    print(evp_rc)
;;;


    ; f2=addfile("/home/yangsong3/data-observation/NCEP_NCAR/monthly/prate.sfc.mon.mean.nc", "r")
    ; date:=cd_calendar(f2->time, 1)
    ; timeind:=ind(date.ge.startmon.and.date.le.endmon)
    ; p_r=f2->prate(timeind,{-90:90},{0:360}) ;;;untis : kg/m^2S

    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    ; f2=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p_r=f2->precip(timeind,{-90:90},{0:360})
    p_r = (/p_r*10^3*0.001/86400/)  ;;;untis : kg/m^2S

    ;;;处理数据
    ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

    yo = fspan(-90.0, 90.0, 73)
    xo = fspan(0.0, 357.5, 144)
    p = linint2_Wrap(p_r&lon, p_r&lat, p_r, True, xo, yo, 0)

    delete(p&LAT)
    delete(p&LON)
    p!1 = "lat"
    p!2 = "lon"
    p&lat = yo
    p&lon = xo
;;;;

    p_sum = p(:,1,1) 

    dx_2d = conform_dims(dimsizes(p(:,{0:20},{250:280})), dx({0:20}), (/1/))

    p_sum = (/dim_sum_n_Wrap(p(:,{0:20},{250:280})*dx_2d*dy,(/1,2/))/)


    p_12 = split_12(p_sum)

    p_annual = month_to_annual(p_sum, 1)
    p_annual = (/dim_avg_n_Wrap(p_12(0:4,:), 0)/)
;;;;


   
;;;;
    p_rc  = regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
    p_rc = (/p_rc*yearnum/)
    print(p_rc)
;;;
;;;

exit()


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
  wks   = gsn_open_wks("pdf", "水汽通量输送的趋势图")
  plot=new(5,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  ; res@cnFillPalette="MPL_RdBu"
  res@cnFillPalette = "precip_diff_12lev"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -40.0
  res@cnMaxLevelValF       =   40.0
  res@cnLevelSpacingF      = 10.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Horizontal"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             = 90
  ; res@mpMaxLonF             = 160
  
  ;;;global region
   res@mpMinLatF             = -10
   res@mpMaxLatF             = 50
   res@mpMinLonF             = 60
   res@mpMaxLonF             = 160
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
 
  res@vcRefLengthF             = 0.025     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.0
  res@vcVectorDrawOrder="PostDraw"
  res@vcRefMagnitudeF          = 20.0
  res@gsnLeftString = "Moisture Transport"


  
  plot(0)=gsn_csm_vector_scalar_map(wks,qu_rc ,qv_rc, q_rc, res)
  res@gsnLeftString = "Dynamic component"
  res@vcRefMagnitudeF          = 20.0
  res@cnMinLevelValF       =   -20.0
  res@cnMaxLevelValF       =   20.0
  res@cnLevelSpacingF      = 5.0
  plot(1)=gsn_csm_vector_scalar_map(wks,qu_dyna_rc ,qv_dyna_rc, q_dyna_rc, res)
  res@gsnLeftString = "Thermaldynamic component"
  res@vcRefMagnitudeF          = 4.0
  res@cnMinLevelValF       =   -4.0
  res@cnMaxLevelValF       =   4.0
  res@cnLevelSpacingF      = 1.0

; res@lbAutoManage = False
; res@lbLabelAutoStride = False
; res@lbLabelStride = 4
; ; res@lbLabelStrings =  sprintf("%2.1e", fspan(res@cnMinLevelValF,res@cnMaxLevelValF,17))
; ; m = res@lbLabelStrings 
; ; m(8) = "0"
; ; res@lbLabelStrings = m


  plot(2)=gsn_csm_vector_scalar_map(wks,qu_thermaldyna_rc ,qv_thermaldyna_rc, q_thermaldyna_rc, res)
  res@gsnLeftString = "Quadratic component"
  res@vcRefMagnitudeF          = 1.0
    
  plot(3)=gsn_csm_vector_scalar_map(wks,qu_quad_rc ,qv_quad_rc, q_quad_rc, res)
  res@vcRefMagnitudeF          = 4.0
  res@gsnLeftString = "Transient component"
    
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
    res_pro@cnLevels = (/0,0.2/)
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
 pres@gsnPanelLabelBar = False
 pres@lbOrientation = "vertical"
 pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)

 gsn_panel(wks, plot, (/3,2/), pres)
 

end 
    
   
   

  





    
   ;;;

   






