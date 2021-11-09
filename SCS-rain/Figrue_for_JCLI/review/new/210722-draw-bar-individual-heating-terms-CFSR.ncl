;;********************************************************
;;draw individual diabatic heat term profile averaged over SCS and 
;;philippines sea
;;diabatic and convection data from CFSR
;;                               lnshuheng 2020.1.2 
;**********************************************************
function vertical_avg(data)
begin
    psfc =1013.0
    plev = data&lev/100.0
    
    ptop = min(plev)
    dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
    q1dp = sum(data*dp)/sum(dp)  ;;;积分
return q1dp
end 

begin
; ==============================================================
; Mission 1. read data
; ==============================================================
;;read diabatic data
  f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/review/new/individual_heating_CFSR_MAM_clm_change_profile-WNP.nc", "r")

  lschr_heat_change = f1->lschr_q(1,{100000:10000})   ;;;large scale condensate heating
  lschr_heat_clm = f1->lschr_q(0,{100000:10000})   ;;;

  level = lschr_heat_clm&lev
  level = level/100.0

  dchr_heat_change = f1->dchr_q(1,{100000:10000})   ;;;Deep convective heating rate 
  dchr_heat_clm = f1->dchr_q(0,{100000:10000})   ;;;

  
  srhr_heat_change = f1->srhr_q(1,{100000:10000})   ;;;solar radiation heating rate 
  srhr_heat_clm = f1->srhr_q(0,{100000:10000})   ;;;
  
  
  lwrhr_heat_change = f1->lwrhr_q(1,{100000:10000})   ;;;long-wave radiative heating rate 
  lwrhr_heat_clm = f1->lwrhr_q(0,{100000:10000})   ;;;


  lschr_heat_clm = lschr_heat_clm*86400.0
  lschr_heat_change = lschr_heat_change*86400.0 
  lschr_heat_change =  lschr_heat_change  - lschr_heat_clm

  dchr_heat_clm = dchr_heat_clm*86400.0
  dchr_heat_change = dchr_heat_change*86400.0 
  dchr_heat_change = dchr_heat_change - dchr_heat_clm

  srhr_heat_clm = srhr_heat_clm*86400.0
  srhr_heat_change = srhr_heat_change*86400.0 
  srhr_heat_change = srhr_heat_change -srhr_heat_clm

  lwrhr_heat_clm = lwrhr_heat_clm*86400.0
  lwrhr_heat_change = lwrhr_heat_change*86400.0 
  lwrhr_heat_change = lwrhr_heat_change - lwrhr_heat_clm
;;;

;;;读入浅对流和扩散
  schr_heat_change = f1->schr_q(1,{100000:10000})   ;;;shallow convection rate
  schr_heat_clm = f1->schr_q(0,{100000:10000})   ;;;

  vdhr_heat_change = f1->vdhr_q(1,{100000:10000})   ;;;l
  vdhr_heat_clm = f1->vdhr_q(0,{100000:10000})   ;;;

  schr_heat_clm = schr_heat_clm*86400.0
  schr_heat_change = schr_heat_change*86400.0 
  schr_heat_change = schr_heat_change - schr_heat_clm

  vdhr_heat_clm = vdhr_heat_clm*86400.0
  vdhr_heat_change = vdhr_heat_change*86400.0 
  vdhr_heat_change = vdhr_heat_change - vdhr_heat_clm
;;;;

;;sum 
  sum_change = dchr_heat_change
  sum_change = dchr_heat_change +  schr_heat_change + lschr_heat_change + lwrhr_heat_change + srhr_heat_change + vdhr_heat_change
;;;


;;;Vertical average
   dchr_heat_change_ver = vertical_avg(dchr_heat_change)
   schr_heat_change_ver = vertical_avg(schr_heat_change)
   lschr_heat_change_ver= vertical_avg(lschr_heat_change)
   lwrhr_heat_change_ver = vertical_avg(lwrhr_heat_change)
   srhr_heat_change_ver = vertical_avg(srhr_heat_change)
   vdhr_heat_change_ver = vertical_avg(vdhr_heat_change)
   sum_change_ver = vertical_avg(sum_change)
;;;

 
number = (/sum_change_ver,dchr_heat_change_ver,schr_heat_change_ver,lschr_heat_change_ver,lwrhr_heat_change_ver,srhr_heat_change_ver,vdhr_heat_change_ver/)

; print(number)
percent = sprintf("%4.2f",100*number/sum_change_ver) + "%"
;;;;


;;;read Q1 monthly data from linshuheng
  fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/q1/q1"+"*") ;;;daily数据计算的q1

  f3 = addfiles(fils, "r")
  date = cd_calendar(f3[:]->time, 1)

  timeind = ind(date.le.201012.and.date.ge.197901)

  q1 = f3[:]->q1(timeind,{1000:100},:,:)

  q1 = (/q1*86400.0/)
  lev = tofloat(q1&lev)
  copy_VarAtts(q1&lev, lev)
  delete(q1&lev)
  lev = (/lev*100.0/)
  q1&lev = lev
  
  
  ; printVarSummary(q1)

  ; q1_annual = month_to_annual(q1, 1)
  q1_annual = month_to_season(q1,"MAM")

  dim = dimsizes(q1_annual)

  q1_annual_rc = (/regCoef_n(ispan(1,dim(0),1), q1_annual,0,0)/)
  copy_VarCoords(q1_annual(1,:,:,:), q1_annual_rc)  

  q1_annual_rc = (/q1_annual_rc*dim(0)/)
  q1_clm = dim_avg_n_Wrap(q1_annual, 0)

  ;;;区域平均
      g_latS=0.
      g_latN=20.
      g_lonW=125.
      g_lonE=160.

      q1_clim_area = wgt_areaave_Wrap(q1_clm(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
      q1_annual_rc_area = wgt_areaave_Wrap(q1_annual_rc(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
      q1_change = q1_clim_area
      q1_change = (/q1_annual_rc_area/)
  ;;;;
  
      q1_change_ver = vertical_avg(q1_change)

 percent_q1 = sprintf("%4.2f",100*q1_change_ver/sum_change_ver) + "%"
print(q1_change_ver)

print(sum_change_ver)


wks = gsn_open_wks("pdf","Barplot-individual-heating-trend-MAM-WNP")                  ; send graphics to PNG file
  
;==============================================================
; Mission 5. bar plot
; ==============================================================
     
     x = ispan(1,7,1)

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15


     sres@trXMinF = x(0)-0.5
     sres@trXMaxF = x(dimsizes(x)-1)+0.5
     sres@trYMinF = -0.2
     sres@trYMaxF = 0.6

     sres@gsnDraw = True
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.3           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = (/1,2,3,4,5,6,7/)
     sres@tmXBLabels = (/"DIABH","DPLH","SHLH","LSCH","LWRH","SWRH","VDH"/)
     ; sres@tmYLPrecision = 1
      ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
      ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01
     
     sres@tmXBLabelFontHeightF = 0.02
     sres@tmXTLabelFontHeightF = 0.02
     sres@tmYLLabelFontHeightF = 0.02
     
     sres@gsnRightString = ""

     sres@tiYAxisString = "K/day"
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnLeftString = "(b) Changes in vertically averaged heating terms "
     sres@gsnLeftStringFontHeightF = 0.02


     sres@gsnYRefLine = 0.      ; set the reference line equal to 0
    

     str = ispan(1,7,1)
     
     color = (/"black","blue","green","orange","pink","red","gray"/)


     sres@gsnXYBarChartColors2 = color


     x_plot = x - 0.25+0.15
     
     plot = gsn_csm_xy(wks,x_plot,number,sres)   ; draw each time series
     
     text = percent
     

     txres = True
     txres@txFontColor = "black"
     txres@txFontHeightF = 0.015

     add= gsn_add_text(wks, plot, text, x_plot+0.015, number+0.02, txres)

     draw(plot)
     frame(wks)      

end

  