begin 

    f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/review/new/Q-term-8.nc", "r")
    
    series_era_divQ_delta = f1->series_Q   
    series_dyna_adv_q = f1->series_dyna_adv_q
    series_dyna_con_q = f1->series_dyna_con_q
    series_thermaldyna_adv_q = f1->series_thermaldyna_adv_q
    series_thermaldyna_con_q =  f1->series_thermaldyna_con_q

    series_quad_adv_q = f1->series_quad_adv_q
    series_quad_con_q = f1->series_quad_con_q 
    series_eddy_q = f1->series_eddy_q 

    yearnum  = dimsizes(series_era_divQ_delta)

    q_delta_rc = regline(ispan(1,yearnum,1), series_era_divQ_delta)
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

    print("q_thermaldyna_adv_rc="+q_thermaldyna_adv_rc)
    
    q_thermaldyna_con_rc = regline(ispan(1,yearnum,1), series_thermaldyna_con_q)
    q_thermaldyna_con_rc = (/q_thermaldyna_con_rc*yearnum/)
    
    print("q_thermaldyna_con_rc= "+ q_thermaldyna_con_rc)
    
    q_quad_adv_rc = regline(ispan(1,yearnum,1), series_quad_adv_q)
    q_quad_adv_rc = (/q_quad_adv_rc*yearnum/)
    
    q_quad_con_rc = regline(ispan(1,yearnum,1), series_quad_con_q)
    q_quad_con_rc = (/q_quad_con_rc*yearnum/)

    
    q_eddy_rc = regline(ispan(1,yearnum,1), series_eddy_q)
    q_eddy_rc = (/q_eddy_rc*yearnum/)





;==============================================================
; Mission 5. bar plot
; ==============================================================
  
  ; x = (/1,2,3,4,5,6,7,8/)

  ; y = (/q_delta_rc,q_dyna_adv_rc,q_dyna_con_rc,q_thermaldyna_adv_rc,\
  ;      q_thermaldyna_con_rc,q_quad_adv_rc,q_quad_con_rc,q_eddy_rc/)
  ; labels = (/"TOTAL","ADV(DY)","CON(DY)","ADV(TH)","CON(TH)",\
  ; "ADV(QD)","CON(QD)","EDDY"/)

  x = (/1,2,3,4,5,6/)

  y = (/q_delta_rc,q_dyna_con_rc,q_dyna_adv_rc,q_thermaldyna_con_rc,q_thermaldyna_adv_rc,\
       q_eddy_rc/)

   labels = (/"DIVQ","DY(CON)","DY(ADV)","TH(CON)","TH(ADV)",\
    "EDDY"/)
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
    
  


