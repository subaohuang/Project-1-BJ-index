;;********************************************************
;;draw diabatic heat and Q1 HEAT rate profile averaged over SCS and 
;;philippines sea
;;diabatic and convection data from CFSR
;;q1 derived from NCEP_DOE calculate by eqution
;;                               lnshuheng 2020.1.1 
;**********************************************************
begin
; ==============================================================
; Mission 1. read data
; ==============================================================
;;read diabatic data 
  f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/diabatic_heating_CFSR_ANN_clm_change_profile.nc", "r")

  diabatic_heat_change = f1->diab_q(1,:)
  diabatic_heat_clm = f1->diab_q(0,:) 
  level = diabatic_heat_change&lev
  level = level/100.0
  
  diabatic_heat_change = (/diabatic_heat_change*86400/)
  diabatic_heat_clm = (/diabatic_heat_clm*86400/)
;;;

;;read deep convection data
  f2 = addfile("/home/ys17-19/lsh/Project/SCS-rain/deepcon_heating_CFSR_ANN_clm_change_profile.nc", "r")

  deepcon_heat_change = f2->deepcon_q(1,:)
  deepcon_heat_clm = f2->deepcon_q(0,:) 
  
  
  deepcon_heat_change = (/deepcon_heat_change*86400/)
  deepcon_heat_clm = (/deepcon_heat_clm*86400/)
;;;

;;;read Q1 monthly data from linshuheng
  fils = systemfunc("ls /home/yangsong3/data-observation/NCEP_NCAR/monthly/q1.1948-2014.nc")
  f3 = addfile(fils, "r")
  date = cd_calendar(f3->time, 1)

  timeind = ind(date.le.201312.and.date.ge.197901)

  q1 = f3->q1(timeind,:,:,:)
  q1 = (/q1*86400/)

  printVarSummary(q1)

  q1_annual = month_to_annual(q1, 1)
  ; q1_annual = month_to_season(q1, "MAM")

  dim = dimsizes(q1_annual)

  q1_annual_rc = regCoef_n(ispan(1,dim(0),1), q1_annual,0,0)
  copy_VarCoords(q1_annual(1,:,:,:), q1_annual_rc)  
  q1_annual_rc = q1_annual_rc*dim(0)

  q1_clm = dim_avg_n_Wrap(q1_annual, 0)


  ;;;区域平均
      g_latS=0.
      g_latN=15.
      g_lonW=110.
      g_lonE=160.

      q1_clim_area = wgt_areaave_Wrap(q1_clm(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
      q1_annual_rc_area = wgt_areaave_Wrap(q1_annual_rc(:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
      q1_change = q1_clim_area
      q1_change = (/q1_clim_area+q1_annual_rc_area/)
  ;;;
  

  q1_clim_area = q1_clim_area
  q1_change = q1_change
  ;;;
  

;;;;

  
    wks = gsn_open_wks("pdf","NCEP-NCAR_Q1,非绝热加热,深对流加热的profile-lindata-era")                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =  (/"blue","blue","red","red"/)          ; change line color
    resxy@xyMarkLineModes    =  (/"MarkLines","MarkLines","MarkLines","MarkLines"/)
    resxy@xyDashPatterns = (/0,1,0,1/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           = (/16,16,9,9/)
    resxy@xyMarkerColors =  (/"blue","blue","red","red"/)

   
    ; resxy@pmLegendDisplayMode = "Always"
    ; resxy@lgBoxLinesOn = False
    ; resxy@xyExplicitLegendLabels  =  (/"DPLH Clm","DPLH Clm+Change","DIABH Clm","DIABH Clm+Change"/)     ;加线的标注
    resxy@trYReverse        = True                       ; reverse Y-axis

    plot1 = gsn_csm_xy(wks, (/deepcon_heat_clm,deepcon_heat_change, \
                             diabatic_heat_clm,diabatic_heat_change/),level ,resxy)
    
    ;;;
    resxy@xyLineColors       :=  (/"black","black"/)          ; change line color
    resxy@xyMarkLineModes := (/"lines","lines"/)    
    ; resxy@xyExplicitLegendLabels := (/"Q1 Clim","Q1 Clim+Change"/)
    resxy@xyDashPatterns := (/0,1/)
    plot2 = gsn_csm_xy(wks, (/q1_clim_area,q1_change/),q1&lev, resxy)

    overlay(plot1, plot2)



   ;;添加legend 
      ;--------------------------------------------------------------------------
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF  = 2.                    ; line thicker

    
    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.01                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    res_lines@gsLineColor       = "black"

    yy = (/25.,25./)
    xx = (/2.5,3.0/)
    text_pos = 3.1
   
    dum3 = gsn_add_polyline(wks,plot1,xx,yy,res_lines)              ; add polyline
    dum4 = gsn_add_text(wks,plot1,"Q1 Clim",text_pos,yy(0),res_text); add text
 
    res_lines@gsLineDashPattern = 1.                    ; solid line
    yy = yy+50
    dum31 = gsn_add_polyline(wks,plot1,xx,yy,res_lines)              ; add polyline
    dum41 = gsn_add_text(wks,plot1,"Q1 Clim+Change",text_pos,yy(0),res_text); add text
    
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsMarkerIndex     = 16
    res_lines@gsMarkerSizeF     = 0.01
    res_lines@gsLineColor       = "blue"
    res_lines@gsMarkerColor = "blue"

    yy = yy+50
    dum5 = gsn_add_polyline(wks,plot1,xx,yy,res_lines)              ; add polyline
    dum6 = gsn_add_text(wks,plot1,"DPLH Clim",text_pos,yy(0),res_text); add text
    dum55 = gsn_add_polymarker(wks,plot1,xx,yy,res_lines)   

    res_lines@gsLineDashPattern = 1.
   
    yy = yy+50
    dum51 = gsn_add_polyline(wks,plot1,xx,yy,res_lines)              ; add polyline
    dum61 = gsn_add_text(wks,plot1,"DPLH Clim+Change",text_pos,yy(0),res_text); add text
    dum551 = gsn_add_polymarker(wks,plot1,xx,yy,res_lines)   

   
    res_lines@gsMarkerIndex     = 9
    res_lines@gsLineColor       = "red"
    res_lines@gsLineDashPattern = 0.
    res_lines@gsMarkerColor     = "red"
    res_lines@gsMarkerThicknessF   = 3.0
    yy = yy+50
    dum7 = gsn_add_polyline(wks,plot1,xx,yy,res_lines)              ; add polyline
    dum8 = gsn_add_text(wks,plot1,"DIABH Clim",text_pos,yy(0),res_text); add text
    dum77 = gsn_add_polymarker(wks,plot1,xx,yy,res_lines)   
    
    res_lines@gsLineDashPattern = 1.
    yy = yy+50
    dum71 = gsn_add_polyline(wks,plot1,xx,yy,res_lines)              ; add polyline
    dum81 = gsn_add_text(wks,plot1,"DIABH Clim+Change",text_pos,yy(0),res_text); add text
    dum771 = gsn_add_polymarker(wks,plot1,xx,yy,res_lines)   
    
    draw(plot1)
    frame(wks)
end 
