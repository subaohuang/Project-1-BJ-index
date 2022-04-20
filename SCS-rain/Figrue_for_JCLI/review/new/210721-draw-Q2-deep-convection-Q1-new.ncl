;;********************************************************
;;draw different season long-term change diabatic heat and Q1 HEAT rate profile averaged over SCS and 
;;philippines sea
;;diabatic and convection data from CFSR
;;q1 derived from NCEP_DOE calculate by eqution
;;                               lnshuheng 2020.1.1 
;**********************************************************
begin

  season = (/"MAM","JJA","SON","DJF"/)
  ; file1 = "/home/ys17-19/lsh/Project/SCS-rain/season/philipian/diabatic_heating_CFSR_" +season+"_clm_change_profile_2010-philipian.nc"
  ; file2 = "/home/ys17-19/lsh/Project/SCS-rain/season/philipian/deepcon_heating_CFSR_" +season+ "_clm_change_profile_2010-philipian.nc"
  file1 = "/home/ys17-19/lsh/Project/SCS-rain/season/philipian/diabatic_heating_CFSR_" +season+"_clm_change_profile_2010-TWNP.nc"
  file2 = "/home/ys17-19/lsh/Project/SCS-rain/season/philipian/deepcon_heating_CFSR_" +season+ "_clm_change_profile_2010-TWNP.nc"
 
  plot1 = new(4,graphic)
  plot2 = new(4,graphic)
  ; wks = gsn_open_wks("pdf","Q1,非绝热加热,深对流加热的profile-lindata-ERA-2010-philipian-season")                  ; send graphics to PNG file
  wks = gsn_open_wks("pdf","Q1,非绝热加热,深对流加热的profile-lindata-ERA-2015-TWNP-season")                  ; send graphics to PNG file

  order  = (/"(a)","(b)","(c)","(d)"/)


do i=0,3
; ==============================================================
; Mission 1. read data
; ==============================================================
;;read diabatic data 
  f1 = addfile(file1(i), "r")

  diabatic_heat_change = f1->diab_q(1,:)
  diabatic_heat_clm = f1->diab_q(0,:) 
  level = diabatic_heat_change&lev
  level = level/100.0
  
  diabatic_heat_change = (/diabatic_heat_change*86400/)
  diabatic_heat_clm = (/diabatic_heat_clm*86400/)
  if i.eq.0 then
  print("*****diabatic-heating=****"+dim_avg(diabatic_heat_change({20000:85000}) - diabatic_heat_clm({20000:85000})))
  end if
;;;

;;read deep convection data
  f2 = addfile(file2(i), "r")

  deepcon_heat_change = f2->deepcon_q(1,:)
  deepcon_heat_clm = f2->deepcon_q(0,:) 
  
  
  deepcon_heat_change = (/deepcon_heat_change*86400/)
  deepcon_heat_clm = (/deepcon_heat_clm*86400/)
 
 if i.eq.0 then
  print("*****deepcon-heating=****"+dim_avg(deepcon_heat_change({20000:85000}) -deepcon_heat_clm({20000:85000})))
end if
;;;

;;;read Q1 monthly data from linshuheng
  fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/q2/q2"+"*") ;;;daily数据计算的q1
  print(fils)
  f3 = addfiles(fils, "r")
  date = cd_calendar(f3[:]->time, 1)

  timeind = ind(date.le.201512.and.date.ge.197901)

  q1 = f3[:]->q2(timeind,:,:,:)

  q1 = (/q1*-86400.0/)
  cp = 1004
  q1 = q1/1004.0

  ; printVarSummary(q1)

  ; q1_annual = month_to_annual(q1, 1)
  q1_annual = month_to_season(q1, season(i))

  if season(i).eq."DJF"
     q1_annual(0:dimsizes(q1_annual&time)-2,:,:,:) = q1_annual(1:,:,:,:) 
     q1_annual(dimsizes(q1_annual&time)-1,:,:,:) = q1(dimsizes(q1&time)-1,:,:,:)
  end if 


  dim = dimsizes(q1_annual)

  q1_annual_rc = (/regCoef_n(ispan(1,dim(0),1), q1_annual,0,0)/)

  ; q1_annual_mk=trend_manken(q1_annual, False, 0)
  ; q1_annual_rc = q1_annual_mk(1,:,:,:)
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
      q1_change = (/q1_clim_area+q1_annual_rc_area/)
if i.eq.0 then
 print("*****Q1-heating=****"+dim_avg(q1_annual_rc_area({200:850})))
end if


  ;;;
  
  q1_clim_area = q1_clim_area
  q1_change = q1_change
  ;;;

;;;;
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/2.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       :=  (/"blue","blue","red","red"/)          ; change line color
    resxy@xyMarkLineModes    :=  (/"MarkLines","MarkLines","MarkLines","MarkLines"/)
    resxy@xyDashPatterns := (/0,1,0,1/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           := (/16,16,7,7/)
    resxy@xyMarkerColors :=  (/"blue","blue","red","red"/)

    resxy@gsnLeftStringFontHeightF = 0.02
    resxy@gsnLeftString = order(i)+" " +season(i) + " DPLH DIABH Q1 Profile"

  

    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
    resxy@tmXBMajorThicknessF=1.0
    resxy@tmYLMajorThicknessF=1.0

    
    resxy@tmXBMode  =  "Explicit"
    resxy@trXMaxF =  6.
    resxy@trXMinF = -2.
    resxy@tmXBValues = fspan(-2,6,5)
    resxy@tmXBLabels = fspan(-2,6,5)
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = fspan(-2,6,17)
    resxy@tiYAxisString = ""
   
    ; resxy@pmLegendDisplayMode = "Always"
    ; resxy@lgBoxLinesOn = False
    ; resxy@xyExplicitLegendLabels  =  (/"DPLH Clm","DPLH Clm+Change","DIABH Clm","DIABH Clm+Change"/)     ;加线的标注
    resxy@trYReverse        = True                       ; reverse Y-axis

    plot1(i) = gsn_csm_xy(wks, (/deepcon_heat_clm,deepcon_heat_change, \
                             diabatic_heat_clm,diabatic_heat_change/),level ,resxy)
    
    ;;;
    resxy@xyLineColors       :=  (/"black","black"/)          ; change line color
    resxy@xyMarkLineModes := (/"lines","lines"/)    
    ; resxy@xyExplicitLegendLabels := (/"Q1 Clim","Q1 Clim+Change"/)
    resxy@xyDashPatterns := (/0,1/)
    plot2(i) = gsn_csm_xy(wks, (/q1_clim_area,q1_change/),q1&lev, resxy)

    overlay(plot1(i), plot2(i))

end do 


pres =True
; drawNDCGrid(wks)
gsn_panel(wks, plot1, (/2,2/), pres)

print("done stript")

  ;  ;;添加legend 
  ;     ;--------------------------------------------------------------------------
    ; res_lines                   = True                  ; polyline mods desired
    ; res_lines@gsLineDashPattern = 0.                    ; solid line
    ; res_lines@gsLineThicknessF  = 2.                    ; line thicker

    
    ; res_text                    = True                  ; text mods desired
    ; res_text@txFontHeightF      = 0.01                 ; change text size
    ; res_text@txJust             = "CenterLeft"          ; text justification

    ; res_lines@gsLineColor       = "black"

    ; yy = (/25.,25./)
    ; xx = (/2.5,3.0/)
    ; text_pos = 3.1
   
    ; dum3 = gsn_polyline_ndc(wks,xx,yy,res_lines)              ; add polyline
    ; dum4 = gsn_text_ndc(wks,"Q1 Clim",text_pos,yy(0),res_text); add text
 
    ; res_lines@gsLineDashPattern = 1.                    ; solid line
    ; yy = yy+50
    ; dum31 = gsn_polyline_ndc(wks,xx,yy,res_lines)              ; add polyline
    ; dum41 = gsn_text_ndc(wks,"Q1 Clim+Change",text_pos,yy(0),res_text); add text
    
    ; res_lines@gsLineDashPattern = 0.                    ; solid line
    ; res_lines@gsMarkerIndex     = 16
    ; res_lines@gsMarkerSizeF     = 0.01
    ; res_lines@gsLineColor       = "blue"
    ; res_lines@gsMarkerColor = "blue"

    ; yy = yy+50
    ; dum5 = gsn_polyline_ndc(wks,xx,yy,res_lines)              ; add polyline
    ; dum6 = gsn_text_ndc(wks,"DPLH Clim",text_pos,yy(0),res_text); add text
    ; dum55 = gsn_polymarker_ndc(wks,xx,yy,res_lines)   

    ; res_lines@gsLineDashPattern = 1.
   
    ; yy = yy+50
    ; dum51 = gsn_polyline_ndc(wks,xx,yy,res_lines)              ; add polyline
    ; dum61 = gsn_text_ndc(wks,"DPLH Clim+Change",text_pos,yy(0),res_text); add text
    ; dum551 = gsn_polymarker_ndc(wks,xx,yy,res_lines)   

   
    ; res_lines@gsMarkerIndex     = 9
    ; res_lines@gsLineColor       = "red"
    ; res_lines@gsLineDashPattern = 0.
    ; res_lines@gsMarkerColor     = "red"
    ; res_lines@gsMarkerThicknessF   = 3.0
    ; yy = yy+50
    ; dum7 = gsn_polyline_ndc(wks,xx,yy,res_lines)              ; add polyline
    ; dum8 = gsn_text_ndc(wks,"DIABH Clim",text_pos,yy(0),res_text); add text
    ; dum77 = gsn_polymarker_ndc(wks,xx,yy,res_lines)   
    
    ; res_lines@gsLineDashPattern = 1.
    ; yy = yy+50
    ; dum71 = gsn_polyline_ndc(wks,xx,yy,res_lines)              ; add polyline
    ; dum81 = gsn_text_ndc(wks,"DIABH Clim+Change",text_pos,yy(0),res_text); add text
    ; dum771 = gsn_polymarker_ndc(wks,xx,yy,res_lines)       


end 