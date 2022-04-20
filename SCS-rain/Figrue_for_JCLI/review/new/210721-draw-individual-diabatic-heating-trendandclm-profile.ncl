;;********************************************************
;;draw individual diabatic heat term profile averaged over SCS and 
;;philippines sea
;;diabatic and convection data from CFSR
;;                               lnshuheng 2020.1.2 
;**********************************************************
begin
; ==============================================================
; Mission 1. read data
; ==============================================================
;;read diabatic data
  f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/review/new/individual_heating_CFSR_MAM_clm_change_profile-WNP.nc", "r")

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

    wks = gsn_open_wks("pdf","非绝热加热的各项的trend以及气候态的profile-MAM-WNP")                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@trXMaxF = 5
    resxy@trXMinF = -5
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =(/"blue","blue","green","green","orange","orange","pink","pink","red","red","gray","gray"/)           ; change line color
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode    =  "Lines"

    resxy@xyDashPatterns = (/0,1,0,1,0,1,0,1,0,1,0,1/)
    ; resxy@xyMonoMarker        = False
    ; resxy@xyMarkers           = (/0,0,1,1,2,2,3,3,4,4,5,5/)
    ; resxy@xyMarkerColors =  (/"blue","blue","red","red","green","green","orange","orange","gray","gray","pink","pink"/)    
    resxy@trYReverse        = True                       ; reverse Y-axis
    resxy@gsnXRefLine = 0



    resxy@tmXBLabelFontHeightF=0.02
    resxy@tmYLLabelFontHeightF=0.02
    ; resxy@tmXBMajorThicknessF=1.0
    ; resxy@tmYLMajorThicknessF=1.0

    
    resxy@tmXBMode  =  "Explicit"
    
    resxy@tmXBValues = fspan(-5,5,11)
    resxy@tmXBLabels = fspan(-5,5,11)
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = fspan(-5,5,21)
    resxy@tiYAxisString = "Pressure level (hPa)"

    resxy@tiXAxisString = "heating rate (k/day)"
    resxy@tiYAxisFontHeightF = 0.02
    resxy@tiXAxisFontHeightF = 0.02

    resxy@tmYRMode  =  "Explicit"
    resxy@trYMaxF = 1000.0
    resxy@trYMinF = 0.0

    resxy@tmYLValues = ispan(0,1000,200)
    resxy@tmYLLabels = fspan(0,1000,200)

    resxy@tmYLMinorOn = True
    resxy@tmYLMinorValues = ispan(0,1000,50)

    resxy@vpYF = 0.95
    resxy@vpXF = 0.15

    resxy@tmYROn = False
    resxy@tmXTOn = False

    resxy@gsnLeftString = "(a) Changes in profile of heating terms "
    resxy@gsnLeftStringFontHeightF = 0.02


    plot1 = gsn_csm_xy(wks, (/ dchr_heat_clm,dchr_heat_change,\
                             schr_heat_clm,schr_heat_change,\
                             lschr_heat_clm,lschr_heat_change,\
                             lwrhr_heat_clm,lwrhr_heat_change,\
                             srhr_heat_clm,srhr_heat_change,\
                             vdhr_heat_clm,vdhr_heat_change/),level ,resxy)
     
;     ;;;draw legend

    ; ;***********************************************
    ; ; legend resources
    ; ;***********************************************  
    res_lines = True
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF = 3.0
    
    color = (/"blue","green","orange","pink","red","gray"/)    
    text = (/"Deep convection","Shallow convection","Large scale condensate","Longwave radiation","Shortwave radiation","Vertical diffusion"/)
    yloc =fspan(0.9, 0.7, 6)  
    x  = (/0.76,0.76+0.05/)
    


    do i=0,5
     y = x
     y = yloc(i)
     res_lines@gsLineColor       = color(i)
     gsn_polyline_ndc(wks, x, y, res_lines)
     delete(y)
    end do

    textres = True
    textres@txFontHeightF = 0.012
    textres@txJust = "CenterLeft"


   do i = 0,5
     y = x
     y = yloc(i)
    gsn_text_ndc(wks, text(i), x(dimsizes(x)-1) + 0.01, y , textres)
   end do
  
    


;;;;

    ; drawNDCGrid(wks)


    draw(plot1)
    frame(wks)


    ;;;
    end



  