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
  f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/individual_heating_CFSR_ANN_clm_change_profile.nc", "r")

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

    wks = gsn_open_wks("pdf","非绝热加热的各项的trend以及气候态的profile-ANN-WP")                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@trXMaxF = 5
    resxy@trXMinF = -5
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =(/"blue","blue","red","red","green","green","orange","orange","gray","gray","pink","pink"/)           ; change line color
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode    =  "MarkLines"

    resxy@xyDashPatterns = (/0,1,0,1,0,1,0,1,0,1,0,1/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           = (/0,0,1,1,2,2,3,3,4,4,5,5/)
    resxy@xyMarkerColors =  (/"blue","blue","red","red","green","green","orange","orange","gray","gray","pink","pink"/)    
    resxy@trYReverse        = True                       ; reverse Y-axis
    resxy@gsnXRefLine = 0

    plot1 = gsn_csm_xy(wks, (/lschr_heat_clm,lschr_heat_change,\
                             dchr_heat_clm,dchr_heat_change,\
                             srhr_heat_clm,srhr_heat_change,\
                             lwrhr_heat_clm,lwrhr_heat_change,\
                             schr_heat_clm,schr_heat_change,\
                             vdhr_heat_clm,vdhr_heat_change/),level ,resxy)
    draw(plot1)
    frame(wks)
    ;;;
    end



  