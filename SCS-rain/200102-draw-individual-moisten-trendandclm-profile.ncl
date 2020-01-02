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
  f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/individual_moistening_CFSR_ANN_clm_change_profile.nc", "r")

  lscmr_heat_change = f1->lscmr_q(1,:)   ;;;large scale condensate heating
  lscmr_heat_clm = f1->lscmr_q(0,:)   ;;;

  level = lscmr_heat_clm&lev
  level = level/100.0

  dcmr_heat_change = f1->dcmr_q(1,:)   ;;;Deep convective heating rate 
  dcmr_heat_clm = f1->dcmr_q(0,:)   ;;;

  
  

  lscmr_heat_clm = lscmr_heat_clm*86400.0
  print(lscmr_heat_clm)
  exit()
  lscmr_heat_change = lscmr_heat_change*86400.0

  dcmr_heat_clm = dcmr_heat_clm*86400.0
  dcmr_heat_change = dcmr_heat_change*86400.0

  ;;

;;;读入浅对流和扩散
  scmr_heat_change = f1->scmr_q(1,:)   ;;;shallow convection rate
  scmr_heat_clm = f1->scmr_q(0,:)   ;;;

  vdmr_heat_change = f1->vdmr_q(1,:)   ;;;l
  vdmr_heat_clm = f1->vdmr_q(0,:)   ;;;

  scmr_heat_clm = scmr_heat_clm*86400.0
  scmr_heat_change = scmr_heat_change*86400.0

  vdmr_heat_clm = vdmr_heat_clm*86400.0
  vdmr_heat_change = vdmr_heat_change*86400.0

  
;;;;

    wks = gsn_open_wks("pdf","凝结的各项的trend以及气候态的profile-ANN")                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    ; resxy@trXMaxF = 5
    ; resxy@trXMinF = -5
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =(/"blue","blue","red","red","green","green","orange","orange"/)           ; change line color
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode    =  "MarkLines"

    resxy@xyDashPatterns = (/0,1,0,1,0,1,0,1,0,1,0,1/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           = (/0,0,1,1,2,2,3,3,4,4,5,5/)
    resxy@xyMarkerColors =  (/"blue","blue","red","red","green","green","orange","orange"/)    
    resxy@trYReverse        = True                       ; reverse Y-axis
    resxy@gsnXRefLine = 0

    plot1 = gsn_csm_xy(wks, (/lscmr_heat_clm,lscmr_heat_change,\
                             dcmr_heat_clm,dcmr_heat_change,\
                             scmr_heat_clm,scmr_heat_change,\
                             vdmr_heat_clm,vdmr_heat_change/),level ,resxy)
    draw(plot1)
    frame(wks)
    ;;;
    end



  