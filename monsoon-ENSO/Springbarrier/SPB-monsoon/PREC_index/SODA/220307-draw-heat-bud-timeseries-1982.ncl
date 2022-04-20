begin 

    f4 = addfile("hb_soda.nc", "r")
;====================================================

    time := f4->time 
    date:=cd_calendar(time, 1)
    print(date(0))
    ; timeind:=ind(date.le.endmon.and.date.ge.startmon)
    timeind:=ind(date.le.198212.and.date.ge.198201)
    ; timeind:=ind(date.eq.1982)
    dTdt_ct = f4->dTdt(timeind,:,:)

    UATXB_ct = f4->UATXB(timeind,:,:)  ;;;; ZA 256
    VATYB_ct = f4->VATYB(timeind,:,:)

    WATZB_ct = f4->WATZB(timeind,:,:)  ;;;EK

    UBTXA_ct = f4->UBTXA(timeind,:,:)  ;;MEAN zonal advection
    VBTYA_ct = f4->VBTYA(timeind,:,:)
    WBTZA_ct = f4->WBTZA(timeind,:,:) ;;; TH
    Q_ct = f4->Q(timeind,:,:)

    UATXA_ct = f4->UATXA(timeind,:,:)
    VATYA_ct = f4->VATYA(timeind,:,:)
    WATZA_ct = f4->WATZA(timeind,:,:)

    NLT = UATXA_ct
    NLT = (/UATXA_ct+VATYA_ct+WATZA_ct/)
;======================================================
;;;nino34 
;=====================================================
  lat1 = -5
  lat2 = 5
  lon1 = 210
  lon2 = 270
  ;;;regional mean 
    ; dTdt_ct_nino3  =  dim_avg_n_Wrap(dTdt_ct(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
 
    dTdt_ct_nino3  :=  avg(dTdt_ct(5:8,{lat1:lat2},{lon1:lon2}))
    UATXB_ct_nino3  =  avg(UATXB_ct(5:8,{lat1:lat2},{lon1:lon2}))
    VATYB_ct_nino3  =  avg(VATYB_ct(5:8,{lat1:lat2},{lon1:lon2}))
    WATZB_ct_nino3 =  avg(VATYB_ct(5:8,{lat1:lat2},{lon1:lon2})) ;;;EK
    UBTXA_ct_nino3 = avg(UBTXA_ct(5:8,{lat1:lat2},{lon1:lon2})) ;;MEAN zonal advection
    VBTYA_ct_nino3 = avg(VBTYA_ct(5:8,{lat1:lat2},{lon1:lon2})) ;;MEAN zonal advection
    WBTZA_ct_nino3 = avg(WBTZA_ct(5:8,{lat1:lat2},{lon1:lon2})) ;;MEAN zonal advection
    Q_ct_nino3 = avg(Q_ct(5:8,{lat1:lat2},{lon1:lon2})) ;;MEAN zonal advection
    NLT_nino3 = avg(NLT(5:8,{lat1:lat2},{lon1:lon2})) ;;MEAN zonal advection
;;;
   hb_mean_mamjja = new(9,float)
   hb_mean_mamjja(0) = dTdt_ct_nino3
   hb_mean_mamjja(1) = UATXB_ct_nino3
   hb_mean_mamjja(2) = VATYB_ct_nino3 
   hb_mean_mamjja(3) =  WATZB_ct_nino3
   
   hb_mean_mamjja(4) =  UBTXA_ct_nino3
   hb_mean_mamjja(5) = VBTYA_ct_nino3
   hb_mean_mamjja(6) = WBTZA_ct_nino3
    hb_mean_mamjja(7) =  Q_ct_nino3
    hb_mean_mamjja(8) =  NLT_nino3 
;;;;


     wks = gsn_open_wks("pdf","El Nina-NOAA-YEAR-BAR-JJASON-MEAN-HEATBUDGET-TERMS")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15


     sres@trXMinF = 0.5
     sres@trXMaxF = 9.5
     sres@trYMinF = -0.6
     sres@trYMaxF = 0.6

     sres@gsnDraw = False
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.3           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = ispan(1, 9, 1)
     sres@tmXBLabels = (/"dTdt","ZA","VA","EK","MAU","MAV","TH","NHF","NDH"/)
     ; sres@tmYLPrecision = 1
      ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
      ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01
     
     sres@tmXBLabelFontHeightF = 0.0205
     sres@tmXTLabelFontHeightF = 0.0205
     sres@tmYLLabelFontHeightF = 0.0225
     sres@tiMainFontHeightF = 0.025
     sres@gsnRightString = ""
     sres@tiYAxisString = "[~S~o~N~C~N~/month]" 
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnYRefLine = 0.      ; set the reference line equal to 0

     str = ispan(1,2,1)
     ; color = (/"black","gray","yellow","green","pink","deeppink3","red","blue"/)
     color = (/"black"/)
    
     sres@gsnXYBarChartColors = color

     sres@gsnLeftStringFontHeightF = 0.02
     ; sres@gsnLeftString = "AMJJASOND mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino1+2 "
     ; sres@gsnLeftString = "JJASON mean Heat Budget Nino3 "
     sres@gsnLeftString = "GODAS JJASON mean Heat Budget"
     ; plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i/),number(:,i),sres)   ; draw each time series

     plot1 = gsn_csm_xy(wks,ispan(1,9,1)-0.15,hb_mean_mamjja(:),sres)   ; draw each time series     
    ;  color = (/"red"/)
    ;  sres@gsnXYBarChartColors = color
    ;  plot2 = gsn_csm_xy(wks,ispan(1,9,1)+0.15,hb_mean_mamjja(1,:),sres)   ; draw each time series     
    ; 
    draw(plot1)
    frame(wks)
end