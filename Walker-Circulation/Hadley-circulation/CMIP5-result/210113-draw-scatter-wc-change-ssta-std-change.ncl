begin 

    read_name = asciiread("/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/info.txt", -1,"string")
    print(read_name)
    model_name = read_name 

    sta_std = asciiread("sta_std_change_30model", -1,"float")

    ; sta_std = where(sta_std.le.0, sta_std@_FillValue, sta_std)
    

    HC_change = asciiread("WC_change_30model_150-210", -1,"float")
    ; HC_change = where(sta_std.eq.sta_std@_FillValue, sta_std@_FillValue, HC_change)

    ; HC_change = (/-1*HC_change/)


    dim_model = dimsizes(HC_change)
    ;;;polyline 


    rc = regline(HC_change,sta_std)
    line = rc*HC_change + rc@yintercept

    r = escorc(sta_std, HC_change)


  
;;;draw plot 
   wks = gsn_open_wks("pdf", "WC_change&sta_std_change_30model")
    res                   = True                     ; plot mods desired
    res@gsnDraw=False
    res@gsnFrame=False
    res@xyLineColors     = (/"red"/)          ; change line color
    res@xyMarkerColor ="white"
    res@xyMarkLineModes   = (/"Markers","Lines"/)
    res@xyLineThicknessF=3.0
    res@xyDashPatterns=1
    ; res@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    res@xyMarkers         =  16                      ; choose type of marker  
    res@xyMarkerSizeF     = 0.0001      ; Marker size (default 0.01)

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@trYMaxF = 0.3
    res@trYMinF = -0.3
    res@tmYLValues    = ispan(-30,30,5)*0.01
    res@tmYLLabels    = ispan(-30,30,5)*0.01
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    ; res@trXMaxF =  4.0
    ; res@trXMinF = -2.0
    ; res@tmXBValues    = ispan(-2,4,1) 
    ; res@tmXBLabels    = ispan(-2,4,1) 
   
    res@trXMaxF =  0.3
    res@trXMinF = -0.3
    res@tmXBValues    =  ispan(-30,30,5)*0.01
    res@tmXBLabels    =ispan(-30,30,5)*0.01


    res@tmXMajorGrid                = True          ; implement x grid 
    res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
    res@tmYMajorGrid                = True          ; implement y grid
    res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
    res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
    res@tmXBLabelFontHeightF=0.015

    res@tmXBMajorLengthF=0.008        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0.008  ;调整刻度的长度，注意2个命令要同时使用
    res@tmXUseBottom=False
    res@tmXBLabelDeltaF=-0.8
    res@tmBorderThicknessF=2.0

    res@tmYLLabelFontHeightF=0.015
    res@vpXF=0.15
    res@vpYF=0.9
    res@vpHeightF=0.6
    res@vpWidthF=0.6
    
    ; res@tiMainString=" Spring Double ITCZ Bias and DJF Enino Intensity"
    res@tiMainString=" WC change and SSTA std change"
    res@tiMainFontHeightF = 0.015
    res@tiYAxisString="~F22~ SSTA std change[k]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~ WC change [kg^-3]  "

    res@tiXAxisFontHeightF=0.012
    
    plot=gsn_csm_xy(wks, HC_change, sta_std, res)

    res_poly=True
    res_poly@gsLineThicknessF=3.0
    res_poly@gsLineColor="black"
    plot1=gsn_add_polyline(wks, plot, HC_change, line, res_poly)

    ; ; 90~S~o~N~E

    ; Attach some text strings next to the markers.
    txres               = True
    txres@txFontHeightF = 0.01
    txres@txFont = "times-roman"
    txres@txFontColor="blue"
    txres@tfPolyDrawOrder="PostDraw"
    txres@txJust = "TopCenter"
    ; txres@txp
    labels            = "  " + ispan(1,dim_model,1)
    ; labels=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","M"/)
    ; labels=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/)
    x1=conform(res@tmYLValues,0,0)
    y1=conform(res@tmXBValues,0,0)

    txid = gsn_add_text(wks,plot,labels,HC_change,sta_std+0.015,txres)

    polyline=True
    polyline@gsLineColor="black"
    polyline@gsLineThicknessF=1.0

    plotline1=gsn_add_polyline(wks, plot, x1, res@tmYLValues, polyline)

    plotline2=gsn_add_polyline(wks, plot, res@tmXBValues,y1, polyline)

    Polymarker=True
    Polymarker@gsMarkerColor="black"
    Polymarker@gsMarkerSizeF=5.0
    Polymarker@gsMarkerIndex=16

    plotmark=gsn_add_polymarker(wks, plot, HC_change, sta_std, Polymarker)
    ; Polymarker@gsMarkerColor="red"
    ; plotmark2=gsn_add_polymarker(wks, plot, rc_width_trend_rea_mean, rc_sst_rea_mean, Polymarker)
     

     ;;;添加最右边模式名称
        xtxt=0.8
        ytxt = fspan(0.3,0.9,dim_model)


        labels = ""+ispan(1,dim_model,1) + ":" +model_name + ""  

        ; alfpe= (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","A","B","C","D"/)
        ; alfpe=(/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r"/)
        ; labels = alfpe+ ":" +models

        txres2               = True
        txres2@txFontHeightF = 0.015
        txres2@txFontColor="black"
        txres2@txJust="TopLeft"
        do i=0,dim_model-1
         gsn_text_ndc (wks,labels(i),xtxt,ytxt(dim_model-1-i),txres2)
        end do
        txres2@txFontThicknessF=3.0
        gsn_text_ndc (wks,"CMIP5",xtxt,0.95,txres2)



        ;;;添加x轴下标
        txres2@txFontColor="black"
        txres2@txJust="CenterCenter"
        gsn_text_ndc (wks,"weakening",0.2,0.26,txres2)
        gsn_text_ndc (wks,"strengthening",0.7,0.26,txres2)

        ; txres2@txFontHeightF = 0.015
        ; txres2@txFontColor="black"
        ; txres2@txJust="CenterCenter"
        ; txres2@txAngleF="90"
        ; gsn_text_ndc (wks,"weaker",0.08,0.40,txres2)
        ; gsn_text_ndc (wks,"stronger",0.08,0.80,txres2)


        ;;;
                                        
        txres2@txFontColor="black"
        txres2@txAngleF="0"
        plot_text3=gsn_add_text(wks, plot, "Correlation="+sprintf ("%4.2f",r),1, 3.35, txres2)



        ;;legend 
        gres = True
        gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
        gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

        lineres = True
        lineres@lgLineColors = (/"black"/) ; line colors
        lineres@lgLineThicknesses = 3.0                       ; line thicknesses
        lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

        textres = True
        textres@lgLabels = "Slope="+sprintf("%4.2f",rc)

        plot = simple_legend(wks,plot,gres,lineres,textres)   

    draw(plot)
    frame(wks)


end 



