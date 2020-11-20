begin 

read_name = asciiread("/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/info.txt", -1,"string")
print(read_name)
model_name = read_name 

enino_intensity = asciiread("enino_intensity_diff_30model", -1,"float")

; dbitcz = asciiread("dbitcz_30model", -1,"float")

dbitcz = asciiread("ninoyear_spring_north_eastern_pacific_bais_30model", -1,"float")
; dbitcz = asciiread("north_eastern_pacific_baise_gradient_30model", -1,"float")

dim_model = dimsizes(enino_intensity)
;;;polyline 

r = escorc(dbitcz, enino_intensity)


rc = regline(dbitcz,enino_intensity)
line = rc*dbitcz + rc@yintercept
  
  
;;;draw plot 
   wks = gsn_open_wks("pdf", "东太平洋厄尔尼诺发展年春季ITCZbias-和冬季Enino强度bias的关系")
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
    res@trYMaxF=1.5
    res@trYMinF=-1.5
    res@tmYLValues    = ispan(-15,15,5)*0.1
    res@tmYLLabels    = ispan(-15,15,5)*0.1
    
    res@tmXTOn=False
    res@tmYROn=False
 
    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@trXMaxF=4.0
    res@trXMinF=-4.0
    res@tmXBValues    = fspan(-0.4,0.4,9) *10.0
    res@tmXBLabels    = fspan(-0.4,0.4,9) *10.0

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
    res@tiMainString=" diff in Spring Eastern Pacific ITCZ  and DJF Enino Intensity "
    res@tiMainFontHeightF = 0.015
    res@tiYAxisString="~F22~ Enino Intensity Bias [k]"
    res@tiYAxisFontHeightF=0.012
    
    res@tiXAxisString="~F22~ Enino develop spring ITCZ Bias [mm/day]  "

    res@tiXAxisFontHeightF=0.012
    
    plot=gsn_csm_xy(wks, dbitcz, enino_intensity, res)

    res_poly=True
    res_poly@gsLineThicknessF=3.0
    res_poly@gsLineColor="black"
    plot1=gsn_add_polyline(wks, plot, dbitcz, line, res_poly)

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

    txid = gsn_add_text(wks,plot,labels,dbitcz,enino_intensity+0.08,txres)

    polyline=True
    polyline@gsLineColor="black"
    polyline@gsLineThicknessF=1.0
    polyline@gsLineDashPattern = 0.0

    plotline1=gsn_add_polyline(wks, plot, x1, res@tmYLValues, polyline)

    plotline2=gsn_add_polyline(wks, plot, res@tmXBValues,y1, polyline)

    Polymarker=True
    Polymarker@gsMarkerColor="black"
    Polymarker@gsMarkerSizeF=10.0
    Polymarker@gsMarkerIndex=16

    plotmark=gsn_add_polymarker(wks, plot, dbitcz, enino_intensity, Polymarker)
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
        gsn_text_ndc (wks,"dryer",0.2,0.26,txres2)
        gsn_text_ndc (wks,"wetter",0.7,0.26,txres2)

        txres2@txFontHeightF = 0.015
        txres2@txFontColor="black"
        txres2@txJust="CenterCenter"
        txres2@txAngleF="90"
        gsn_text_ndc (wks,"weakening",0.08,0.40,txres2)
        gsn_text_ndc (wks,"strengthening",0.08,0.80,txres2)


        ;;;
                                        
        txres2@txFontColor="black"
        txres2@txAngleF="0"
        plot_text3=gsn_add_text(wks, plot, "Correlation="+sprintf ("%4.2f",r),2.5, 1.3, txres2)

        
        R2  = r^2
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
        textres@lgLabels = "regression"+"("+"R~S~2~N~="+sprintf("%4.2f",R2)+")"
        ; textres@lgLabels = "regression"+"("+"R~S~2~N~="+sprintf("%4.2f",R2)+")"

        ; plot = simple_legend(wks,plot,gres,lineres,textres)   

    draw(plot)
    frame(wks)


end 



