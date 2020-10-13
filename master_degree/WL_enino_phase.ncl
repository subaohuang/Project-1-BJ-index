
;;************************************
;;WL index 在不同Enino相位中的变化
;**************************************



begin 
;;---------------------------------------
;;read data
;----------------------------------------
locationxi = asciiread("/mnt/f/share/rea_msf/era_1.5wlocationxi.ascii", 456, "float")
locationxi = (/dtrend_msg_n(ispan(1,456,1),locationxi,False,True,0)/)
locationxi = rmAnnCycle1D(locationxi)
locationxi!0 = "time"
locationxi&time = ispan(1,456,1)
; locationxi_anu = month_to_annual(locationxi, 1)

Enino_devo_year = (/1982,1991,1997/) 
Enino_decay_year = (/1983,1992,1998/) 

locationxi_com = new((/dimsizes(Enino_devo_year),24/),typeof(locationxi))


do i = 0,2
  locationxi_com(i,:) = locationxi((Enino_devo_year(i)-1979)*12:(Enino_decay_year(i)-1979+1)*12-1)
end do 

x = array_append_record(ispan(1,12,1), ispan(1, 12, 1), 0)

;;;read sst and calculate correlation
    f2=addfile("/mnt/f/share/ersst/HadISST_sst.nc", "r")
    time = f2->time
    date = cd_calendar(time, 1)
    timeind=ind(date.ge.197901.and.date.le.201612)
    sst=lonFlip(f2->sst(timeind,{30:-30},:))

    sst = (/dtrend_msg_n(ispan(1,456,1),sst,False,True,0)/)
    sst_a = rmMonAnnCycTLL(sst)
    nino3 = dim_avg_n_Wrap( sst_a(:,{-5:5},{210:270}) , (/1,2/))
    
    nino3_com = new((/dimsizes(Enino_devo_year),24/),typeof(nino3))

    do i = 0,2
     nino3_com(i,:) = nino3((Enino_devo_year(i)-1979)*12:(Enino_decay_year(i)-1979+1)*12-1)
    end do 
;;;

;; corr 
   
;*************************************************** 超前滞后相关
    mxlag    = 12
    x_Lead_y = esccr(nino3_com(2,:),locationxi_com(2,:), mxlag)
    y_Lead_x = esccr(locationxi_com(2,:),nino3_com(2,:), mxlag)    ; switch the order of the series

    ccr = new ( 2*mxlag+1, double)    
    ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
    ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
    x1 = ispan(-12,12,1)     
    wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/不同ENSO相位下位置和强度相关性")
    res=True
    res@gsnXRefLine = 0
    plot1=gsn_csm_xy(wks1,x1,ccr,res)

; exit()



;;---------------------------------------
;;draw plot 
;----------------------------------------
     

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/不同Enino相位下变化的nio3指数")

    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@xyLineThicknesses = (/ 2.0,2.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black"  ,"black","black"/)          ; change line color
    resxy@xyDashPatterns = (/2,1,0/)
    resxy@xyMarkerColor = "Black"
    resxy@xyMonoMarkLineMode = True
    resxy@xyMarkLineMode   = (/"MarkLines"/) 
    

    resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    resxy@xyMarkers         =  (/4,5,16/)                    ; choose type of marker  
    resxy@xyMarkerSizeF     = 0.01        ; Marker size (default 0.01)
    
  

    
    ; resxy@gsnRightString=""+slp_data(3,1)+""
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01

    resxy@trXMaxF = 24
    resxy@trXMinF = 1
    
    resxy@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    resxy@tmXBValues=ispan(1,24,1)
    resxy@tmXBLabels = ""+x

    
    resxy@tmXBLabelFontAspectF=1.0
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  

    resxy@tmXBMinorLengthF=0.003
    resxy@tmXBMinorOutwardLengthF=0.003
    
    resxy@tmYLMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    resxy@tmYLLabelFont=30
    resxy@tmYLLabelFontHeightF=0.008
    resxy@tmYLMinorOn=False
    resxy@tmXTOn=False
    resxy@tmYROn=False
    

    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
    resxy@pmLegendSide           = "Top"               ; Change location of 
    resxy@pmLegendParallelPosF   = .2                  ; move units right
    resxy@pmLegendOrthogonalPosF = -0.3                ; move units down
    resxy@pmLegendWidthF         = 0.15                ; Change width and
    resxy@pmLegendHeightF        = 0.18                ; height of legend.
    resxy@lgPerimOn              = False               ; turn off box around
    resxy@lgLabelFontHeightF     = .025                ; label font height
    resxy@xyExplicitLegendLabels = (/"82-83","91-92","97-98"/)         ; create explicit labels


    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.6
    ; resxy@vpXF=0.15
    ; resxy@vpYF=0.95

    plot = gsn_csm_xy(wks,ispan(1,24,1),nino3_com,resxy)


    draw(plot)
    frame(wks)

    end 