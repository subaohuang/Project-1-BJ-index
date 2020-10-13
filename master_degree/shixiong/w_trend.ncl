begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_iod_pacific_w_u.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)
w_r=f1->w(tind,{500},{10:-10},{40:280})
w=(w_r@scale_factor*w_r+w_r@add_offset)*(-100.0)
copy_VarCoords(w_r, w)


w_mean=dim_avg_n_Wrap(w, (/0,1/))
w&time=ispan(1,432,1)
w_trend=regCoef_n(w&time, dim_avg_n_Wrap(w, 1), 0, 0)
printVarSummary(w_trend)
w_trend!0="longitude"
printVarSummary(w)

;;;********************************************************
;;;画图
;;；********************************************************
 resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnYRefLine=(/0/)
    resxy@gsnYRefLineThicknesses=1.5
    resxy@tmYLMode="Explicit"
    resxy@trYMaxF=0.01
    resxy@trYMinF=-0.01
    resxy@tmYLValues=fspan(-0.01,0.01,5)
    resxy@tmYLLabels=fspan(-0.01,0.01,5)
    resxy@tmYLMinorOn=False
    
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=280
    resxy@tmXBValues=ispan(40,280,60)
    resxy@tmXBLabels=(/"40~S~o~N~E","100~S~o~N~E","160~S~o~N~E","140~S~o~N~W","80~S~o~N~W"/)
    resxy@tmXBMinorValues=ispan(40,280,20)
    
    resxy@vpHeightF=0.3
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9
    ; resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColor      = (/"blue"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines"/)
    resxy@gsnLeftString="(a)"
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    ; resxy@xyDashPatterns    =(/0,0,2/)
    
    ; resxy@tiYAxisString ="nino3.4" 
    resxy@tmXBLabelFontAspectF=2.5
    ; resxy@tmYLLabelFontAspectF=2.5
    resxy@tmXTOn=False
    resxy@tmBorderThicknessF=3.0     ;调整边框粗细
    resxy@tmXBLabelFontAspectF=1.0 
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.01        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.01  ;调整刻度的长度，注意2个命令要同时使用

    resxy@tmXBMinorLengthF  =0.005        ;調整刻度的长度
    resxy@tmXBMinorOutwardLengthF =0.005  ;调整刻度的长度，注意2个命令要同时使用
    
    resxy@tmYLLabelFontHeightF=0.01
    resxy@tmYLLabelFont=30     ;黑体
    resxy@tmYLLabelFontThicknessF=2.0

    resxy@xyLineThicknesses=(/3.0/)
    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
 
  ; resxy@pmLegendParallelPosF   = 0.45                 ; move units right
  ; resxy@pmLegendOrthogonalPosF = -0.4                ; move units down

  ; resxy@pmLegendWidthF         = 0.015                ; Change width and
  ; resxy@pmLegendHeightF        = 0.018                ; height of legend.
  ; resxy@lgLabelFontHeightF     = .03                 ; change font height
  ; resxy@lgTitleOn              = True                ; turn on legend title

  ; resxy@lgTitleFontHeightF     = .025                ; font of legend title
 
  ; resxy@lgPerimOn              = False               ; turn off box around
  ; resxy@lgLabelFontHeightF     = .008                 ; label font height

  resxy@pmLegendOrthogonalPosF = -1.28    ; Move into plot
  resxy@pmLegendParallelPosF   = 0.76      ; Move to right
  resxy@pmLegendWidthF         = 0.1       ; Decrease width
  resxy@pmLegendHeightF        = 0.08       ; Decrease height
  resxy@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
  resxy@lgLabelFontHeightF     = 0.01    ; Change the font size
  resxy@lgPerimOn=False

  resxy@xyExplicitLegendLabels = (/"Trend of vertical velocity at 500hpa"/)         ; create explicit labels
    
    delete(resxy@gsnYRefLine)

    resxy2=resxy
    delete(resxy2@gsnLeftString)
    resxy2@gsnYRefLine=0
    resxy2@gsnYRefLineDashPattern=1
    resxy2@xyLineColor      = (/"red"/)  
    resxy2@xyLineThicknesses=(/3.0/)
    resxy2@tmYRMode ="Explicit"
    resxy2@trYMaxF=10
    resxy2@trYMinF=-10
    resxy2@tmYRValues=ispan(-10,10,5)
    resxy2@tmYRLabels=ispan(-10,10,5)
    resxy2@tmYRMinorOn=False
    resxy2@tmYRLabelFontHeightF=0.01
    resxy2@tmYRLabelFont=30     ;黑体
    resxy2@tmYRLabelFontThicknessF=2.0
    
    ; resxy2@tiYAxisString ="dslp"


    resxy2@pmLegendOrthogonalPosF = -1.1    ; Move into plot
    resxy2@pmLegendParallelPosF   = 0.78      ; Move to right
    resxy2@pmLegendWidthF         = 0.1       ; Decrease width
    resxy2@pmLegendHeightF        = 0.08       ; Decrease height
    resxy2@lgPerimOn=False
    resxy2@lgBoxMinorExtentF      = 0.2      ; Shorten the legend lines
    resxy2@lgLabelFontHeightF     = 0.01     ; Change the font size
    resxy2@xyExplicitLegendLabels = (/"Mean state of vertical velocity at 500hpa"/)         ; create explicit labels




   wks=gsn_open_wks("pdf", ""+diri+"500hpa垂直速度经向图")

   plot2=gsn_csm_xy2(wks, w&longitude, w_trend, w_mean, resxy, resxy2)
draw(plot2)
frame(wks)




end 
