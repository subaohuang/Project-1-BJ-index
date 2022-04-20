begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_iod_pacific_w.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind1=ind(date.le.201412.and.date.ge.197901)
tind2=ind(date.le.201501.and.date.ge.197901)
w_r=f1->w(:,:,:)

w=w_r@scale_factor*w_r+w_r@add_offset
copy_VarMeta(w_r, w)
w_i=w(tind1,{20:-20},{60:105})
w_p=w(tind2,{20:-20},{110:180})

w_i_in=dim_avg_n_Wrap(w_i, (/1,2/))
w_p_in=dim_avg_n_Wrap(w_p, (/1,2/))
printVarSummary(w_i_in)
printVarSummary(w_p_in)

;;;选择印度上升支强度定义，挑选出7月到11月作为JASON
jason=new(180,typeof(w_i_in))
do i=0,35
jason(i*5:4+i*5)=w_i_in((6+i*12):(10+i*12))*(-100)
end do
printVarSummary(jason)

;;;选择太平洋上升支强度定义，挑选出9月到1月作为sondj
sondj=new(180,typeof(w_i_in))
do i=0,35
sondj(i*5:4+i*5)=w_p_in((8+i*12):(12+i*12))*(-100)
end do

printVarSummary(sondj)
sondj!0="time"
sondj&time=ispan(1,180,1)

printMinMax(sondj, 1)
printMinMax(jason, 1)

wks2= gsn_open_wks("pdf", ""+diri+"上升强度对比")
;*******************************************************
;画图
;********************************************************
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnYRefLine=(/27/)
    resxy@gsnYRefLineThicknesses=1.5
    resxy@tmYLMode="Explicit"
    resxy@trYMaxF=4
    resxy@trYMinF=0
    resxy@tmYLValues=ispan(0,4,1)
    resxy@tmYLLabels=ispan(0,4,1)
    resxy@tmYLMinorOn=False
    
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(sondj&time)
    resxy@tmXBValues=ispan(6,dimsizes(sondj&time),25)
    resxy@tmXBLabels=ispan(1980,2015,5)
    resxy@tmXBMinorValues=ispan(1,180,5)
    
    resxy@vpHeightF=0.3
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9

    ; resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColor      = (/"black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    ; resxy@xyDashPatterns    =(/0,0,2/)
    
    resxy@tiYAxisString ="SONDJ" 
    resxy@tmXBLabelFontAspectF=2.5
    resxy@tmYLLabelFontAspectF=2.5
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
    resxy@tmYLLabelFontThicknessF=3.0

    resxy@xyLineThicknesses=(/3.0/)

      resxy@pmLegendDisplayMode    = "Always"            ; turn on legend
      resxy@pmLegendOrthogonalPosF = -1.28    ; Move into plot
      resxy@pmLegendParallelPosF   = 0.75       ; Move to right
      resxy@pmLegendWidthF         = 0.1       ; Decrease width
      resxy@pmLegendHeightF        = 0.08       ; Decrease height
      resxy@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
      resxy@lgLabelFontHeightF     = 0.015     ; Change the font size
      resxy@lgPerimOn=False

    resxy@xyExplicitLegendLabels = (/"SONDJ"/)         ; create explicit labels
    
   
    resxy2=resxy
    resxy2@xyLineColor      = (/"blue"/)  
    resxy2@xyLineThicknesses=(/3.0/)
    resxy2@tmYRMode ="Explicit"
    resxy2@trYMinF=0
    resxy2@trYMaxF=5
    resxy2@tmYRValues=ispan(0,5,1)
    resxy2@tmYRLabels=ispan(0,5,1)
    resxy2@tmYRMinorOn=False

    resxy2@tmYRLabelFontHeightF=0.01
    resxy2@tmYRLabelFont=30     ;黑体
    resxy2@tmYRLabelFontThicknessF=2.0
    
    resxy2@tiYAxisString ="JASON"
    

    resxy2@pmLegendOrthogonalPosF = -1.23   ; Move into plot
    resxy2@pmLegendParallelPosF   = 0.50     ; Move to right
    resxy2@pmLegendWidthF         = 0.1       ; Decrease width
    resxy2@pmLegendHeightF        = 0.08       ; Decrease height
    resxy2@lgPerimOn=False
    resxy2@lgBoxMinorExtentF      = 0.2      ; Shorten the legend lines
    resxy2@lgLabelFontHeightF     = 0.015     ; Change the font size
    resxy2@xyExplicitLegendLabels = (/"JASON"/)         ; create explicit labels
   


plot=gsn_csm_xy2(wks2, sondj&time, runave_n(sondj,3,0,0), runave_n(jason,3,0,0), resxy, resxy2)
   
draw(plot)
frame(wks2)

end
;;；

