function gain_data_an(data1)
 begin
 
 intesion_an=runave_n_Wrap(dim_rmvmean(month_to_annual(data1(0:455), 1)),3,0,0)
 intesion_an!0="time"
 intesion_an&time=ispan(1979,2016,1)


 rc1=regCoef(intesion_an&time,intesion_an)

 data_fuc=new((/3,dimsizes(intesion_an&time)/),double)


 data_fuc(0,:) = intesion_an(:)
 data_fuc(1,:) = rc1*intesion_an&time+rc1@yintercept
 data_fuc(2,:) =avg(intesion_an)

return(data_fuc)
end 
function  lead_lag(data1,data2)
begin
 mxlag    = 12
 x_Lead_y = esccr(data1,data2,mxlag)
 y_Lead_x = esccr(data2,data1,mxlag)    ; switch the order of the series
 
 ccr = new ( 2*mxlag+1, double)    
 ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
 ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"

return(ccr)

end
;---------------------

begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_iod_sst_slp.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)
sst_r=f1->sst(tind,{10:-10},:)
msl_r=f1->msl(tind,{5:-5},:)

sst=sst_r@scale_factor*sst_r+sst_r@add_offset-274.15
copy_VarMeta(sst_r, sst)
sst_ano=rmMonAnnCycTLL(sst)

msl=(msl_r@scale_factor*msl_r+msl_r@add_offset)/100.0
copy_VarMeta(msl_r, msl)

msl_ano=rmMonAnnCycTLL(msl)

dmi=dim_avg_n(sst_ano(:,:,{50:70}),(/1,2/)) - dim_avg_n(sst_ano(:,{0:-10},{90:110}),(/1,2/))
dmi!0="time"
dmi&time=ispan(1,432,1)


dslp=dim_avg_n(msl(:,:,{50:70}),(/1,2/)) - dim_avg_n(msl(:,{0:-10},{90:110}),(/1,2/))
dslp!0="time"
dslp&time=ispan(1,432,1)
printMinMax(dslp,1)


params=systemfunc("ls "+diri+"nino3.4.txt")
ncols =5
nrows=812
nhead=1
rnino3_4=readAsciiTable(params, ncols, "float", nhead)

nino3_4=rnino3_4(348:779,2)



delete(time)
delete(date)

diri1="/mnt/f/share/era/data/"
fils1=systemfunc("ls "+diri1+"era_u_sst_slp.nc")
f2=addfile(fils1, "r")
time=f2->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)

msl_r2=f2->msl(tind,{5:-5},:)

msl_p=(msl_r2@scale_factor*msl_r2+msl_r2@add_offset)/100.0
copy_VarMeta(msl_r2, msl_p)


rc=regCoef(dmi&time,dmi)
dmi_reg=dmi
dmi_reg = rc*dmi&time+rc@yintercept


dslp_p=dim_avg_n(msl_p(:,:,{195:290}),(/1,2/))-dim_avg_n(msl_p(:,:,{100:195}),(/1,2/))


;--------------------------------------------
;画DMI与DSLP的时间序列图：
;---------------------------------------
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnYRefLine=(/27/)
    resxy@gsnYRefLineThicknesses=1.5
    resxy@tmYLMode="Explicit"
    resxy@trYMaxF=1.5
    resxy@trYMinF=-1.0

    resxy@tmYLValues=fspan(-1, 1.5, 6)
    resxy@tmYLLabels=fspan(-1, 1.5, 6)

    resxy@tmYLMinorOn=False
    
    
    resxy@tmXBMode="Explicit"
    
    resxy@trXMaxF=dimsizes(dslp&time)
    resxy@tmXBValues=ispan(13,dimsizes(dslp&time),60)
    resxy@tmXBLabels=ispan(1980,2015,5)
    resxy@tmXBMinorValues=ispan(1,432,12)
    
    resxy@vpHeightF=0.4
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9

    ; resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColor      = (/"blue"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    ; resxy@xyDashPatterns    =(/0,0,2/)
    
    resxy@tiYAxisString ="DMI" 
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
    resxy@tmYLLabelFontThicknessF=2.0

    resxy@xyLineThicknesses=(/3.0/)
    resxy@pmLegendDisplayMode    = "Always"            ; turn on legend

      resxy@pmLegendOrthogonalPosF = -1.22   ; Move into plot
      resxy@pmLegendParallelPosF   = 0.858       ; Move to right
      resxy@pmLegendWidthF         = 0.1       ; Decrease width
      resxy@pmLegendHeightF        = 0.08       ; Decrease height
      resxy@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
      resxy@lgLabelFontHeightF     = 0.015     ; Change the font size
      resxy@lgPerimOn=False

      resxy@xyExplicitLegendLabels = (/"DMI"/)         ; create explicit labels
      resxy@gsnLeftString="(a)"

    resxy2=resxy
    delete(resxy2@trYMaxF)
    delete(resxy2@trYMinF)
   

    resxy2@xyLineColor      = (/"black"/)  
    resxy2@xyLineThicknesses=(/3.0/)
    resxy2@tmYRMode ="Explicit"
    resxy2@trYMaxF=3.0
    resxy2@trYMinF=-1.0
    resxy2@tmYRValues=fspan(-1.0, 3.0, 5)
    resxy2@tmYRLabels=fspan(-1.0,3.0,5)
    resxy2@tmYRMinorOn=False
    resxy2@tmYRLabelFontHeightF=0.01
    resxy2@tmYRLabelFont=30     ;黑体
    resxy2@tmYRLabelFontThicknessF=2.0
    
    resxy2@tiYAxisString ="dslp"


    resxy2@pmLegendOrthogonalPosF = -1.18 ; Move into plot
    resxy2@pmLegendParallelPosF   = 0.75      ; Move to right
    resxy2@pmLegendWidthF         = 0.1       ; Decrease width
    resxy2@pmLegendHeightF        = 0.08       ; Decrease height
    resxy2@lgPerimOn=False
    resxy2@lgBoxMinorExtentF      = 0.2      ; Shorten the legend lines
    resxy2@lgLabelFontHeightF     = 0.015     ; Change the font size
    resxy2@xyExplicitLegendLabels = (/"dslp"/)         ; create explicit labels

    



   wks2=gsn_open_wks("pdf", ""+diri+"DMI&dslp")

   plot2=gsn_csm_xy2(wks2, dslp&time, runave_n(dmi,7,0,0), runave_n(dslp,7,0,0), resxy, resxy2)

   p_res=True
   p_res@gsLineColor="red"
   ploy_line=gsn_add_polyline(wks2, plot2, dmi&time, dmi_reg, p_res)

draw(plot2)
frame(wks2)



;-----------------------------------------------------------
;画超前滞后相关
;-----------------------------------------------------------

diri3=diri+"locationxi.ascii"
diri4=diri+"locationdong.ascii"

locationxi=asciiread(diri3, 432, "float") 
locationdong=asciiread(diri4, 432, "float")

 
ccr1=lead_lag(locationxi,nino3_4)
ccr2=lead_lag(locationxi,dslp_p)
ccr3=lead_lag(locationxi,dmi)

 x=ispan(-12,12,1)     
;;; 显著性检验
  df=20
  p=0.975
  t_value=cdft_t(p, df)
  r=sqrt(t_value^2/(430+t_value^2))
;;;

 wks=gsn_open_wks("pdf", ""+diri+"dmi_slp_nino3_4")
 res=True
 res@gsnDraw=False
 res@gsnFrame=False
 res@gsnXRefLine=0
 ; res@
  
 res@gsnLeftString="(b)"
 res@vpHeightF=0.4
 res@vpWidthF=0.7
 res@xyLineColors=(/"blue","orange"/)
 res@xyMarkLineModes   = (/"lines","lines"/)
 res@xyDashPatterns    =(/0,0/)
 res@xyLineThicknesses=(/2.0,2.0/)
  res@tmXTOn=False
  res@tmYROn=False

  res@tmYLMode="Explicit"
  res@trYMaxF=1.0
  res@trYMinF=-1.0
  res@vpXF=0.15
  res@vpYF=0.9
  res@tmYLValues=fspan(-1, 1, 5)
  res@tmYLLabels=fspan(-1, 1, 5)

  res@tmYLMinorOn=False
  
  res@tmXBMode="explicit"
  res@tmXBValues=ispan(-12,12,1)
  res@tmXBLabels=ispan(-12,12,1)
  res@tmXBMinorOn=False
        
   res@tmXBLabelFontAspectF=2.5
    res@tmYLLabelFontAspectF=2.5
    res@tmXTOn=False

    res@tmBorderThicknessF=3.0     ;调整边框粗细
    res@tmXBLabelFontAspectF=1.0 
    res@tmXBLabelFontHeightF=0.008
    res@tmXBLabelFont=30     ;黑体
    res@tmXBLabelFontThicknessF=2.0
    res@tmXBMajorLengthF=0.01        ;調整刻度的长度
    res@tmXBMajorOutwardLengthF=0.01  ;调整刻度的长度，注意2个命令要同时使用

    res@tmXBMinorLengthF  =0.005        ;調整刻度的长度
    res@tmXBMinorOutwardLengthF =0.005  ;调整刻度的长度，注意2个命令要同时使用
    
    res@tmYLLabelFontHeightF=0.01
    res@tmYLLabelFont=30     ;黑体
    res@tmYLLabelFontThicknessF=2.0
     res@pmLegendDisplayMode    = "Always"            ; turn on legend
      res@pmLegendOrthogonalPosF = -1.15      ;Move into plot
      res@pmLegendParallelPosF   = 0.88      ; Move to right
      res@pmLegendWidthF         = 0.1       ; Decrease width
      res@pmLegendHeightF        = 0.08       ; Decrease height
      res@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
      res@lgLabelFontHeightF     = 0.015     ; Change the font size
      res@lgPerimOn=False
  
      res@xyExplicitLegendLabels = (/"nino3.4","dslp"/)         ; create explicit labels




;;;画bar图

 plot=gsn_csm_xy(wks,x,(/ccr1,ccr2/),res)
 
;;;;;; Create the bars and attach to appropriate plot
  gres             = True
  gres@gsFillColor = "lightgray"
  gres@tfPolyDrawOrder="PreDraw"
  xbar1=(/-1,0,0,-1,-1/)
  ybar1=(/-1,-1,1,1,-1/)
 poly= gsn_add_polygon(wks,plot,xbar1,ybar1,gres)


  xbar2=(/2,3,3,2,2/)
  ybar2=(/-1,-1,1,1,-1/)
 poly2= gsn_add_polygon(wks,plot,xbar2,ybar2,gres)

;;;画0线

   xline=ispan(-12, 12, 1)
   yline1=xline
   yline2=new(dimsizes(xline),typeof(r))
   yline3=yline2
   yline1=0
   yline2=r
   yline3=-r

  pres             = True
  pres@gsLineColor = "black"
  pres@tfPolyDrawOrder="PreDraw"
  ployline1=gsn_add_polyline(wks, plot, xline, yline1, pres)
  pres@gsLineDashPattern=1
  ployline2=gsn_add_polyline(wks, plot, xline, yline2, pres)

  ployline3=gsn_add_polyline(wks, plot, xline, yline3, pres)


draw(plot)
frame(wks)








end