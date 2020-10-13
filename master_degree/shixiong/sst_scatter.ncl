;;
begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_20cm_sst.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.195001)
sst_r=f1->sst(tind,:,:)
sst=(sst_r@scale_factor*sst_r+sst_r@add_offset)-273.15
copy_VarCoords(sst_r,sst)
delete(sst_r)
printVarSummary(sst)
sst&time=ispan(1,dimsizes(sst&time),1)
sst_pool=sst(:,{10:-10},{80:150})
sst_east=sst(:,{10:-10},{210:280})
printMinMax(sst_pool, 1)
printMinMax(sst_east, 1)

sst_pool_avg=dim_avg_n_Wrap(sst_pool, (/1,2/))
sst_east_avg=dim_avg_n_Wrap(sst_east, (/1,2/))

printMinMax(sst_pool_avg, 1)

tind2=ind(date.eq.197901)

 wks  = gsn_open_wks("pdf",diri+"scatter")             ; send graphics to PNG file
 
    res                     = True                   ; plot mods desired
    res@gsnFrame            =False
    res@gsnDraw             =False
    res@gsnMaximize         = True                   ; maximize plot in frame
    res@gsnXRefLine         =ind(date.eq.197901)+1
    res@gsnXRefLineDashPattern = 1
    res@gsnXRefLineThicknessF =1.5
    res@gsnLeftString="(b)"
    res@xyMarkLineModes     = (/"Markers","Markers"/)  ; choose which have markers
    res@xyMarkers           = (/10,11/)                     ; choose type of marker 
    res@xyMarkerColors       = (/"red","lightseagreen"/)                  ; Marker color
    res@xyMarkerSizeF       = 0.005                  ; Marker size (default 0.01)
    ; res@xyDashPatterns      = 1                      ; solid line 
    ; res@xyLineThicknesses   = (/1,2/)                ; set second line to 2

    res@tiMainString        = "Output from regline"  ; title
    res@tmXBMode="Explicit"
    
    res@trXMaxF=dimsizes(sst&time)
    res@tmXBValues=ispan(1,dimsizes(sst&time),120)
    res@tmXBLabels=ispan(1950,2015,10)
    ; res@tmXBMinorValues=ispan(1,dimsizes(sst&time),12)
    res@tmXBMinorOn=False
    res@tmBorderThicknessF=3.0     ;调整边框粗细
    res@tmYLMode="Explicit"
    res@trYMaxF=30.0

    res@tmYLValues=ispan(24,30,1)
    res@tmYLLabels=ispan(24,30,1)
    res@tmYLMinorOn=False
    res@tmXTOn=False
    res@tmYROn=False
    ;;;设置legend
    res@pmLegendDisplayMode    = "Always"            ; turn on legend
    res@pmLegendOrthogonalPosF = -1.15    ; Move into plot
    res@pmLegendParallelPosF   = 0.2      ; Move to right
    res@pmLegendWidthF         = 0.1       ; Decrease width
    res@pmLegendHeightF        = 0.08       ; Decrease height
    res@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines
    res@lgLabelFontHeightF     = 0.01    ; Change the font size
    res@lgPerimOn=False
    res@xyExplicitLegendLabels = (/"Warm pool SST","Eastern Pacific SST"/)         ; create explicit labels
  




 plot  = gsn_csm_xy (wks,sst&time,(/sst_pool_avg,sst_east_avg/),res)        ; create plot

  
;    mono = 1                            ; ascending=1 , descending=-1
;    ii   =  dim_pqsort_n(sst_pool_avg,mono,0)      ; ascending order 
;    xx  = sst&time(ii)                        
;    yy   = sst_pool_avg(ii)
   
   n = 5
   c1 = lspoly(sst&time, sst_pool_avg, 1,n)    ; all weights are set to one 
   xx=sst&time
   Y1 = c1(0) + c1(1)*xx + c1(2)*xx^2+c1(3)*xx^3+c1(4)*xx^4 
   
   c2 = lspoly(sst&time, sst_east_avg, 1,4)    ; all weights are set to one 
   
   Y2= c2(0) + c2(1)*xx + c2(2)*xx^2+c2(3)*xx^3 

  plres=True
  plres@gsLineThicknessF=5.0
  plres@gsLineDashPattern=1

 plot_line=gsn_add_polyline(wks, plot, sst&time, Y1, plres)
 
 plot_line2=gsn_add_polyline(wks, plot, sst&time, Y2, plres)

 
  c3= lspoly(sst&time(tind2:), sst_east_avg(tind2:), 1,2)    ; all weights are set to one 
   
   Y3= c3(0) + c3(1)*sst&time(tind2:)
    
  c4= lspoly(sst&time(tind2:), sst_pool_avg(tind2:), 1,2)    ; all weights are set to one 
   
   Y4= c4(0) + c4(1)*sst&time(tind2:)
   delete(plres@gsLineDashPattern)
    plres@gsLineColor="red"
    plres@gsLineThicknessF=10.0
 plot_line3=gsn_add_polyline(wks, plot, sst&time(tind2:), Y3, plres)

 plot_line4=gsn_add_polyline(wks, plot, sst&time(tind2:), Y4, plres)




draw(plot)
frame(wks)





end