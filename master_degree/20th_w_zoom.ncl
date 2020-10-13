begin
diri="/mnt/f/share/era/data/"
fils=systemfunc("ls /mnt/f/share/ncep/ncep2/omega.mon.mean.nc") 
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201012.and.date.ge.195001)



; w=f1->omega(:,{500},{30:-30},:)
; w=f1->omega(:,:,{30:-30},:)
w_r=f1->omega(:,{500},{30:-30},:)

w=w_r*w_r@scale_factor + w_r@add_offset
copy_VarCoords(w_r, w)
delete(w_r)

delete(w&time)
w&time=date

;;;选择时间段 1950-2010 ；分成1950-1979  1979-2010

tind1=ind(date.le.199412.and.date.ge.197901)
tind2=ind(date.le.201012.and.date.ge.199501)

; tind1=ind(date.le.197812.and.date.ge.195001)
; tind2=ind(date.le.201012.and.date.ge.197901)

;;;
 
printVarSummary(w)



; w_period1=dim_avg_n_Wrap(w(tind1,{850:200},:,{30:230}),(/0,1/))
; w_period2=dim_avg_n_Wrap(w(tind2,{850:200},:,{30:230}),(/0,1/))

w_period1=dim_avg_n_Wrap(w(tind1,:,{30:230}),0)
w_period2=dim_avg_n_Wrap(w(tind2,:,{30:230}),0)

printVarSummary(w_period1)
;---------------------------------------------------
;画图 
;-------------------------------------------------------
wks=gsn_open_wks("pdf", diri+"w上升区域图变化")
    
    res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="(c)"
    ; res@cnFillPalette  = "temp_diff_18lev"
  
    
    res@cnFillOn = False
    res@cnLinesOn = True
   
    res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
    res@cnLevels=-0.04
    res@cnMonoLineThickness=False
    res@cnLineThicknesses=5.0
    res@cnLineLabelsOn=True

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 230.0
    res@mpCenterLonF          = (230+30)/2

    res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"gray60","transparent"/) ;填充位置相应的具体颜色
    res@mpOutlineOn = True; turn off outline   
    res@mpNationalLineThicknessF= 2.               ; double the thickness of national boundaries 
    res@mpOutlineSpecifiers = (/"Land"/)

    res@cnLineDrawOrder="PreDraw" 
    res@cnMonoLineColor=False
    res@cnLineColors="darkgoldenrod"
    res@cnInfoLabelOn=False
    res@cnLineLabelsOn=False
    res@tmXTOn=False
    res@tmYROn=False

    ; res@tmXBMode="Explicit"
    ; res@tmXBValues=ispan(30,190,40)
    ; res@tmXBLabels=(/"30~S~o~N~E","70~S~o~N~E","110~S~o~N~E","150~S~o~N~E","170~S~o~N~W"/)
    ; res@tmXBMinorOn=False

    ; res@tmYLMode="Explicit"
    ; res@tmYLValues=ispan(-30,30,10)
   
   
    ; res@tmYLMinorOn=False
    
    res@mpShapeMode="FreeAspect"
    
    
    plot = gsn_csm_contour_map(wks, w_period1, res)

    res2=True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@cnFillOn=False
    res2@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
    res2@cnLevels=-0.04
    res2@cnLineLabelsOn=False
    res2@cnMonoLineThickness=False
    res2@cnLineThicknesses=5.0
    res2@cnMonoLineColor=False
    res2@cnInfoLabelOn=False
    res2@cnLineColors="darkgreen"
    res2@cnLineDashPattern=1
    res2@cnLineDrawOrder="PreDraw"
    res2@tmBorderThicknessF=3.0     ;调整边框粗细
    plot1 = gsn_csm_contour(wks, w_period2, res2)
    overlay(plot, plot1)
   draw(plot)
   frame(wks)

;;;******************************************************
;;;画时间序列图
;;;******************************************************


; w_time=month_to_annual(w(tind,{10:-10},{70:190}), 1)
; printVarSummary(w_time)
; w_time_avg=dim_avg_n_Wrap(w_time, (/1,2/))*(-100)
; delete(w_time)
; printMinMax(w_time_avg, 1)

; w_time_avg!0="time"
; w_time_avg&time=ispan(1,36,1)
; rc=regCoef(w_time_avg&time, w_time_avg)
; w_trend= rc*w_time_avg&time+rc@yintercept

;     resxy=True
;     resxy@gsnDraw=False
;     resxy@gsnFrame=False

;     resxy@tmXBMode="Explicit"
  
;     resxy@tmXTOn=False
;     resxy@tmYLOn=False
;     resxy@trXMaxF=36
;     resxy@tmXBValues=ispan(2,36,5)
;     resxy@tmXBLabels=ispan(1980,2010,5)
;     resxy@tmXBLabelStride=2
;     resxy@tmXBMinorOn=False
;     resxy@tmBorderThicknessF=3.0     ;调整边框粗细
;     resxy@tmYUseLeft=False
;     resxy@tmYRMode="Explicit"
;     resxy@tmYRLabelsOn         = True    ; Turn on right labels
;     resxy@trYMaxF=5
;     resxy@trYMinF=3
;     resxy@tmYRValues=fspan(3,5,5)
;     resxy@tmYRLabels=fspan(3,5,5)
;     resxy@tmYRMinorOn=False
    
;     resxy@xyLineColors      = (/"blue","red"/)          ; change line color
;     resxy@xyMarkLineModes   = (/"lines","lines"/)
;     resxy@xyDashPatterns    =(/0,0/)
;     resxy@xyLineThicknesses=(/10.0,10.0/)

;     plot_line=gsn_csm_xy(wks, w_time_avg&time, (/w_time_avg,w_trend/), resxy)

;     pres1=True
;     pres2=True


;     plot3 = gsn_attach_plots(plot,plot_line,pres1,pres2)
;     draw(plot)
;     frame(wks)


end