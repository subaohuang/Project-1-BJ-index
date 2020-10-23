
;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201612.and.date.ge.190001)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ;;;
;;;;

;;;;计算标准偏差(5S-5N)

    ssta_mean = dim_avg_n_Wrap(ssta(:,{-5:5},{120:280}), 1)

    std1 = dim_stddev_n(ssta_mean(0:49,:), 0)

    std2 = dim_stddev_n(ssta_mean(59:,:), 0)

    x = ssta_mean&longitude

;;;;
data_all = (/std1,std2/)

;;;;draw 
;;*************************************************
; create plot
;**************************************************
    wks = gsn_open_wks("pdf","STD-HadISST")         ; send graphics to PNG file

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "Standard deviations of SSTA in HadISST"   ; title
    res@tiYAxisString = "SST anomaly s.d."           ; y axis title
   
    res@tiMainFontHeightF = 0.015
    res@xyLineColors      = (/"blue","red"/)  ; line colors
    res@xyLineThicknesses = (/4.0,4.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0/)        ; line patterns

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 280.0
    res@trXMinF = 120.0

    res@trYMinF = 0.0
    res@trYMaxF = 1.0

    


    res@tmXTOn = False
    res@tmYROn  = False
    ; res@tmYLOn  = False

    ; res@tmXBValues=(/0,90,180,270,357.5/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(0,360,10)
    res@vpWidthF = 0.6
    res@vpHeightF = 0.4

    res@tmXBValues=ispan(135,270,45)
    res@tmXBLabels=(/"135~S~o~N~E","180~S~o","135~S~o~N~W","90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(120,280,15)
    res@tmYRMode="Automatic"




  plot = gsn_csm_xy(wks,x,data_all,res)          ; Draw an XY plot with 1 curve.

         ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"Blue","red"/) ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"1900-1950","1960-2014"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 


    mid = simple_legend(wks,plot,gres,lineres,textres)   

draw(plot)
frame(wks)
  
end
