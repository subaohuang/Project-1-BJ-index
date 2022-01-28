load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
;;**************************************************************
;;;main program
;;***************************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      ; simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 
begin
 ;;;time span
    startmon =25101
    endmon = 30012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;


;;;读入控制实验组的海温数据
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.sst.0251-0300_fvgrid.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)


    
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl_r = f1->TEMP(timeind,0,:,:)  

    ; sst_ctrl = PopLatLon(sst_ctrl_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    sst_ctrl = sst_ctrl_r
    ;   TLAT = f1->TLAT
    ; TLONG = f1->TLONG
   
        
;;;读入B22000 WNPCLM海温数据
    

    year = "0257"
    
    fils2 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.pop.h."+year+"*")
    f2 = addfiles(fils2, "r")
  
    ; sst_wnp = f2->TEMP(0,0,:,:)  
    sst_r = f2[:]->TEMP(:,0,:,:)  
     ;;;;差值TEMP 
    sst_wnp = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    
      
;;;;
    ; sst_ctrl_mam = month_to_season(sst_ctrl, "MAM")
    ; sst_wnp_mam = month_to_season(sst_wnp, "MAM")
    ; sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl_mam, 0)
    ; sst_wnp_clm = dim_avg_n_Wrap(sst_wnp_mam, 0)

    sst_ctrl_clm = clmMonTLL(sst_ctrl)
    sst_wnp_clm = sst_wnp


    SST_diff = sst_wnp_clm
    SST_diff = (/sst_wnp_clm - sst_ctrl_clm/)


    SST_diff_nino3_4 = dim_avg_n_Wrap(SST_diff(:,{-5:5},{210:270}), (/1,2/))
    wks = gsn_open_wks("pdf", year+"_sst_evolution_diff")
;;**************************************************************
;画图
;;***************************************************************
    res = True
    res@tmBorderThicknessF = 2.0
    res@gsnDraw        = False
    res@gsnFrame       = False

    res@gsnLeftString = ""
    res@gsnRightString = ""

    res@vpWidthF  = 0.6
    res@vpHeightF = 0.4
    ; res@vpXF      = 0.06
    ; res@vpYF      = 0.75

    res@tmXTLabelsOn   = False
    res@tmXTOn         = False
    res@tmYROn = False


    res@trXMinF        = 0
    res@trXMaxF        = 11



    res@gsnYRefLine = 0.0


    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@gsnLeftStringFontHeightF   = 0.02
    res@gsnRightStringFontHeightF  = 0.02

    res@tmXBMode = "Explicit"
    res@tmXBValues = ispan(0,11,1)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
    res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
    res@tmYLMode = "Explicit"

 
    ; res@trYMinF        = -2.0
    ; res@trYMaxF        = 2.0

    ; res@tmYLValues := fspan(-2.0, 2.0, 11)
    ; res@tmYLLabels := fspan(-2.0, 2.0, 11)
    
    res@trYMinF        = -2.0
    res@trYMaxF        = 2.0

    res@tmYLValues := fspan(-2.0, 2.0, 11)
    res@tmYLLabels := fspan(-2.0, 2.0, 11)


    res@tiYAxisString = "SST [~S~o~N~C]" 
    res@tiYAxisFontHeightF=0.02

    res@tiXAxisString = "Month"
    res@tiXAxisFontHeightF = 0.02

    res@xyMonoMarkLineMode  = True
    res@xyMarkLineMode = (/"MarkLines"/)
    res@xyMonoMarker = True
    res@xyMarker = (/1/)
    res@xyMonoMarkerSize = True
    res@xyMarkerSizeF = (/0.03/)

    res@xyMarkerColors= (/"red","blue","black","green","pink"/)


    res@xyMarkerThicknessF = 0.04
    res@xyLineThicknessF  = (/2.0/)
    res@xyLineColors       =  res@xyMarkerColors
    res@xyMonoDashPattern = True
    res@xyDashPattern     = (/1/)

    res@lgPerimOn = False
    res@lgLabelFontHeightF = 0.018

    plot  = gsn_csm_xy (wks, ispan(0,11,1), SST_diff_nino3_4, res)
    ;;添加legend 
    gres = True

    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPattern
    textres = True

    textres@lgLabels = "ensemble"+ ispan(1,5, 1)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 

    midplot = simple_legend(wks,plot,gres,lineres,textres)   
draw(plot)
frame(wks)
end






