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
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
;;;;
begin
 ;;;time span
    startmon =25001
    endmon = 32512



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
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0200-0280.nc"
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.sst.0250-0325.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl_r= f1->TEMP(timeind,0,:,:)  

    sst_ctrl = PopLatLon(sst_ctrl_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
;;;;


;;;series

    ; ssta = rmMonAnnCycTLL(sst_ctrl)
    ssta = sst_ctrl
    sst_clm = clmMonTLL(sst_ctrl(0:(300-250)*12-1,:,:))

    ssta = calcMonAnomTLL(sst_ctrl, sst_clm)

    ssta_djf = month_to_season(ssta, "DJF")
    ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:) = ssta_djf(1:,:,:) 
    ssta_djf(dimsizes(ssta_djf&time)-1,:,:) = ssta(dimsizes(ssta&time)-1,:,:) 

    ssta_nino3_4 = dim_avg_n_Wrap(ssta_djf(:,{-5:5},{190:240}), (/1,2/))
;;;;


;;;cam data 
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/Cam/B2000_f19g16_CTRL.cam.h1.0250-0300_chazhi.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    
    PRECL_CT  =  f1->PRECL(timeind,:,:)
    PRECC_CT  =  f1->PRECC(timeind,:,:)

    PRE_CT = PRECC_CT
    PRE_CT = (/PRECC_CT+PRECL_CT/)
    ; PRE_CT = (/PRECL_CT/)  ;;;大尺度降水
    PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day
    
    prea = rmMonAnnCycTLL(PRE_CT)
    prea_indian = dim_avg_n_Wrap(prea(:,{10:30},{70:90}), (/1,2/))

    prea_indian_12  = split_12(prea_indian)
    prea_indian_jjas = dim_avg_n_Wrap(prea_indian_12(5:8,:),0)

; corr = escorc_n(prea_indian_jjas, ssta_nino3_4,0, 0)
; print(corr)
; exit()
   ; prea_indian = dim_avg_n_Wrap(prea_djf(:,{-5:5},{150:180}), (/1,2/))

wks = gsn_open_wks("pdf", "ssta_nino3_4_index")


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


    res@trXMinF        = 250
    res@trXMaxF        = 300



    res@gsnYRefLine = 0.0


    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@gsnLeftStringFontHeightF   = 0.02
    res@gsnRightStringFontHeightF  = 0.02

    res@tmXBMode = "Explicit"
    res@tmXBValues = ispan(250,300,10)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
    res@tmXBLabels = ispan(250,300,10) +""
    res@tmYLMode = "Explicit"

 
    ; res@trYMinF        = -2.0
    ; res@trYMaxF        = 2.0

    ; res@tmYLValues := fspan(-2.0, 2.0, 11)
    ; res@tmYLLabels := fspan(-2.0, 2.0, 11)
    
    res@trYMinF        = -3.0
    res@trYMaxF        = 3.0

    res@tmYLValues := ispan(-30, 30, 10)*0.1
    res@tmYLLabels := ispan(-30, 30, 10)*0.1


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
    
    x = ispan(startyear,endyear, 1)
    ; plot  = gsn_csm_xy (wks, x, (/ssta_nino3_4,prea_indian_jjas/), res)
    plot  = gsn_csm_xy (wks, x, ssta_nino3_4, res)
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







