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
    startmon =25601
    endmon = 28012

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
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl = dim_avg_n_Wrap(f1->TEMP(timeind,{0:40000},{-30:30},{120:280}),3)

    sst_ctrl&z_t = (/sst_ctrl&z_t/100.0/)
 
;;;读入B22000 WNPCLM海温数据
    
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    
    f2 = addfile(fils2, "r")
    time := f2->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)


    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    sst_wnp = dim_avg_n_Wrap(f2->TEMP(timeind,{0:40000},{-30:30},{120:280}) ,3)
     sst_wnp&z_t = (/sst_wnp&z_t/100.0/)
    
    ; sst_r = f2->TEMP(timeind,:,:)  
     ;;;;差值TEMP 
    ; sst_wnp = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    
    
;;;;
    sst_ctrl_mam = month_to_season(sst_ctrl, "MAM")
    sst_wnp_mam = month_to_season(sst_wnp, "MAM")
    sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl_mam, 0)
    sst_wnp_clm = dim_avg_n_Wrap(sst_wnp_mam, 0)
    
    ; sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl, 0)
    ; sst_wnp_clm = dim_avg_n_Wrap(sst_wnp, 0)

    SST_diff = sst_wnp_clm
    SST_diff = (/sst_wnp_clm - sst_ctrl_clm/)

    ; prob_p = t_test(sst_ctrl_mam ,sst_wnp_mam)
    ; SST_diff = where(prob_p.le.0.1, SST_diff, SST_diff@_FillValue)
    ;;;

   wks = gsn_open_wks("pdf", "Meridional-Ocean-temp-diff-WNPheating-minus-control")

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False

    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
   
    ; res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnFillPalette = "BlRe"
    res@cnInfoLabelOn  =False


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False

    ; res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.6
    res@vpHeightF = 0.4
    ; ;;
    res@gsnRightString = ""
    res@trYReverse               = True       ; reverses y-axis  


    plot= gsn_csm_contour(wks,SST_diff,res)
    
    res@cnFillOn = False 
    res@cnLinesOn = True
    res@cnLineColor = "black"
    ; res@cnMinLevelValF       =   10.0
    ; res@cnMaxLevelValF       =   28.0
    ; res@cnLevelSpacingF      =  2.0
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = 20.0

    plot_clm = gsn_csm_contour(wks,sst_ctrl_clm,res)
    res@cnLineColor = "red"
     plot_clm1 = gsn_csm_contour(wks,sst_wnp_clm,res)
    overlay(plot, plot_clm)
    overlay(plot, plot_clm1)
    ; plot(i) = gsn_csm_contour_map(wks,PRE_diff,res)


    
    ;;;打点
    ; res_pro=True
    ; res_pro@gsnAddCyclic = False
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnFillOn = True

    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLevels = (/0,0.01/)
    ; res_pro@lbLabelBarOn = False
    ; res_pro@cnLinesOn=False
    ; ; res_pro@cnMonoLineColor =False
    ; ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    ; res_pro@cnMonoFillPattern = True
    ; res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.003
  
    ; res_pro@cnMonoFillColor = False
    ; res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False
   
    ; res_pro@gsnCenterString = ""
    ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; ; res_pro@cnLineThicknessF = 3.
    ; ; res_pro@cnFillOn        = True
    ; iplot(i)=gsn_csm_contour(wks,prob_p,res_pro)


    ;  ;;;打点
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.01     ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False


    ; iplot= gsn_csm_contour(wks,prob_p,res_pro)

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; opt@gsnShadeMid     = 6         ; stipple
    ; ; opt@gsnShadeFillDotSizeF =  0.0025
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    ; iplot = gsn_contour_shade(iplot , 0,0.01,opt)
    ; overlay(plot ,iplot)
    ; print("finish subset figure")
draw(plot)
frame(wks)
end




