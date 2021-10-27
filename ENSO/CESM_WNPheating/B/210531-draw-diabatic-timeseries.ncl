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
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_ctrl = f1->PTENDT(timeind,{1000:100},:,:)
    sst_ctrl_avg = dim_avg_n_Wrap(sst_ctrl, 1)

    sst_ctrl_spring = month_to_season(sst_ctrl, "MAM")


    ; sst_ctrl_DTCOND = f1->DTCOND(timeind,{1000:100},:,:)  ;;;k/s
    ; sst_ctrl_DTCOND_avg = dim_avg_n_Wrap(sst_ctrl_DTCOND, 1)
    ; sst_ctrl_DTCOND_mam = month_to_season(sst_ctrl_DTCOND_avg, "MAM")
    ; sst_ctrl_DTCOND_clm = dim_avg_n_Wrap(sst_ctrl_DTCOND_mam, 0)
 
;;;读入B22000 WNPCLM海温数据
    
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.PTENDT.cam.h1.0251-0280_chazhi.nc"
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils2 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    
    f2 = addfile(fils2, "r")
    time := f2->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    sst_wnp = f2->PTENDT(timeind,{1000:100},:,:)  
    ; sst_wnp = f2->DTCOND(timeind,{1000:100},:,:)  
    sst_wnp_avg = dim_avg_n_Wrap(sst_wnp, 1)
;;;;

;;读入B22000 WNPCLM海温数据
    
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils2 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    
    f2 = addfile(fils2, "r")
    time := f2->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    sst_idl = f2->PTENDT(timeind,{1000:100},:,:)  
    ; sst_wnp = f2->DTCOND(timeind,{1000:100},:,:)  
    sst_idl_avg = dim_avg_n_Wrap(sst_idl, 1)


    sst_ctrl_mam = month_to_season(sst_ctrl_avg, "MAM")
    sst_wnp_mam = month_to_season(sst_wnp_avg, "MAM")
    sst_idl_mam = month_to_season(sst_idl_avg, "MAM")


    sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl_mam(:,{0:20},{125:160}), (/1,2/))
    sst_wnp_clm = dim_avg_n_Wrap(sst_wnp_mam(:,{0:20},{125:160}), (/1,2/))
    sst_idl_clm = dim_avg_n_Wrap(sst_idl_mam(:,{0:20},{125:160}), (/1,2/))

    sst_ctrl_clm = (/86400*sst_ctrl_clm/1004.0/)
    sst_wnp_clm = (/86400*sst_wnp_clm/1004.0/)
    sst_idl_clm = (/86400*sst_idl_clm/1004.0/)
    ; SST_diff = (/86400*SST_diff/1004.0/)

 std1 = dim_stddev_n(sst_ctrl_clm, 0)
 std2 = dim_stddev_n(sst_wnp_clm, 0)
 std3 = dim_stddev_n(sst_idl_clm, 0)
 
 print(std1)
 print(std2)
 print(std3)
 exit()
;;; get variables: temperature, longitude, latitude

x = ispan(251,280,1)

  
  wks = gsn_open_wks("pdf", "diabatic_timeseries_diff_experiments")

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "Annual mean MSE anomalies (80W-150W,0-20N)"   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02
    res@tiYAxisString = "MSE anomaly"           ; y axis title
    res@tiXAxisString = "Year"

    ; res@xyLineColors      = (/"black","red","blue","green"/)  ; line colors
    res@xyLineColors      = (/"black","red","blue"/)  ; line colors
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/1.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0,0.0/)        ; line patterns
    res@gsnYRefLine = 0.0

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = max(x)
    ; res@trXMinF = 1950
    res@trXMinF = min(x)
    ; res@trYMinF = 0.2
    ; res@trYMaxF = 1.8

    res@tmXTOn = False
    res@tmYROn  = False

    ; res@tmXBValues=ispan(1950,2020,20)
    ; res@tmXBLabels=ispan(1950,2020,20)+""
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(1950,2020,5)
    ; res@tmYRMode="Automatic"


    plot = gsn_csm_xy(wks, x, (/sst_ctrl_clm,sst_wnp_clm,sst_idl_clm/), res)

    ; ;;添加legend 
    ; gres = True
    ; gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;       ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ; gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ; lineres = True
    ; ; lineres@lgLineColors = (/"black"，/) ; line colors
    ; lineres@lgLineColors       = res@xyLineColors ; colors chosen
    ; lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    ; lineres@LineLengthPercent = 9.
    ; lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    ; textres = True

    ; textres@lgLabels = (/"MSE","CpT","Lq","Hgt"/)
    ; textres@lgLabelFontHeightF = 0.01
    ; ; textres@lgLabelOpacity = (/0/) 
   
    ; plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    ; plot(1) = gsn_csm_xy(wks, ispan(1950, 2019, 1), (/MSE_nep_MAM,MSE_nep_JJA,MSE_nep_SON,MSE_nep_DJF/), res)
    



    ; textres@lgLabels = (/"MAM","JJA","SON","DJF"/)
    ; res@xyLineThicknessF = (/1.0/)        ; line thicknesses
    ; res@xyDashPatterns    = (/1.0,1.0,1.0,1.0/)        ; line patterns



   ; plot(1) = simple_legend(wks,plot(1),gres,lineres,textres)   

    ; pres = True
    ; gsn_panel(wks, plot, (/1,2/), pres)
    draw(plot)
    frame(wks)

end



