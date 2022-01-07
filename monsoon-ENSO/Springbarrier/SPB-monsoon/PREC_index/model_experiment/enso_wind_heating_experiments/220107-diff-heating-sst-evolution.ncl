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
    startmon =25001
    endmon = 30012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
;;;;


;;;读入控制实验组的海温数据
    fils1 = "/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/enso_wind_heating_experiments/pop.temp.fvgrid.0250-0300.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl = f1->TEMP(timeind,:,:) 
    
    sst_nino34 = dim_avg_n_Wrap(sst_ctrl(:,{lat1:lat2},{lon1:lon2}), (/1,2/))



    std_nino34 = dim_stddev_n_Wrap(sst_nino34, 0)
    
    print(0.5*std_nino34)
   

    sst_ctrl_12 = clmMonTLL(sst_ctrl)
    
    dim = dimsizes(sst_ctrl_12)
    sst_ctrl_24 = new((/24,dim(1),dim(2)/),typeof(sst_ctrl_12))
    sst_ctrl_24(0:11,:,:) = sst_ctrl_12 
    sst_ctrl_24(12:,:,:) = sst_ctrl_12 

    sst_ctrl_24_nino3_4 = dim_avg_n_Wrap(sst_ctrl_24(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
;;;



;;;读入控制实验组的海温数据
  ; modelname  = "year286_enso_wind_forcing_coolingsen_5m"
  ; modelname  = "year286_enso_wind_forcing_dbcoolsen_5m"
    modelname  = "year286_enso_wind_forcing_ctrl_5m"
    fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/"+modelname+"/"+modelname+".sst.0286-0287_ensemblemean.nc")

    f1 = addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    print(date)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    td = f1->TEMP(timeind,0,:,:)

    sst_wind_24 = td

    sst_wind = PopLatLon(sst_wind_24,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_wind)
;;;;


;;;;
    modelname  = "year286_enso_wind_forcing_heatsen_5m"
    ; modelname  = "year286_enso_wind_forcing_dbheatsen_5m"
    fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/"+modelname+"/"+modelname+".sst.0286-0287_ensemblemean.nc")

    f1 = addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    print(date)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    td = f1->TEMP(timeind,0,:,:)

    sst_heat_24 = td

    sst_heat = PopLatLon(sst_heat_24,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_heat)
;;;;


;;;;
    modelname  = "year286_enso_wind_forcing_dbcoolsen_5m"
    fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/"+modelname+"/"+modelname+".sst.0286-0287_ensemblemean.nc")
    ; fils1 = systemfunc("ls /WORK/sysu_hjkx_ys/linshh/B/ensemble_run/"+modelname+"/ensemble/"+modelname+".sst.0286-0287_ensemblemean_11-20.nc")

    f1 = addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    print(date)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    td = f1->TEMP(timeind,0,:,:)
    sst_cool_24 = td
    sst_cool = PopLatLon(sst_cool_24,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_cool)
;;;;
   
;;;nino3.4

    sst_wind_nino34 = dim_avg_n_Wrap(sst_wind(:,{lat1:lat2},{lon1:lon2}), (/1,2/))

    sst_heat_nino34 = dim_avg_n_Wrap(sst_heat(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    
    sst_cool_nino34 = dim_avg_n_Wrap(sst_cool(:,{lat1:lat2},{lon1:lon2}), (/1,2/))



    sst_wind_nino34 = (/sst_wind_nino34 - sst_ctrl_24_nino3_4/)

    sst_heat_nino34 = (/sst_heat_nino34 - sst_ctrl_24_nino3_4/)  

    sst_cool_nino34 = (/sst_cool_nino34 - sst_ctrl_24_nino3_4/)    
    
   plot = new(3,graphic)
    
   ; wks = gsn_open_wks("pdf", "dbheating-SST-evolution-286year-SST-diff-5mWBforcing-minus-control")
   wks = gsn_open_wks("pdf", "SST-evolution-286year-SST-diff-5mWBforcing-minus-control")
      ; wks = gsn_open_wks("pdf", "286year-SST-diff-minus-control")

 
;;**************************************************************
;画图
;;***************************************************************
    res = True
    res@tmBorderThicknessF = 2.0
    res@gsnDraw        = False
    res@gsnFrame       = False

    res@gsnLeftString = ""
    res@gsnRightString = ""

    res@vpWidthF  = 0.7
    res@vpHeightF = 0.3
    ; res@vpXF      = 0.06
    ; res@vpYF      = 0.75

    res@tmXTLabelsOn   = False
    res@tmXTOn         = False
    res@tmYROn = False


    res@trXMinF        = 0
    res@trXMaxF        = 19



    res@gsnYRefLine = 0.0


    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@gsnLeftStringFontHeightF   = 0.02
    res@gsnRightStringFontHeightF  = 0.02

    res@tmXBMode = "Explicit"
    ; res@tmXBValues = ispan(0, 17, 1)
    res@tmXBValues = ispan(0, 19, 1)
    res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
    res@tmYLMode = "Explicit"

 
    ; res@trYMinF        = -2.0
    ; res@trYMaxF        = 2.0

    ; res@tmYLValues := fspan(-2.0, 2.0, 11)
    ; res@tmYLLabels := fspan(-2.0, 2.0, 11)
    
    res@trYMinF        = -0.6
    res@trYMaxF        = 1.6

    res@tmYLValues := ispan(-6,16,2)*0.1
    res@tmYLLabels := ispan(-6,16,2)*0.1


    res@tiYAxisString = "SST [~S~o~N~C]" 
    res@tiYAxisFontHeightF=0.02

    res@tiXAxisString = "Month"
    res@tiXAxisFontHeightF = 0.02

    res@xyMonoMarkLineMode  = True
    res@xyMarkLineMode = (/"Lines"/)
    res@xyMonoMarker = True
    res@xyMarker = (/1/)
    res@xyMonoMarkerSize = True
    res@xyMarkerSizeF = (/0.03/)

    res@xyMarkerColors= (/"black","red","blue"/)


    res@xyMarkerThicknessF = 0.04
    res@xyLineThicknesses  = (/3.0,1.0,1.0/)
    res@xyLineColors       =  res@xyMarkerColors
    res@xyMonoDashPattern = True
    res@xyDashPattern     = (/1/)

    res@lgPerimOn = False
    res@lgLabelFontHeightF = 0.018
    

    
    ; plot(0)  = gsn_csm_xy (wks, ispan(0,17,1), (/sst_wind_nino34(2:19),sst_heat_nino34(2:19),sst_cool_nino34(2:19)/), res)
    plot(0)  = gsn_csm_xy (wks, ispan(0,19,1), (/sst_wind_nino34(0:19),sst_heat_nino34(0:19),sst_cool_nino34(0:19)/), res)

     
    ; ;;添加legend 
    ; gres = True

    ; gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;     ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ; gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ; lineres = True
    ; ; lineres@lgLineColors = (/"black"，/) ; line colors
    ; lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    ; lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    ; lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    ; lineres@lgDashIndexes = res@xyDashPattern
    ; textres = True

    ; textres@lgLabels = "ensemble"+ ispan(1,5, 1)
    ; textres@lgLabelFontHeightF = 0.01
    ; ; textres@lgLabelOpacity = (/0/) 

    ; midplot = simple_legend(wks,plot,gres,lineres,textres)   

draw(plot)
frame(wks)
  ; resp = True

  ; gsn_panel(wks, plot, (/3,1/), resp)

end





