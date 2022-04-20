load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
;;**************************************************************
;;;main program
;;***************************************************************
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
 


;;;;读入控制实验结果
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    ; fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
   sst = sst_r
    ; sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})

    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)
    
;;;读入B22000-SST海温数据

    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_IDLthree/lsh_B2000_WNP_heating_IDLthree.TEMP.fvgrid.0251-0280.nc"

    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst := f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 

    sst_sen = sst(:,{lat1:lat2},{lon1:lon2})
    sst_sen = (/rmMonAnnCycTLL(sst_sen)/)
;;;;

;;;读入B22000 WNPCLM海温数据

    ;;;;model name 
     ; modelname = "B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer"
    ;  modelname = "lsh_B2000_alt_north_year_WNP_nowindstress_new"
    ; ;;;;
    ; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    ; fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")

    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/HEAT_B_BP/HEAT_B_BP.TEMP.fvgrid.0251-0280.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"

    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    print(date)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    ; sst_r = f1->TEMP(timeind,0,:,:)  
    sst := f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    ; sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    sst_wnp = sst(:,{lat1:lat2},{lon1:lon2})
    sst_wnp = (/rmMonAnnCycTLL(sst_wnp)/)
;;;;

   nino3_ct =  dim_avg_n_Wrap(sst_ct(:,{-5:5},{210:270}), (/1,2/))   
   nino3_sen =  dim_avg_n_Wrap(sst_sen(:,{-5:5},{210:270}), (/1,2/))   
   nino3_wnp =  dim_avg_n_Wrap(sst_wnp(:,{-5:5},{210:270}), (/1,2/))   
   
   std1  = dim_stddev_n(nino3_ct, 0)
   std2  = dim_stddev_n(nino3_sen, 0)
   std3  = dim_stddev_n(nino3_wnp, 0)

   print(std1)
   print(std2)
   print(std3)
  
   
 
   x =  ispan(1, dimsizes(nino3_ct&time),1)

;;;

wks = gsn_open_wks("pdf", "Nino3_timeseries_diff_experiments")

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
    res@trXMinF = 1
    ; res@trYMinF = 0.2
    ; res@trYMaxF = 1.8

    res@tmXTOn = False
    res@tmYROn  = False

    ; res@tmXBValues=ispan(1950,2020,20)
    ; res@tmXBLabels=ispan(1950,2020,20)+""
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(1950,2020,5)
    ; res@tmYRMode="Automatic"


    plot = gsn_csm_xy(wks, x, (/nino3_ct,nino3_sen,nino3_wnp/), res)

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





