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
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
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

    sst_ctrl_12 = clmMonTLL(sst_ctrl)
    
    dim = dimsizes(sst_ctrl_12)
    sst_ctrl_24 = new((/24,dim(1),dim(2)/),typeof(sst_ctrl_12))
    sst_ctrl_24(0:11,:,:) = sst_ctrl_12 
    sst_ctrl_24(12:,:,:) = sst_ctrl_12 

    sst_ctrl_24_nino3_4 = dim_avg_n_Wrap(sst_ctrl_24(:,{-5:5},{190:240}), (/1,2/))
    ;   TLAT = f1->TLAT
    ; TLONG = f1->TLONG
   
        
;;;读入B22000 WNPCLM海温数据
    
    ensemble = "Ensemble_fixinitial_" + ispan(1,20,1)
    ; ensemble = "Ensemble_" + ispan(1,10,1)

    ; startmon = 28601
    ; endmon = 28712
    ;  startmon = 28001
    ; endmon = 28112

     startmon = 25701
    endmon = 25812

    sst_wnp_nino3_4_ensemble = new((/dimsizes(ensemble),24/),typeof(sst_ctrl_24_nino3_4))
    
    ; modelname = "year257_enso_wind_forcing_heatsen"
    modelname = "year257_enso_wind_forcing_ctrl"
    ; modelname="year286_enso_wind_forcing_dbcoolsen_5m"
;     do i=0,dimsizes(ensemble)-1
 
;      diri = "/home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/year257/"+modelname+"/"
;      ; fils2=systemfunc("ls "+diri+"/"+ensemble(i)+"/"+modelname+".fixinitial.sst.0286-0287.nc")

;      ;  diri = "/WORK/sysu_hjkx_ys/array/linshh/B/enso_mon_ensemble/year286_enso_wind_forcing_dbcoolsen"
;      ; fils2=systemfunc("ls "+diri+"/"+ensemble(i)+"/"+modelname+".sst.0286-0287.nc")


;      fils2=systemfunc("ls "+diri+"/"+ensemble(i)+"/"+modelname+".pop.h.0*.nc")
;     ; print(fils2)
;     f2 = addfiles(fils2, "r")
;     ; f2 = addfile(fils2, "r")

;     time := f2[:]->time 
;     ; time := f2->time 
;     time = (/time-1/)
;     date:=cd_calendar(time, 1)
    
;     timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
;     sst_r = f2[:]->TEMP(timeind,0,:,:)
;     ; sst_r = f2->TEMP(timeind,0,:,:)  
;      ;;;差值TEMP 
;     sst_wnp = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
;     printVarSummary(sst_wnp)
;     sst_wnp_nino3_4_ensemble(i,:) = dim_avg_n_Wrap(sst_wnp(:,{-5:5},{190:240}), (/1,2/))

;     print(sst_wnp_nino3_4_ensemble(i,11))

;     end do

;     system("rm "+modelname+"_sst_nino3_4_ensemble.nc")
;     fi_out = addfile(modelname+"_sst_nino3_4_ensemble.nc", "c")
;     fi_out->sst_nino3_4_ensemble   = sst_wnp_nino3_4_ensemble
; exit()

    fi_in = addfile(modelname+"_sst_nino3_4_ensemble.nc", "r")
    ; fi_in = addfile("sst_nino3_4_ensemble.nc", "r")
    sst_wnp_nino3_4_ensemble = fi_in->sst_nino3_4_ensemble

    sst_total =  new((/dimsizes(ensemble)+1,24/),typeof(sst_ctrl_24_nino3_4))
    sst_total(1:,:) = sst_wnp_nino3_4_ensemble
    sst_total(0,:) = dim_avg_n_Wrap(sst_wnp_nino3_4_ensemble, 0)
;;;;


;;;;
    sst_ctrl_24_nino3_4_5 = conform_dims(dimsizes(sst_total), sst_ctrl_24_nino3_4, 1)
    SST_diff = (/sst_total - sst_ctrl_24_nino3_4_5/)
    print(SST_diff)
    ; exit()

   wks = gsn_open_wks("pdf", modelname+"_257year-SST-anoma")
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

    res@vpWidthF  = 0.6
    res@vpHeightF = 0.4
    ; res@vpXF      = 0.06
    ; res@vpYF      = 0.75

    res@tmXTLabelsOn   = False
    res@tmXTOn         = False
    res@tmYROn = False


    res@trXMinF        = 0
    res@trXMaxF        = 23



    res@gsnYRefLine = 0.0


    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@gsnLeftStringFontHeightF   = 0.02
    res@gsnRightStringFontHeightF  = 0.02

    res@tmXBMode = "Explicit"
    res@tmXBValues = ispan(0, 23, 1)
    ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
    res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D"/)
    res@tmYLMode = "Explicit"

 
    ; res@trYMinF        = -2.0
    ; res@trYMaxF        = 2.0

    ; res@tmYLValues := fspan(-2.0, 2.0, 11)
    ; res@tmYLLabels := fspan(-2.0, 2.0, 11)
    
    res@trYMinF        = -4.0
    res@trYMaxF        = 4.0

    res@tmYLValues := fspan(-4.0, 4.0, 11)
    res@tmYLLabels := fspan(-4.0, 4.0, 11)


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

    ; res@xyMarkerColors= (/"red","gray"/)


    res@xyMarkerThicknessF = 0.04
    res@xyLineThicknesses  = (/3.0,1.0/)
    ; res@xyLineColors       =  res@xyMarkerColors
    res@xyMonoDashPattern = True
    res@xyDashPattern     = (/1/)

    res@lgPerimOn = False
    res@lgLabelFontHeightF = 0.018

    plot  = gsn_csm_xy (wks, ispan(0,23,1), SST_diff, res)
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
end





