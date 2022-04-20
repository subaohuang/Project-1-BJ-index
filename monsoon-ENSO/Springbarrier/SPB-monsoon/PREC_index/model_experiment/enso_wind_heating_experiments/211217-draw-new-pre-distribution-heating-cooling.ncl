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
  fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/year286_enso_wind_forcing_ctrl_5m/year286_enso_wind_forcing_ctrl_5m.cam.h1.0286-0287_ensemblemean.nc")
  f1 = addfile(fils1, "r")
  time := f1->time 
  time = (/time-1/)
  date:=cd_calendar(time, 1)
  print(date)
  timeind:=ind(date.le.endmon.and.date.ge.startmon)

  td = f1->PRECC(timeind,:,:)
  td2 = f1->PRECL(timeind,:,:)

  td = (/td + td2/)

  td = (/10^3*86400*td/)  ;;;;k/day

  
  td_ctrl_24 = td
  
  dim =  dimsizes(td)
 
;;;
        
;;;read sensitive run experiment
    startmon = 28601
    endmon = 28712
;;;;
  modelname1 = "year286_enso_wind_forcing_heatsen_5m"
  ; modelname1 = "year286_enso_wind_forcing_cooltenper_5m"


;;;读入控制实验组的海温数据
  fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/"+modelname1+"/"+modelname1+".cam.h1.0286-0287_ensemblemean.nc")
  f1 = addfile(fils1, "r")
  time := f1->time 
  time = (/time-1/)
  date:=cd_calendar(time, 1)
  print(date)
  timeind:=ind(date.le.endmon.and.date.ge.startmon)

  td_heating = f1->PRECC(timeind,:,:)
  td1_heating = f1->PRECL(timeind,:,:)

  td_heating = (/td_heating + td1_heating/)

  td_heating = (/10^3*86400*td_heating/)  ;;;;k/day


  modelname2 = "year286_enso_wind_forcing_dbcoolsen_5m"

  fils2 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/"+modelname2+"/"+modelname2+".cam.h1.0286-0287_ensemblemean_21-30.nc")

  
  f2 = addfile(fils2, "r")

  td_cooling = f2->PRECC(timeind,:,:)
  td1_cooling = f2->PRECL(timeind,:,:)

  td_cooling = (/td_cooling + td1_cooling/)

  td_cooling = (/10^3*86400*td_cooling/)  ;;;;mm/day
;;;

;;;;anomaly 
  td_heating_ano = td_heating
  td_heating_ano = (/td_heating - td_ctrl_24/)

  td_cooling_ano = td_cooling
  td_cooling_ano = (/td_cooling - td_ctrl_24/)

  td_heating_ano_jjas  = dim_avg_n_Wrap(td_heating_ano(5:8,:,:), 0)
  td_cooling_ano_jjas  = dim_avg_n_Wrap(td_cooling_ano(5:8,:,:), 0)
  ; td_heating_ano_jjas = (/td_heating_ano_jjas - td_cooling_ano_jjas/)
;;;
    
; ;;; ensemble
;     ensemble_num = 10
;     ensemble = "Ensemble_" + ispan(1,ensemble_num,1)

;     fi_out1 = addfile("windctrl_pre_indian_ensemble.nc", "r")
;     PRE_sen1_nino3_4_ensemble  = fi_out1->sen1_24_pre_indian   

;     fi_out2 = addfile("heating_pre_indian_ensemble.nc", "")
;     PRE_sen2_nino3_4_ensemble =  fi_out2->sen2_24_pre_indian   

;     PRE_total1 =  new((/dimsizes(ensemble)+1,24,dim(1),dim(2)/),typeof(PRE_sen1_nino3_4_ensemble))
;     PRE_total2 =  new((/dimsizes(ensemble)+1,24,dim(1),dim(2)/),typeof(PRE_sen1_nino3_4_ensemble))

;     PRE_total1(1:,:,:,:) = PRE_sen1_nino3_4_ensemble
;     PRE_total1(0,:,:,:) = dim_avg_n_Wrap(PRE_sen1_nino3_4_ensemble, 0)

;     PRE_total2(1:,:,:,:) = PRE_sen2_nino3_4_ensemble
;     PRE_total2(0,:,:,:) = dim_avg_n_Wrap(PRE_sen2_nino3_4_ensemble, 0)

;     diff  = PRE_sen1_nino3_4_ensemble(1,1,:,:) 
;     diff =  dim_avg_n_Wrap(PRE_total1(0,0:4,:,:),0) -  dim_avg_n_Wrap(td_ctrl_24(0:4,:,:),0)
;     ; diff =  dim_avg_n_Wrap(PRE_total2(0,0:4,:,:),0) -  dim_avg_n_Wrap(td_heating(0:4,:,:),0)
; ;;;
; printMinMax(diff, 1)
; exit()
;;;;


  ; wks = gsn_open_wks("pdf", "Distribution-PREDIFF-heating-minus-windcontrol")
  ; wks = gsn_open_wks("pdf", "Distribution-PREDIFF-coolingper-minus-windcontrol")
  wks = gsn_open_wks("pdf", "Distribution-PREDIFF-dbcoolingsen-minus-windcontrol")
  ; wks = gsn_open_wks("pdf", "286year-SST-diff-minus-control")
  plot = new(2,graphic)

;;**************************************************************
;画图
;;***************************************************************
   gsn_define_colormap(wks, "NCV_blue_red")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
; 
  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 300

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.01
  res@tmYLLabelFontHeightF=0.01
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.4

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; ; res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False
  res@gsnLeftString = "Heating"
  plot(0) = gsn_csm_contour_map(wks,td_heating_ano_jjas, res)
 res@gsnLeftString = "Cooling"
  plot(1) = gsn_csm_contour_map(wks,td_cooling_ano_jjas, res) 

  pres = True

  gsn_panel(wks, plot, (/2,1/), pres)


end

