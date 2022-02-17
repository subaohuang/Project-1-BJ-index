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


;;读入控制实验组的海温数据
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0200-0280.nc"
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"
    
    ; fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.pop.h.0{250..300}"+"*")
    ; f1 = addfiles(fils1, "r")
    ; time := f1[:]->time 
    ; time = (/time-1/)
    ; date:=cd_calendar(time, 1)
    ; timeind:=ind(date.le.endmon.and.date.ge.startmon)
    ; ; sst_ctrl_r=  f1[:]->TEMP(timeind,{0:31000},:,:)  
    ; sst_ctrl_r=  dim_avg_n_Wrap(f1[:]->TEMP(timeind,{0:31000},:,:),1)

    ; sst_ctrl = PopLatLon(sst_ctrl_r,"gx1v6","fv1.9x2.5","bilin","da","090206")

    ; system("rm /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/pop.temp_upperlevel.0250-0300_levelmean.nc")
    ; pmm_out = addfile("/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/pop.temp_upperlevel.0250-0300_levelmean.nc", "c")
    ; pmm_out->TEMP = sst_ctrl
    ; exit()

   
    ; fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/pop.temp_upperlevel.0250-0300.nc")
    fils1 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/pop.temp_upperlevel.0250-0300_levelmean.nc")
    f1 = addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)

    date:=cd_calendar(time, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    ; sst_ctrl_r = f1->TEMP(timeind,:,:,:)  ;;;0-300m
    sst_ctrl = f1->TEMP(timeind,:,:)  ;;;0-300m
    ; sst_ctrl = dim_avg_n_Wrap(f1->TEMP(timeind,:,:,:),1)  ;;;0-300m

    ; printVarSummary(sst_ctrl_r)
    ; sst_ctrl = dim_avg_n_Wrap(sst_ctrl_r,1)  ;;;0-300
    printVarSummary(sst_ctrl)

  
        
;;;读入B22000 WNPCLM海温数据
    
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_IDLthree/lsh_B2000_WNP_heating_IDLthree.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/HEAT_B_BP/HEAT_B_BP.TEMP.fvgrid.0251-0280.nc"
    ; fils2 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"

    ; year = "0265"
    ; fils2 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.pop.h."+year+"*")
     year = "year286_windforcing"
    fils2 =  systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/year286_enso_wind_forcing_ctrl_5m/year286_enso_wind_forcing_ctrl_5m.upperlevel.temp.0286-0287_ensemblemean_1-30.nc")
    ; f2 = addfiles(fils2, "r")
    f2 = addfile(fils2, "r")
  
    ; sst_wnp = f2->TEMP(0,0,:,:)  
    ; sst_r = dim_avg_n_Wrap(f2[:]->TEMP(:,{0:31000},:,:),1)
    sst_r = dim_avg_n_Wrap(f2->TEMP(0:11,{0:31000},:,:),1)
    ; sst_r  = f2[:]->TEMP(:,{0:31000},:,:)
    
     ;;;;差值TEMP 
    sst_wnp_r = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")

    sst_wnp = sst_wnp_r  
    ; sst_wnp  = dim_avg_n_Wrap(sst_wnp_r(:,{0:31000},:,:),1)

      
;;;;
    ; sst_ctrl_mam = month_to_season(sst_ctrl, "MAM")
    ; sst_wnp_mam = month_to_season(sst_wnp, "MAM")
    ; sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl_mam, 0)
    ; sst_wnp_clm = dim_avg_n_Wrap(sst_wnp_mam, 0)

    sst_ctrl_clm = clmMonTLL(sst_ctrl)
    sst_wnp_clm = sst_wnp


    SST_diff = sst_wnp_clm
    SST_diff = (/sst_wnp_clm - sst_ctrl_clm/)


    SST_diff_nino3_4 = dim_avg_n_Wrap(SST_diff(:,{-5:5},{120:280}), (/1/))

    ; SST_diff_nino3_4 = dim_avg_n_Wrap(sst_wnp_clm(:,{-5:5},{120:280}), (/1/))

    ; SST_diff_nino3_4 = dim_avg_n_Wrap(SST_diff(:,{-5:5},:), (/1/))



    SST_diff_nino3_4&time = ispan(1,12,1)

    SST_diff_nino3_4 = where( ismissing(SST_diff_nino3_4),0.0, SST_diff_nino3_4)

    wks = gsn_open_wks("pdf", year+"_ouc_evolution_hov_diff")

;;**************************************************************
;画图
;;***************************************************************

  ; gsn_define_colormap(wks, "NCV_blue_red")

    res            = True                         ; plot mods desired
    res@gsnDraw    = False                        ; don't draw
    res@gsnFrame   = False                        ; don't advance frame
    res@gsnLeftString=""
    res@gsnRightString=""

    res@lbLabelBarOn=True
    res@lbOrientation="Vertical"

    res@cnFillOn=True
    res@cnFillDrawOrder="PreDraw"
    res@cnLinesOn=False

    res@cnInfoLabelOn=False

    
 ;    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    ; res@cnMinLevelValF       = -2.8               ; min level
    ; res@cnMaxLevelValF       =  2.8              ; max leve
    ; res@cnLevelSpacingF      =0.2

    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = -0.9               ; min level
    res@cnMaxLevelValF       =  0.9              ; max leve
    res@cnLevelSpacingF      = 0.1

    res@cnFillPalette = "BlueWhiteOrangeRed"
 
    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    ; res@tmYLLabelsOn=True
    res@tmYUseLeft=True
    ; res@trYReverse               = True       ; reverses y-axis  

    res@tmYLValues    = ispan(1,12,1)
    res@tmYLLabels    = ispan(1,12,1) +""

    res@trYMaxF        = 12
    res@trYMinF        = 1


    res@tmXTOn=False
    res@tmYROn=False

    res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
    res@tmXBValues= ispan(120, 280, 40)
    res@tmXBLabels = ispan(120, 280, 40) +"E"
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(0,360,10)

    
     plot = gsn_csm_hov(wks, SST_diff_nino3_4, res)   ; create plot
     draw(plot)
     frame(wks)

end






