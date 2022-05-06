load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
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
;;;

;;;
    pmm_out = addfile("/home/ys17-23/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/pop.temp_upperlevel.0250-0300_levelmean.nc", "r")

    time := pmm_out->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_ctrl_hc = pmm_out->TEMP(timeind,:,:)
;;;;


year1 = 26501
year2 = 26512

;;;sst data 
  fils1 = "/home/ys17-23/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.sst.0250-0325_fvgrid.nc"
  
  f1=addfile(fils1, "r")
  time := f1->time 
  time = (/time-1/)
  date:=cd_calendar(time, 1)

  timeind:=ind(date.le.endmon.and.date.ge.startmon)
  print(date(timeind))
  sst_ctrl = f1->TEMP(timeind,0,:,:)  

  timeind2 = ind(date.le.year2.and.date.ge.year1)

  sst_265 = f1->TEMP(timeind2,0,:,:)  
;;;


;;;读入B22000 WNPCLM海温数据
    year = "0265"
    fils2 = systemfunc("ls /home/ys17-23/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.pop.h."+year+"*")
    f2 = addfiles(fils2, "r")
  
    hc_r = dim_avg_n_Wrap(f2[:]->TEMP(:,{0:31000},:,:),1)
    sst_wnp_hc = PopLatLon(hc_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    delete(hc_r)  

;;;;
    hc_ctrl_clm = clmMonTLL(sst_ctrl_hc)
    hc_diff = hc_ctrl_clm
    hc_diff = (/sst_wnp_hc - hc_ctrl_clm/)
;;;;


;;;;sst diff 
    sst_ctrl_clm = clmMonTLL(sst_ctrl)
    printVarSummary(sst_ctrl_clm)
    printVarSummary(sst_wnp)
    sst_diff = sst_ctrl_clm
    sst_diff = (/sst_wnp - sst_ctrl_clm/)
   
    sst_diff_hov = dim_avg_n_Wrap(sst_diff(:,{-5:5},{120:280}), (/1/))

    hc_diff_hov = dim_avg_n_Wrap(hc_diff(:,{-5:5},{120:280}), (/1/))

    ; SST_diff_nino3_4&time = ispan(1,12,1)

    ; SST_diff_nino3_4 = where(ismissing(SST_diff_nino3_4),0.0, SST_diff_nino3_4)


    wks = gsn_open_wks("pdf", year+"_ouc_evolution_hov_diff")
    
    plot = new(2,graphic)
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

  res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
  res@cnMinLevelValF       = -1              ; min level
  res@cnMaxLevelValF       =  1.              ; max leve
  res@cnLevelSpacingF      = 0.1

  res@cnFillPalette = "BlueWhiteOrangeRed"

  res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  ; res@tmYLLabelsOn=True
  res@tmYUseLeft=True

  res@tmYLValues    = ispan(0,11,1)
  res@tmYLLabels = (/"M","A","M","J","J","A","S","O","N","D","J(+1)","F(+1)"/)
  


  res@trYMaxF        = 12
  res@trYMinF        = 1


  res@tmXTOn=False
  res@tmYROn=False

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues= ispan(120, 280, 40)
  res@tmXBLabels = (/"120E","160E","160W","120W","80W"/)
  res@tmXBMinorOn=False
  res@gsnLeftString = "HC evolution"

  res@tmBorderThicknessF = 2.0
  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  ; res@tmXBMinorValues=ispan(0,360,10)


  plot = gsn_csm_hov(wks, SST_diff_nino3_4, res)   ; create plot
  draw(plot)
  frame(wks)

end






