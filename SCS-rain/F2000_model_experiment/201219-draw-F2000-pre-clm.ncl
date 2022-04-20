;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 
begin
    fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012_chazhi.nc"
 
    f1=addfile(fils1, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.601.and.date.le.4012)

    PRECL_CT  =  f1->PRECL(timeind,:,:)
    PRECC_CT  =  f1->PRECC(timeind,:,:)

    PRE_CT = PRECC_CT
    PRE_CT = (/PRECC_CT+PRECL_CT/)
    ; PRE_CT = (/PRECL_CC/)  ;;;大尺度降水
    PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day


    pre_mam = dim_avg_n_Wrap(month_to_season(PRE_CT, "MAM"),0)
    pre_jja = dim_avg_n_Wrap(month_to_season(PRE_CT, "JJA"),0)
    pre_son = dim_avg_n_Wrap(month_to_season(PRE_CT, "SON"),0)
    pre_djf = dim_avg_n_Wrap( month_to_season(PRE_CT, "DJF"),0)


    
    wks=gsn_open_wks("pdf", "CESM-F2000-降水的气候态")
    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =True
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette  = "MPL_Blues"




 


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   0.0
    res@cnMaxLevelValF       =   1
    res@cnLevelSpacingF      =  0.1
    

   

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    ; res@mpMinLatF             = -30.0
    ; res@mpMaxLatF             = 30.0
    ; res@mpMinLonF             = 30.0
    ; res@mpMaxLonF             = 160.0

    res@mpMinLatF             = -40.0
    res@mpMaxLatF             = 40.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0




    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5


    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3


    plot = new(4,graphic)
    plot(0) = gsn_csm_contour_map(wks, pre_mam, res)
    plot(1) = gsn_csm_contour_map(wks, pre_jja, res)
    plot(2) = gsn_csm_contour_map(wks, pre_son, res)  
    plot(3) = gsn_csm_contour_map(wks, pre_djf, res)

 
    pres = True
    pres@gsnPanelLabelBar = True
    gsn_panel(wks, plot, (/2,2/), pres)

    
    print("done script")
  end
  

