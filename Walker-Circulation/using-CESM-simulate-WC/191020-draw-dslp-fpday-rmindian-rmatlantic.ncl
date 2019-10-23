  ;;;************************************************************
  ;; 查看F-pday控制试验以及去掉印度洋海温增温和大西洋增温的敏感性试验的
  ;; 用SLP表示的沃克环流强度
  ;；随时间的变化情况
  ;; linshuheng 2019.10.20
  ;;;*********************************************************
  begin 

 ;;------------------read data---------------------------------
	

    fils1 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY/F_1850-PDAY.cam.h0.197801-201212.nc"
    
    fils2 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmindian/F_1850-PDAY-rmindian.cam.h0.197801-201212.nc"
   
    fils3 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmatlantic/F_1850-PDAY-rmatlan.cam.h0.197801-201212.nc"

    fils4 = "/home/ys17-19/lsh/CESM-data/F_1850-PDAY-cam4-rmpacific/F_1850-PDAY-rmpacific.cam.h0.197801-201212.nc"


  
    f1  =  addfile(fils1, "r")
    f2  =  addfile(fils2, "r")
    f3  =  addfile(fils3, "r")
    f4  =  addfile(fils4, "r")


    date  =  cd_calendar(f1->time, 1)
    timeind  =  ind(date.le.201212.and.date.ge.197901)

    
    timestart = 197901
    timeend = 201212

    startyear_r = str_split_by_length(timestart, 4)
    endyear_r = str_split_by_length(timeend, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))


    PSL_CT  =  f1->PSL(timeind,:,:)
    PSL_rmin  =  f2->PSL(timeind,:,:)    ;;rmin denotes rm indian warming 
    PSL_rmal  =  f3->PSL(timeind,:,:)    ;;rmal denotes rm atlantic warming 
    PSL_rmpa  =  f4->PSL(timeind,:,:)    ;;rmal denotes rm atlantic warming 

     
    PSL_CT_a  =  rmMonAnnCycTLL(PSL_CT)
    PSL_rmin_a  =  rmMonAnnCycTLL(PSL_rmin)
    PSL_rmal_a  =  rmMonAnnCycTLL(PSL_rmal)
    PSL_rmpa_a  =  rmMonAnnCycTLL(PSL_rmpa)

 

;;;--------------------- index of PWC intensity ------------------
      ;;Tahiti region (5S–5N, 160–80W) and the Darwin region (5S–5N, 80–160E)

    dslp_CT  =  dim_avg_n_Wrap(PSL_CT_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_CT_a(:,{-5:5},{80:160}), (/1,2/)) 

    dslp_rmin  =  dim_avg_n_Wrap(PSL_rmin_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_rmin_a(:,{-5:5},{80:160}), (/1,2/)) 

    dslp_rmal  =  dim_avg_n_Wrap(PSL_rmal_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_rmal_a(:,{-5:5},{80:160}), (/1,2/)) 
 
    dslp_rmpa  =  dim_avg_n_Wrap(PSL_rmpa_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_rmpa_a(:,{-5:5},{80:160}), (/1,2/)) 
 

     dslp_CT_y  =  (/month_to_annual(dslp_CT, 0)/)

     dslp_rmin_y  =  (/month_to_annual(dslp_rmin, 0)/)

     dslp_rmal_y  =  (/month_to_annual(dslp_rmal, 0)/)
     
      dslp_rmpa_y  =  (/month_to_annual(dslp_rmpa, 0)/)
    

     time_year = ispan(startyear,endyear,1)
     printVarSummary(time_year)

     ;;;计算趋势

      PSL_CT_trend = regCoef_n(ispan(1,dimsizes(time_year),1),dslp_CT_y,0,0)  
      PSL_CT_trend = (/PSL_CT_trend*10*12/)  ;;; convert to decade trend 
      

      PSL_in_trend = regCoef_n(ispan(1,dimsizes(time_year),1),dslp_rmin_y,0,0)  
      PSL_in_trend = (/PSL_in_trend*10*12/)  ;;; convert to decade trend 


      PSL_al_trend = regCoef_n(ispan(1,dimsizes(time_year),1),dslp_rmal_y,0,0)  
      PSL_al_trend = (/PSL_al_trend*10*12/)  ;;; convert to decade trend 


      PSL_pa_trend = regCoef_n(ispan(1,dimsizes(time_year),1),dslp_rmpa_y,0,0)  
      PSL_pa_trend = (/PSL_pa_trend*10*12/)  ;;; convert to decade trend 
      


      print("PSL_CT_trend="+PSL_CT_trend)
      print("PSL_in_trend="+PSL_in_trend)
      print("PSL_al_trend="+PSL_al_trend)
      print("PSL_pa_trend="+PSL_pa_trend)
exit()

;;;----------------------- draw time series --------------------------------------

	resxy = True
	resxy@gsnDraw = False
	resxy@gsnFrame = False
    
	resxy@tmXBMode  =  "Explicit"
	resxy@trXMaxF = 2012
	resxy@trXMinF = 1979
	resxy@tmXBValues = ispan(1979,2012,5)
	resxy@tmXBLabels = ispan(1979,2012,5)
	resxy@tmXBMinorOn = True
	resxy@tmXBMinorValues = ispan(1979,2012,1)

	resxy@tmYLMinorOn = False

	resxy@vpHeightF = 0.4
	resxy@vpWidthF = 0.6
	resxy@vpXF = 0.15
	resxy@vpYF = 0.95

	resxy@xyLineThicknesses  =  (/ 3.0,3.0,3.0,3.0/)          ; make second line thicker
	resxy@xyLineColors       =  (/"black","red","blue","gray"/)          ; change line color
	resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines"/)
	resxy@xyExplicitLegendLabels  =  (/"AMIP CT",  "AMIP RMINDIAN","AMIP RMATLANTIC","AMIP RMPACIFIC"/)     ;加线的标注

	
  wks=gsn_open_wks("pdf", "CESM-DSLP-Compare-F-pday-rmindian-rmalt-rmpaci")

  plot = gsn_csm_xy(wks,time_year,(/dslp_CT_y,dslp_rmin_y,dslp_rmal_y,dslp_rmpa_y/), resxy)
     draw(plot)
     frame(wks)
     end 

