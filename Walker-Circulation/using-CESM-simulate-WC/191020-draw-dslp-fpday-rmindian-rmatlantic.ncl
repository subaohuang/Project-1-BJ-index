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

  
    f1  =  addfile(fils1, "r")
    f2  =  addfile(fils2, "r")
    f3  =  addfile(fils3, "r")


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

     
    PSL_CT_a  =  rmMonAnnCycTLL(PSL_CT)
    PSL_rmin_a  =  rmMonAnnCycTLL(PSL_rmin)
    PSL_rmal_a  =  rmMonAnnCycTLL(PSL_rmal)
 

;;;--------------------- index of PWC intensity ------------------
      ;;Tahiti region (5S–5N, 160–80W) and the Darwin region (5S–5N, 80–160E)

    dslp_CT  =  dim_avg_n_Wrap(PSL_CT_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_CT_a(:,{-5:5},{80:160}), (/1,2/)) 

    dslp_rmin  =  dim_avg_n_Wrap(PSL_rmin_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_rmin_a(:,{-5:5},{80:160}), (/1,2/)) 

    dslp_rmal  =  dim_avg_n_Wrap(PSL_rmal_a(:,{-5:5},{200:280}), (/1,2/)) -\
    dim_avg_n_Wrap(PSL_rmal_a(:,{-5:5},{80:160}), (/1,2/)) 
 
     
     dslp_CT_y  =  (/month_to_annual(dslp_CT, 0)/)

     dslp_rmin_y  =  (/month_to_annual(dslp_rmin, 0)/)
     dslp_rmal_y  =  (/month_to_annual(dslp_rmal, 0)/)
    

     time_year = ispan(startyear,endyear,1)
     printVarSummary(time_year)


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

	resxy@xyLineThicknesses  =  (/ 3.0,3.0,3.0/)          ; make second line thicker
	resxy@xyLineColors       =  (/"black","red","blue"/)          ; change line color
	resxy@xyMarkLineModes    =  (/"lines","lines","lines"/)
	resxy@xyExplicitLegendLabels  =  (/"AMIP CT",  "AMIP RMINDIAN","AMIP RMATLANTIC"/)     ;加线的标注

	
  wks=gsn_open_wks("pdf", "CESM-DSLP-Compare-F-pday-rmindian-rmalt")

  plot = gsn_csm_xy(wks,time_year,(/dslp_CT_y,dslp_rmin_y,dslp_rmal_y/), resxy)
     draw(plot)
     frame(wks)
     end 

