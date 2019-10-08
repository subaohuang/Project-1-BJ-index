  ;;;************************************************************
  ;; 查看AMIP控制试验以及去掉印度洋海温增温的AMIP敏感性试验的SLP表示的沃克环流强度
  ;；随时间的变化情况
  ;; linshuheng 2019.9.25
  ;;;*********************************************************
  begin 

 ;;------------------read data---------------------------------
	
	fils1  =  "/home/ys17-19/lsh/CESM-data/f19_f19_FAMIP/\
    f19_f19_FAMIP.cam.h0.197901-200611.nc"
 
    fils2 = "/home/ys17-19/lsh/CESM-data/f19_f19_FAMIP_rmindian/\
    f19_f19_FAMIP_rmindian.cam.h0.197901-200611.nc"
 
	f1  =  addfile(fils1, "r")
	f2  =  addfile(fils2, "r")

	date  =  cd_calendar(f1->time, 1)
	timeind  =  ind(date.le.201212.and.date.ge.197901)

	PSL_CT  =  f1->PSL(timeind,:,:)
    PSL_rmin  =  f2->PSL(timeind,:,:)    ;;rmin denotes rm indian warming 
    
    PSL_CT_a  =  rmMonAnnCycTLL(PSL_CT)
    PSL_rmin_a  =  rmMonAnnCycTLL(PSL_rmin)
;;;--------------------- index of PWC intensity ------------------
      ;;Tahiti region (5S–5N, 160–80W) and the Darwin region (5S–5N, 80–160E)

     dslp_CT  =  dim_avg_n_Wrap(PSL_CT_a(:,{-5,5},{200,280}), (/1,2/)) -\
          dim_avg_n_Wrap(PSL_CT_a(:,{-5,5},{80,160}), (/1,2/)) 
       
     dslp_rmin  =  dim_avg_n_Wrap(PSL_rmin_a(:,{-5,5},{200,280}), (/1,2/)) -\
          dim_avg_n_Wrap(PSL_rmin_a(:,{-5,5},{80,160}), (/1,2/)) 
     
     dslp_CT_y  =  (/month_to_annual(dslp_CT, 0)/)

     dslp_rmin_y  =  (/month_to_annual(dslp_rmin, 0)/)

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

	resxy@vpHeightF = 0.13
	resxy@vpWidthF = 0.7
	resxy@vpXF = 0.15
	resxy@vpYF = 0.95

	resxy@xyLineThicknesses  =  (/ 3.0,2.0/)          ; make second line thicker
	resxy@xyLineColors       =  (/"black","red"/)          ; change line color
	resxy@xyMarkLineModes    =  (/"lines","lines"/)
	resxy@xyExplicitLegendLabels  =  (/"AMIP CT",  "AMIP RMINDIAN/"/)     ;加线的标注

	
     wks=gsn_open_wks("pdf", "~/lsh/CESM-DSLP-Compare-amip-amiprm")

     plot = gsn_csm_xy(wks,time,data1 , resxy)
     draw(plot)
     frame(wks)
     end 

