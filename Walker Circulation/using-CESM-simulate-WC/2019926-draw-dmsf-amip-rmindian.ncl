  ;;;************************************************************
  ;; 查看AMIP控制试验以及去掉印度洋海温增温的AMIP敏感性试验的msf表示的沃克环流强度
  ;；随时间的变化情况
  ;; linshuheng 2019.9.25
  ;;;*********************************************************
  begin 

 ;;------------------read data---------------------------------
	
	fils1  =  "/home/ys17-19/lsh/CESM-data/f19_f19_FAMIP/f19_f19_FAMIP.cam.h0.197901-200611.nc"
 
  fils2 = "/home/ys17-19/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611.nc"

  f1  =  addfile(fils1, "r")
  f2  =  addfile(fils2, "r")

  time = f2->time
  time  = (/time-1/)   ;;;

  date  =  cd_calendar(time,1)

  timeind  =  ind(date.le.200512.and.date.ge.197901)

  msf_CT  =  f1->msf_time(timeind,:,:)
  msf_rmin  =  f2->msf_time(timeind,:,:)    ;;rmin denotes rm indian warming 

  printVarSummary(msf_CT)


  msf_CT_a  =  rmMonAnnCycTLL(msf_CT)
  msf_rmin_a  =  rmMonAnnCycTLL(msf_rmin)
;;;--------------------- index of PWC intensity by msf ------------------
    ;;150°E-150°W ,all level 

    
     

     dmsf_CT = dim_avg_n_Wrap(msf_CT_a(:,:,{150:210}), 2)
     dmsf_rmin= dim_avg_n_Wrap(msf_rmin_a(:,:,{150:210}), 2)
 
     dmsf_CT_y  =  (/month_to_annual(dmsf_CT, 0)/)
     dmsf_rmin_y  =  (/month_to_annual(dmsf_rmin, 0)/)
     
     printVarSummary(dmsf_CT_y)
     time_year = ispan(1979,2005,1)
     printVarSummary(time_year)





     
;;;----------------------- draw time series --------------------------------------

	resxy = True
	resxy@gsnDraw = False
	resxy@gsnFrame = False
    
	resxy@tmXBMode  =  "Explicit"
	resxy@trXMaxF = 2005
	resxy@trXMinF = 1979
	resxy@tmXBValues = ispan(1979,2005,5)
	resxy@tmXBLabels = ispan(1979,2005,5)
	resxy@tmXBMinorOn = True
	resxy@tmXBMinorValues = ispan(1979,2005,1)

	resxy@tmYLMinorOn = False

	resxy@vpHeightF = 0.13
	resxy@vpWidthF = 0.7
	resxy@vpXF = 0.15
	resxy@vpYF = 0.95

	resxy@xyLineThicknesses  =  (/ 3.0,2.0/)          ; make second line thicker
	resxy@xyLineColors       =  (/"black","red"/)          ; change line color
	resxy@xyMarkLineModes    =  (/"lines","lines"/)
	resxy@xyExplicitLegendLabels  =  (/"AMIP CT",  "AMIP RMINDIAN/"/)     ;加线的标注

	
     wks=gsn_open_wks("pdf", "~/lsh/CESM-dmsf-Compare-amip-amiprm")

     plot = gsn_csm_xy(wks,time_year,(/dmsf_CT_y,dmsf_rmin_y/), resxy)
     draw(plot)
     frame(wks)
     end 

