;************************************************
; 读入数据
;************************************************
	diri="/mnt/f/share/era/data/"
	fils1=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc") 
	f1=addfile(fils1, "r")
	date=cd_calendar(f1->time, 1)
	timeind=ind(date.le.200812.and.date.ge.197901)
	sst=lonFlip(f1->sst(timeind,{20:-20},:))
	delete(sst&time)
	sst&time=date(timeind)

	ssta=rmMonAnnCycTLL(sst) 
	printVarSummary(ssta)
    ssta_aunnal=month_to_annual(ssta, 1)
;;;选择区域grad和mean看有什么变化
     
    ssta_xi=dim_avg_n_Wrap(ssta_aunnal(:,{5:-5},{130:160}), (/1,2/))
    ssta_dong=dim_avg_n_Wrap(ssta_aunnal(:,{5:-5},{170:200}), (/1,2/))
    
    ssta_gra=(/ssta_xi-ssta_dong/)


    ssta_mean=dim_avg_n_Wrap(ssta_aunnal(:,{5:-5},{120:270}), (/1,2/))

  diri_pi="/mnt/f/share/era/data/"
  wks1=gsn_open_wks("pdf", ""+diri_pi+"hadisst_gra_mean_time")
;;;画图-----------------------------------------------

  lineres=True 
  lineres@gsnDraw=False
  lineres@gsnFrame=False
  lineres@gsnYRefLine=0.0
  lineres@xyLineColors=(/"blue"/)
  lineres@xyDashPatterns=(/0,0/)
  lineres@xyLineThicknesses=3.0


  lineres@tmXTOn=False
  lineres@tmYROn=False
  lineres@tmXBMode="Explicit"

  lineres@tmXBMinorValues=ispan(1979,2008,1)
  lineres@tmXBValues=ispan(1979,2008,5)
  lineres@tmXBLabels=ispan(1979,2008,5)

  lineres@tmYLLabelFontHeightF=0.01
  lineres@tmXBLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  lineres@vpHeightF=0.4
  lineres@vpWidthF=0.6
     
  plot=new(2,graphic)


	plot(0)=gsn_csm_xy(wks1, ispan(1979,2008,1),ssta_gra,lineres)

	plot(1)=gsn_csm_xy(wks1, ispan(1979,2008,1),ssta_mean, lineres)


 pres=True

 gsn_panel(wks1, plot, (/2,1/), pres)

