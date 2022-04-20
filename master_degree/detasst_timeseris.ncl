;;set
  begin
    diri="/mnt/f/share/ersst/"
	fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
	f1=addfile(fils1, "r")
	date=cd_calendar(f1->time, 1)
	timeind=ind(date.le.201212.and.date.ge.197901)
	sst=f1->sst(timeind,{10:-10},:)
	sst!0="time"
	sst!1="lat"
	sst!2="lon"
     
    sst_annual = month_to_annual(sst, 1)
    sst_annual!0="time"
    sst_annual&time = ispan(1, 34, 1)

      
    ssta = dim_avg_n_Wrap(sst_annual(:,{5:-5},{80:160}), (/1,2/)) - \
          dim_avg_n_Wrap(sst_annual(:,{5:-5},{200:280}), (/1,2/))
      
    sst_id = dim_avg_n_Wrap(sst_annual(:,{5:-5},{40:100}), (/1,2/))

    sst_alt = dim_avg_n_Wrap(sst_annual(:,{5:-5},{330:360}), (/1,2/))
;;;
   
 

;;;;
     nino3_4=dim_avg_n(sst_annual(:,{5:-5},{190:240}), (/1,2/))  
     nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)


    ;;;remove enso cycle
     a1=regCoef_n(nino3_4, ssta, 0, 0)
     ssta_enso=a1*nino3_4+a1@yintercept
     ssta_rm = ssta
     ssta_rm = (/ssta- ssta_enso/)

    ; ssta = runave((/dim_rmvmean_n(ssta, 0)/), 3 ,1)
    ; sst_id = runave((/dim_rmvmean_n(sst_id, 0)/), 3, 1)
    ; sst_alt =runave((/dim_rmvmean_n(sst_alt, 0)/) ,3 ,1)
 
   
    ssta = (/dim_rmvmean_n(ssta, 0)/)
    sst_id = (/dim_rmvmean_n(sst_id, 0)/)
    sst_alt = (/dim_rmvmean_n(sst_alt, 0)/)
    ssta_rm = (/dim_rmvmean_n(ssta_rm, 0)/)
;;;

    
    t=ispan(1979,2012,1)

;;; trend line

     rc_pa = regCoef_n(sst_annual&time,ssta_rm, 0, 0 )

     line_pa = rc_pa*sst_annual&time + rc_pa@yintercept
   

     t1 = ispan(2008,2012,1)
     rc_pa1 = regCoef_n(t1,ssta_rm(29:33), 0, 0 )
 
     line_pa1 = rc_pa1*t1 + rc_pa1@yintercept

   

;;;

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/detaSST的时间图")

;;;画图
	res                    = True              ; plot mods desired
	res@gsnDraw            = False             ; don't draw yet
	res@gsnFrame           = False             ; don't advance frame yet

	res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
	res@vpWidthF 	         = 0.6

	res@trXMaxF            = 2012
	res@trXMinF	         = 1979              ; set x-axis minimum
	res@trYMaxF          =1.6
	res@trYMinF          =-1.6


	res@tmXTOn=False
	res@tmYROn=False
	res@gsnYRefLine     =0.0
	res@xyMonoLineColor    = False             ; want colored lines
	res@xyLineColors       = (/"green","blue","red","purple","purple"/) ; colors chosen
	res@xyLineThicknesses	 = (/4.,4.,4.,4.,2.0/)      ; line thicknesses
	res@xyDashPatterns	 = (/0,0.,0.,0.,1.0/)      ; make all lines solid
    
  



   
    plot = gsn_csm_xy(wks,t, (/sst_id,ssta,sst_alt,ssta_rm ,line_pa/), res)
    

    line_res = True
    line_res@gsLineColor = "black"
    line_res@gsLineDashPattern = 1
    line_res@gsLineThicknessF = 3.0

    addline = gsn_add_polyline(wks, plot, t1, line_pa1, line_res)

    
   ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                        ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"green","Blue","red","purple"/) ; colors chosen
    lineres@lgLineThicknesses = 3.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True
    textres@lgLabels = (/"Indian Ocean SST","Pacific SST Gradient","Atlantic Ocean SST","Pacific SST Gradient(remove enso)"/)
   
    plot= simple_legend(wks,plot,gres,lineres,textres)   

    draw(plot)
    frame(wks)

end





 