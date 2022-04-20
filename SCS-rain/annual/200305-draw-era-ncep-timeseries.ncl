begin 
;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160
 
;;;read NCEP-NCAR heat flux data

    fils2 = "/home/ys17-19/lsh/Project/SCS-rain/NCEP-NCAR-heatflux-1degree.nc"
    f2 = addfile(fils2, "r")

    date:=cd_calendar(f2->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)


    lhtfl = f2->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f2->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f2->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f2->nlwrs(timeind,{lat1:lat2},{lon1:lon2})

    lhtfl = (/-1*lhtfl/)
    shtfl = (/-1*shtfl/)
    nswrs = (/-1*nswrs/)
    nlwrs = (/-1*nlwrs/)
    

   lhtfl_a = dim_avg_n_Wrap(lhtfl(:,{0:15},{110:155}), (/1,2/))
   shtfl_a = dim_avg_n_Wrap(shtfl(:,{0:15},{110:155}), (/1,2/))
   nswrs_a = dim_avg_n_Wrap(nswrs(:,{0:15},{110:155}), (/1,2/))
   nlwrs_a = dim_avg_n_Wrap(nlwrs(:,{0:15},{110:155}), (/1,2/))
  
   

    lhtfl_annual = month_to_annual(lhtfl_a, 1)
    shtfl_annual = month_to_annual(shtfl_a, 1)
    nswrs_annual = month_to_annual(nswrs_a, 1)
    nlwrs_annual = month_to_annual(nlwrs_a, 1)


  
    
    

    total_annual = shtfl_annual
    total_annual =  nswrs_annual +nlwrs_annual + lhtfl_annual + shtfl_annual
    total_annual = total_annual
;;;;

;;;;read ERAI data
    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/ERAI-heatflux-1degree.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.197901)


    lhtfl_era = f3->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl_era = f3->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs_era = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs_era = f3->nlwrs(timeind,{lat1:lat2},{lon1:lon2})

    lhtfl_era = (/lhtfl_era/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    shtfl_era = (/shtfl_era/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nswrs_era = (/nswrs_era/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2
    nlwrs_era = (/nlwrs_era/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2

    

    lhtfl_era_a = dim_avg_n_Wrap(lhtfl_era(:,{0:15},{110:155}), (/1,2/))
    shtfl_era_a = dim_avg_n_Wrap(shtfl_era(:,{0:15},{110:155}), (/1,2/))
    nswrs_era_a = dim_avg_n_Wrap(nswrs_era(:,{0:15},{110:155}), (/1,2/))
    nlwrs_era_a = dim_avg_n_Wrap(nlwrs_era(:,{0:15},{110:155}), (/1,2/))

   
    lhtfl_annual_era = month_to_annual(lhtfl_era_a, 1)
    shtfl_annual_era = month_to_annual(shtfl_era_a, 1)
    nswrs_annual_era = month_to_annual(nswrs_era_a, 1)
    nlwrs_annual_era = month_to_annual(nlwrs_era_a, 1)
    
    total_annual_era = shtfl_annual_era
    total_annual_era = nswrs_annual_era + nlwrs_annual_era + lhtfl_annual_era + shtfl_annual_era
   
;;;;


print("***********画时间序列***************")

;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "Q-ncepandQ-era-total")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = 2013
  resxy@trXMinF = 1979
  resxy@tmXBValues = ispan(1979,2013,5)
  resxy@tmXBLabels = ispan(1979,2013,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(1979,2013,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknesses =  (/ 5,3/)          ; make second line thicker
  
  resxy@xyLineColors       =  (/"black","red"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0/)
  resxy@xyMarkLineModes    =  (/"lines","lines"/)
  
  time = ispan(1979, 2013, 1)
  plot = gsn_csm_xy(wks,time,(/total_annual,total_annual_era/), resxy)

  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       = resxy@xyLineColors
  lineres@lgLineThicknesses = (/ 5,5,5,5/)
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = resxy@xyDashPatterns


  textres = True
  textres@lgLabels = (/"total_annual","total_annual_era"/)     ;加线的标注

  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)

  end


;;;