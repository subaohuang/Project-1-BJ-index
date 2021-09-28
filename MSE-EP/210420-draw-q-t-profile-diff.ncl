begin

yearstart = 197901
yearend = 201912
yearstart = 197901
yearend = 201912
year1 = 1979
yearmid = 1999
year2 = 2019

lat1 = 0
lat2 = 20
lon1 = 210
lon2 = 280

;;;;read data
    var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(ind(date.eq.yearstart)))
    delete(date)     
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q = q/1000.0  ;;;convert to kg/kg

    hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;

    Cp = 1004
    Lv = 2.5*10^6
    g = 9.8
    
    t_a_r = rmMonAnnCycLLLT(t(level|:,lat|:,lon|:,time|:))
    t_a = t_a_r(time|:,level|:,lat|:,lon|:)
    
    q_a_r = rmMonAnnCycLLLT(q(level|:,lat|:,lon|:,time|:))
    q_a = q_a_r(time|:,level|:,lat|:,lon|:)

    hgt_a_r = rmMonAnnCycLLLT(hgt(level|:,lat|:,lon|:,time|:))
    hgt_a = hgt_a_r(time|:,level|:,lat|:,lon|:)

    t = t_a
    q= q_a
    hgt = hgt_a


    ;;;MSE 
    plev = t&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)
    ;;;

    t  = Cp*t
    q  = Lv*q
    
    t = t/10^3
    q = q/10^3
    hgt = hgt/10^3


    season  = (/"MAM","JJA","SON","DJF"/)
    order = (/"(a)","(b)","(c)","(d)"/)

    plot = new(4,graphic)

    plot_q = new(4,graphic)
    wks = gsn_open_wks("pdf", "00-19 minus 79-99diff-MSE-CPT-Lvq_change")

do i=0,3
   

 t_MAM = month_to_season(t, season(i))
 q_MAM = month_to_season(q, season(i))
  hgt_MAM = month_to_season(hgt, season(i))



 t_nep_MAM = dim_avg_n_Wrap(t_MAM, (/2,3/))

 q_nep_MAM = dim_avg_n_Wrap(q_MAM, (/2,3/))

 hgt_nep_MAM = dim_avg_n_Wrap(hgt_MAM, (/2,3/))

 t_nep_MAM&time = ispan(year1, year2, 1)
 q_nep_MAM&time = ispan(year1, year2, 1)
 hgt_nep_MAM&time = ispan(year1, year2, 1)



 t_nep_clm_mam = dim_avg_n_Wrap(t_nep_MAM, 0)


 t_nep_mam_change = t_nep_MAM(1,:)
 t_nep_mam_change = dim_avg_n_Wrap(t_nep_MAM({yearmid:year2},:), 0) - dim_avg_n_Wrap(t_nep_MAM({year1:yearmid},:), 0)

q_nep_mam_change = q_nep_MAM(1,:)
q_nep_mam_change = dim_avg_n_Wrap(q_nep_MAM({yearmid:year2},:), 0) - dim_avg_n_Wrap(q_nep_MAM({year1:yearmid},:), 0)


 hgt_nep_mam_change = hgt_nep_MAM(1,:)
 hgt_nep_mam_change = dim_avg_n_Wrap(hgt_nep_MAM({yearmid:year2},:), 0) - dim_avg_n_Wrap(hgt_nep_MAM({year1:yearmid},:), 0)
 
 mse_change = hgt_nep_mam_change + t_nep_mam_change + q_nep_mam_change

   

   
;;;;
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/2.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       :=  (/"red","blue","green","black"/)          ; change line color
    ; resxy@xylineColors :=  (/"red","blue","green","black"/)
    resxy@xyMarkLineModes    :=  (/"MarkLines","MarkLines","MarkLines","MarkLines"/)
    resxy@xyDashPatterns := (/0,0,0,0/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           := (/16,4,10,15/)
    resxy@xyMarkerSizes = (/0.01,0.01,0.01,0.01/)
    ; resxy@xyMarkerColors :=  (/"blue","blue","red","red"/)
    resxy@xyMarkerColors :=  (/"red","blue","green","black"/)

    resxy@gsnLeftStringFontHeightF = 0.02
    resxy@gsnLeftString = order(i)+" " +season(i)

  

    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
    resxy@tmXBMajorThicknessF=1.0
    resxy@tmYLMajorThicknessF=1.0

    
    resxy@tmXBMode  =  "Explicit"
    resxy@trXMaxF =  3.
    resxy@trXMinF = -1.
    resxy@tmXBValues = ispan(-10,30,5)*0.1
    resxy@tmXBLabels = ispan(-10,30,5)*0.1
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = ispan(-10,30,1)*0.1
    resxy@tiYAxisString = ""
   
    ; resxy@pmLegendDisplayMode = "Always"
    ; resxy@lgBoxLinesOn = False
    ; resxy@xyExplicitLegendLabels  =  (/"DPLH Clm","DPLH Clm+Change","DIABH Clm","DIABH Clm+Change"/)     ;加线的标注
    resxy@trYReverse        = True                       ; reverse Y-axis

    ; plot(i) = gsn_csm_xy(wks, (/t_nep_clm_mam,t_nep_mam_change, \
    ;                          q_nep_clm_mam,q_nep_mam_change/),plev ,resxy)

    plot(i) = gsn_csm_xy(wks, (/t_nep_mam_change,q_nep_mam_change, \
                             hgt_nep_mam_change,mse_change/),plev ,resxy)
    
    ; ;;;
    ; resxy@xyLineColors       :=  (/"black","black"/)          ; change line color
    ; resxy@xyMarkLineModes := (/"lines","lines"/)    
    ; ; resxy@xyExplicitLegendLabels := (/"Q1 Clim","Q1 Clim+Change"/)
    ; resxy@xyDashPatterns := (/0,1/)
    ; plot2(i) = gsn_csm_xy(wks, (/q1_clim_area,q1_change/),q1&lev, resxy)

    ; overlay(plot1(i), plot2(i))
      ;;添加legend 
    gres = True
    gres@YPosPercent = 20.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 80    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = resxy@xyDashPatterns
    lineres@lgMonoMarkerIndex = False 
    lineres@lgMarkerIndexes = resxy@xyMarkers
    lineres@lgItemTypes =  resxy@xyMarkLineModes
    lineres@lgMarkerIndex = resxy@xyMarkers
                        ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"CpT","Lq","hgt","MSE"/)
    textres@lgLabelFontHeightF = 0.01
    plot(i) = simple_legend(wks,plot(i),gres,lineres,textres)   
end do 




    pres = True
    gsn_panel(wks, plot, (/2,2/), pres)

end




