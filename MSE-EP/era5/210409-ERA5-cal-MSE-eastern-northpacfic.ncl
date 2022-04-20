begin

yearstart = 199101
yearend = 201912

lat1 = 0
lat2 = 20
lon1 = 210
lon2 = 280

;;;;read data
    var= (/"t_mon_1x1_1979_2019.nc","q_mon_1x1_1979_2019.nc","hgt_mon_1x1_1979_2019.nc"/)

    diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
       
    t = short2flt(f1->t(timeind,{1000:300},{lat1:lat2},{lon1:lon2}))

    q = short2flt(f2->q(timeind,{1000:300},{lat1:lat2},{lon1:lon2}))

    ; q= q/1000.0

    hgt = short2flt(f3->z(timeind,{1000:300},{lat1:lat2},{lon1:lon2}))
;;;


Cp = 1004
Lv = 2.5*10^6
g = 1.0


;;;MSE 
 plev = t&level
 ptop = min(plev)
 psfc = 1013.0
 dz = dpres_plevel(plev, psfc, ptop, 0)
 

 A = t
 A =Cp*t

 B = q
 B = Lv*q
 
 C = hgt
 C = hgt*g
 

 MSE = t

 MSE = (/Cp*t + Lv*q + g*hgt/)


 MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

 MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

 MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))

 MSE_nep = (/MSE_nep/10^6/)


 A_col =  dim_sum_wgt_n_Wrap(A, dz, 1, 1)

 A_col = (/rmMonAnnCycTLL(A_col)/)

 A_nep = dim_avg_n_Wrap(A_col, (/1,2/))

 A_nep = (/A_nep/10^6/)

 B_col =  dim_sum_wgt_n_Wrap(B, dz, 1, 1)

 B_col = (/rmMonAnnCycTLL(B_col)/)

 B_nep = dim_avg_n_Wrap(B_col, (/1,2/))

 B_nep = (/B_nep/10^6/)

 C_col =  dim_sum_wgt_n_Wrap(C, dz, 1, 1)

 C_col = (/rmMonAnnCycTLL(C_col)/)

 C_nep = dim_avg_n_Wrap(C_col, (/1,2/))

 C_nep = (/C_nep/10^6/)

 MSE_nep_annual = month_to_annual(MSE_nep, 1)
 A_nep_annual = month_to_annual(A_nep, 1)
 B_nep_annual = month_to_annual(B_nep, 1)
 C_nep_annual = month_to_annual(C_nep, 1)

 MSE_nep_MAM = month_to_season(MSE_nep, "MAM")
 MSE_nep_JJA = month_to_season(MSE_nep, "JJA")
 MSE_nep_SON = month_to_season(MSE_nep, "SON")
 MSE_nep_DJF = month_to_season(MSE_nep, "DJF")

   plot = new(2,graphic)

   printMinMax(A_nep, 1)
   printMinMax(B_nep, 1)
   printMinMax(C_nep, 1)
   
;;;;

wks = gsn_open_wks("pdf", "MSE_anomaly_series-ERA5")

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "MSE anomalies (80W-150W,0-20N)"   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02
    res@tiYAxisString = "MSE anomaly"           ; y axis title
    res@tiXAxisString = "Year"

    res@xyLineColors      = (/"black","blue","green","red"/)  ; line colors
     res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
    res@gsnYRefLine = 0.0

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 2019
    ; res@trXMinF = 1979
    res@trXMinF = 1979
    ; res@trYMinF = 0.2
    ; res@trYMaxF = 1.8

    res@tmXTOn = False
    res@tmYROn  = False

    res@tmXBValues=ispan(1979,2020,20)
    res@tmXBLabels=ispan(1979,2020,20)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(1979,2020,5)
    res@tmYRMode="Automatic"

 plot(0) = gsn_csm_xy(wks, ispan(1979, 2019, 1), (/MSE_nep_annual,A_nep_annual,B_nep_annual,C_nep_annual/), res)


    res@xyLineColors      = (/"black","blue","green","red"/)  ; line colors
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/1.0/)        ; line thicknesses
    
    res@xyDashPatterns    = (/1.0,1.0,1.0,1.0/)        ; line patterns

    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"MSE","CpT","Lq","Hgt"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
   
    plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    plot(1) = gsn_csm_xy(wks, ispan(1979, 2019, 1), (/MSE_nep_MAM,MSE_nep_JJA,MSE_nep_SON,MSE_nep_DJF/), res)
    



    textres@lgLabels = (/"MAM","JJA","SON","DJF"/)
    res@xyLineThicknessF = (/1.0/)        ; line thicknesses
    res@xyDashPatterns    = (/1.0,1.0,1.0,1.0/)        ; line patterns



   plot(1) = simple_legend(wks,plot(1),gres,lineres,textres)   

    pres = True
    gsn_panel(wks, plot, (/1,2/), pres)

end




