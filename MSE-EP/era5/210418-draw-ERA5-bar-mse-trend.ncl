
function trend(data)
begin
  rc1 =regCoef_n(ispan(1, dimsizes(data&time),1), data,0,0)
  rc1 = rc1*dimsizes(data&time)
  return rc1
end 

begin

    yearstart = 199101
    yearend = 201912

    lat1 = 0
    lat2 = 20
    lon1 = 210
    lon2 = 280

;;;;read data
    var= (/"t_mon_1x1_1979_2020.nc","q_mon_1x1_1979_2019.nc","hgt_mon_1x1_1979_2019.nc"/)

    diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

  

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
 
    t = short2flt(f1->t(timeind,{1000:100},{lat1:lat2},{lon1:lon2}))

    q = short2flt(f2->q(timeind,{1000:100},{lat1:lat2},{lon1:lon2}))
    ; q = q/1000.0  ;;;convert to kg/kg
    hgt = short2flt(f3->z(timeind,{1000:100},{lat1:lat2},{lon1:lon2}))

    
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

  A_nep_MAM = month_to_season(A_nep, "MAM")
 A_nep_JJA = month_to_season(A_nep, "JJA")
 A_nep_SON = month_to_season(A_nep, "SON")
 A_nep_DJF = month_to_season(A_nep, "DJF")

  B_nep_MAM = month_to_season(B_nep, "MAM")
 B_nep_JJA = month_to_season(B_nep, "JJA")
 B_nep_SON = month_to_season(B_nep, "SON")
 B_nep_DJF = month_to_season(B_nep, "DJF")

 C_nep_MAM = month_to_season(C_nep, "MAM")
 C_nep_JJA = month_to_season(C_nep, "JJA")
 C_nep_SON = month_to_season(C_nep, "SON")
 C_nep_DJF = month_to_season(C_nep, "DJF")

 rc_mse_mam = trend(MSE_nep_MAM)
 rc_mse_jja = trend(MSE_nep_JJA)
 rc_mse_son = trend(MSE_nep_SON)
 rc_mse_djf = trend(MSE_nep_DJF)

 rc_A_mam = trend(A_nep_MAM)
 rc_A_jja = trend(A_nep_JJA)
 rc_A_son = trend(A_nep_SON)
 rc_A_djf = trend(A_nep_DJF)
 
 rc_B_mam = trend(B_nep_MAM)
 rc_B_jja = trend(B_nep_JJA)
 rc_B_son = trend(B_nep_SON)
 rc_B_djf = trend(B_nep_DJF)

 rc_C_mam = trend(C_nep_MAM)
 rc_C_jja = trend(C_nep_JJA)
 rc_C_son = trend(C_nep_SON)
 rc_C_djf = trend(C_nep_DJF)


  ;;;

  number = new((/4,4/),float)
  number(:,0) = (/rc_mse_mam,rc_A_mam,rc_B_mam,rc_C_mam/) 
  number(:,1) = (/rc_mse_jja,rc_A_jja,rc_B_jja,rc_C_jja/) 
  number(:,2) = (/rc_mse_son,rc_A_son,rc_B_son,rc_C_son/) 
  number(:,3) = (/rc_mse_djf,rc_A_djf,rc_B_djf,rc_C_djf/) 
 

;==============================================================
; Mission 5. bar plot
; ==============================================================
  
    wks = gsn_open_wks("pdf","ERA5-MSE变化项的bar图1991-2019")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15
     sres@trXMinF = 0.5
     sres@trXMaxF = 4.5
     sres@trYMinF = 0.0
     sres@trYMaxF = 1.5


     sres@gsnDraw = True
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.15           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = (/1,2,3,4/)
     sres@tmXBLabels = (/"MAM","JJA","SON","DJF"/)
     ; sres@tmYLPrecision = 1
      ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
      ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01
     
     sres@tmXBLabelFontHeightF = 0.0205
     sres@tmXTLabelFontHeightF = 0.0205
     sres@tmYLLabelFontHeightF = 0.0225
     sres@tiMainFontHeightF = 0.025
     sres@gsnRightString = ""
     sres@tiYAxisString = "(10~S~6~N~J/kg/28yr)"
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnYRefLine = 0.      ; set the reference line equal to 0
     sres@tiMainString = "Linear trend of MSE and it's componnents (91-19) "

     str = ispan(1,4,1)
     color = (/"black","red","blue","green"/)
     plot = new(4,graphic)
     do i = 0,3
     sres@gsnXYBarChartColors = color(i)
     plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i/),number(i,:),sres)   ; draw each time series
     end do

     lbres                    = True          ; labelbar only resources
     lbres@vpWidthF           = 0.3           ; labelbar width
     lbres@vpHeightF          = 0.1           ; labelbar height
     lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
     lbres@lbFillColors       = (/"red","black"/)
     lbres@lbMonoFillPattern  = True          ; Solid fill pattern
     lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
     lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
     lbres@lbPerimOn          = False
     lbres@lgPerimColor    = "white"
     labels = (/"MSE","CpT"/)
     gsn_labelbar_ndc(wks,2,labels(::-1),0.17,0.23,lbres) ; draw right labelbar column
  
     lbres@lbFillColors       := (/"green","blue"/)
     labels := (/"Lq","gz"/)
     gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.23,lbres) ; draw left labelbar column

     frame(wks)      

end


