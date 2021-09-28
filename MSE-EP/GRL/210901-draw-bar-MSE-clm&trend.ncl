
function trend(data)
begin
  rc1 =regCoef_n(ispan(1, dimsizes(data&time),1), data,0,0)
  rc1 = rc1*dimsizes(data&time)
  return rc1
end 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 

begin

    yearstart = 198001
    yearend = 202012

    lat1 = 0
    lat2 = 20
    lon1 = 250
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
   
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q = q/1000.0  ;;;convert to kg/kg

    hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;

Cp = 1004
Lv = 2.5*10^6
g = 9.8


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

 MSE_ne = dim_avg_n_Wrap(MSE_col, (/1,2/))

 MSE_ne_12 = split_12(MSE_ne)

 MSE_clm = avg(MSE_ne_12(0:4,:))

 MSE_clm = (/MSE_clm/10^6/) ;;;气候态


 MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

 MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))

 MSE_nep = (/MSE_nep/10^6/)
;;;;




;;;;;;;;;;;;;;;;;;;
 A_col =  dim_sum_wgt_n_Wrap(A, dz, 1, 1)
 
 A_ne = dim_avg_n_Wrap(A_col, (/1,2/))
 A_ne_12 = split_12(A_ne)
 A_clm = avg(A_ne_12(0:4,:))

 A_clm = (/A_clm/10^6/)

 A_col = (/rmMonAnnCycTLL(A_col)/)

 A_nep = dim_avg_n_Wrap(A_col, (/1,2/))

 A_nep = (/A_nep/10^6/)
 ;;;;;


 B_col =  dim_sum_wgt_n_Wrap(B, dz, 1, 1)
 
 B_ne = dim_avg_n_Wrap(B_col, (/1,2/))
 B_ne_12 = split_12(B_ne)
 B_clm = avg(B_ne_12(0:4,:))


 B_clm = (/B_clm/10^6/)


 B_col = (/rmMonAnnCycTLL(B_col)/)

 B_nep = dim_avg_n_Wrap(B_col, (/1,2/))

 B_nep = (/B_nep/10^6/)



 ;;;;

 C_col =  dim_sum_wgt_n_Wrap(C, dz, 1, 1)

 C_ne = dim_avg_n_Wrap(C_col, (/1,2/))
 C_ne_12 = split_12(C_ne)
 C_clm = avg(C_ne_12(0:4,:))

 C_clm = (/C_clm/10^6/)

 C_col = (/rmMonAnnCycTLL(C_col)/)

 C_nep = dim_avg_n_Wrap(C_col, (/1,2/))

 C_nep = (/C_nep/10^6/)

 ;;;;;

 
 MSE_nep_12 = split_12(MSE_nep)
 A_nep_12 = split_12(A_nep)
 B_nep_12 = split_12(B_nep)
 C_nep_12 = split_12(C_nep)


 MSE_nep_annual = month_to_season(MSE_nep,"MAM" )
 A_nep_annual = month_to_season(A_nep, "MAM")
 B_nep_annual = month_to_season(B_nep, "MAM")
 C_nep_annual = month_to_season(C_nep, "MAM")
 
 MSE_nep_annual = (/dim_avg_n_Wrap(MSE_nep_12(0:4,:), 0)/)
 A_nep_annual = (/dim_avg_n_Wrap(A_nep_12(0:4,:), 0)/)
 B_nep_annual = (/dim_avg_n_Wrap(B_nep_12(0:4,:), 0)/)
 C_nep_annual = (/dim_avg_n_Wrap(C_nep_12(0:4,:), 0)/)


 rc_mse_mam = trend(MSE_nep_annual)

 rc_A_mam = trend(A_nep_annual)
 
 rc_B_mam = trend(B_nep_annual)

 rc_C_mam = trend(C_nep_annual)
 

 std_mse = dim_stddev_n_Wrap(MSE_nep_annual, 0)
 std_A_mam = dim_stddev_n_Wrap(A_nep_annual, 0)
 std_B_mam = dim_stddev_n_Wrap(B_nep_annual, 0)
 std_C_mam = dim_stddev_n_Wrap(C_nep_annual, 0)
  
 std_mse = std_mse
 std_A_mam = std_A_mam
 std_B_mam = std_B_mam
 std_C_mam = std_C_mam

 
 percent1 = rc_mse_mam/std_mse
 percent2 = rc_A_mam/std_A_mam
 percent3 = rc_B_mam/std_B_mam
 percent4 = rc_C_mam/std_C_mam
 
 print(percent1)
 print(percent2)
 print(percent3)
 print(percent4)

 ;;;

  number = new((/4,2/),float)
  number(:,0) = (/MSE_clm,A_clm,B_clm,C_clm/) 
  number(:,1) = (/rc_mse_mam,rc_A_mam,rc_B_mam,rc_C_mam/) 
 
print(MSE_clm)


;==============================================================
; Mission 5. bar plot
; ==============================================================
  
     wks = gsn_open_wks("pdf","MSE气候态和变化项的bar图1980-2020")          ; send graphics to PNG file

     sres = True
     sres@vpWidthF = 0.7
     sres@vpHeightF = 0.5
     sres@vpXF = .15


     sres@trXMinF = 0.5
     sres@trXMaxF = 4.5
     sres@trYMinF = 0.0
     sres@trYMaxF = 2.5*100

     sres@gsnDraw = False
     sres@gsnFrame = False
     sres@gsnXYBarChart = True
     sres@gsnXYBarChartBarWidth = 0.15           ; change bar widths
     sres@tmXBMode          = "Explicit"         ; explicit labels
     sres@tmXBValues        = (/1,2,3,4/)
     sres@tmXBLabels = (/"MSE","CpT","Lq","Gz"/)
     ; sres@tmYLPrecision = 1
      ; sres@tmYLMode   = "Explicit"        ; explicitly set Y-axis labels
      ; sres@tmYLValues = ispan(-8, 8, 2)*0.01
      ; sres@tmYLLabels = ispan(-8, 8, 2)*0.01
     
     sres@tmXBLabelFontHeightF = 0.0205
     sres@tmXTLabelFontHeightF = 0.0205
     sres@tmYLLabelFontHeightF = 0.0225
     sres@tiMainFontHeightF = 0.025
     sres@gsnRightString = ""
     sres@tiYAxisString = "(10~S~6~N~J/kg)"
     sres@tiYAxisFontHeightF  =0.02
     sres@gsnYRefLine = 0.      ; set the reference line equal to 0
    

     str = ispan(1,2,1)
     color = (/"black","red","blue","green"/)
     ; color = (/"black","gray"/)
     plot = new(2,graphic)
     
     sres@gsnXYBarChartColors = color
     sres@gsnLeftStringFontHeightF = 0.02
     sres@gsnLeftString = "(a) JFMAM Climatology"
     ; plot(i) = gsn_csm_xy(wks,(/0.75+0.15*i,1.75+0.15*i,2.75+0.15*i,3.75+0.15*i/),number(:,i),sres)   ; draw each time series
     plot(0) = gsn_csm_xy(wks,(/0.75+0.25,1.75+0.25,2.75+0.25,3.75+0.25/),number(:,0),sres)   ; draw each time series


     sres@trYMaxF = 1.4
     sres@gsnLeftString = "(a) JFMAM Change"
     plot(1) = gsn_csm_xy(wks,(/0.75+0.25,1.75+0.25,2.75+0.25,3.75+0.25/),number(:,1),sres)   ; draw each time series

     

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
     ; gsn_labelbar_ndc(wks,2,labels(::-1),0.17,0.33,lbres) ; draw right labelbar column
     gsn_labelbar_ndc(wks,2,labels(::-1),0.17,0.23,lbres) ; draw right labelbar column
  
     lbres@lbFillColors       := (/"green","blue"/)
     labels := (/"Lvq","Gz"/)
     ; gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.33,lbres) ; draw left labelbar column
     gsn_labelbar_ndc(wks,2,labels(::-1),0.52,0.23,lbres) ; draw left labelbar column
     
     ; pres = True
     ; gsn_panel(wks, plot, (/1,2/), pres)
   

     Markres = True
     Markres@gsMarkerColor = "gray"
     Markres@gsMarkerSizeF = 0.01
     Markres@gsMarkerIndex = 16
     plot@maker = gsn_add_polymarker(wks, plot(1),ispan(1,4,1), (/std_mse,std_A_mam,std_B_mam,std_C_mam/), Markres)
       draw(plot(1))
     frame(wks)

     ; frame(wks)      
    
end


