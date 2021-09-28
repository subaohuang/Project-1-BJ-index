;;********************************************************
;;计算SST倾向以及SST值与GPCP降水的趋势
;;其中热通量的方向都是向上为正，向下为负
;;的相关系数场，
;;NCEP-NCAR的精度为1.87×1.9
;;
;;                                         lnshuheng 2020.02.18
;**********************************************************
function chazhi(data)
 begin
   xi = data&lon
   yi = data&lat
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
   data_chazhi&lat = yo
   data_chazhi&lon = xo
   copy_VarAtts(data, data_chazhi)
 
   return data_chazhi
end 
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end



;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =199101
    endmon = 201912

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
    order  = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/)
 



 ;;;;read data

  lat1 = 0
  lat2 = 20
  lon1 = 240
  lon2 = 280

    var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc","omega.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")
     f4 = addfile(diri+var(3), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.endmon.and.date.ge.startmon)
     
    t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    t=t+273.15

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q = q/1000.0  ;;;convert to kg/kg

    hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    
    time = f4->time 
    date := cd_calendar(time, 1)
    timeind :=ind(date.le.endmon.and.date.ge.startmon)
    omega = f4->omega(timeind,{700},{lat1:lat2},{lon1:lon2})
    omega = (/omega*-100/)

    omega_a =(/rmMonAnnCycTLL(omega)/)
    copy_VarCoords(omega, omega_a)
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
;;;;


    ; wks   = gsn_open_wks("pdf", "CMAP&olr回归到印度洋SST的回归系数")
    wks   = gsn_open_wks("pdf", "NCEP_NCAR q-t-gz change profile 回归局地omega并乘上趋势")
  
    plot=new(4,graphic)
    plot1 = plot

    str = ispan(1,4,1)+ ""
    str2 = ispan(5,8,1)+ ""
 
    iplot = new(4,graphic)

    season =  (/"MAM","JJA","SON","DJF"/)
;;;;;



i = 0


  ssta_annual = month_to_season(omega_a, season(i))
  printVarSummary(ssta_annual)


  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), dim_avg_n_Wrap(ssta_annual, (/1,2/)), 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum


  print("*******ssta_index_rc="+ssta_index_rc)


  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
 
  ;;;;;要分析的变量


    t_MAM = month_to_season(t, season(i))
    q_MAM = month_to_season(q, season(i))
    hgt_MAM = month_to_season(hgt, season(i))



    t_nep_MAM = dim_avg_n_Wrap(t_MAM, (/2,3/))

    q_nep_MAM = dim_avg_n_Wrap(q_MAM, (/2,3/))

    hgt_nep_MAM = dim_avg_n_Wrap(hgt_MAM, (/2,3/))

  
    ssta_index = dim_avg_n_Wrap(ssta_annual, (/1,2/)) 


  ; ;;;计算各变量的实际趋势
  ;   p_rc := regCoef_n(ispan(1,yearnum,1),t_nep_MAM,0,0)
  ;   u_rc := regCoef_n(ispan(1,yearnum,1),q_nep_MAM,0,0)
  ;   v_rc := regCoef_n(ispan(1,yearnum,1),hgt_nep_MAM,0,0)
  ;   copy_VarCoords(t_nep_MAM(0,:), p_rc)
  ;   copy_VarCoords(q_nep_MAM(0,:), u_rc)
  ;   copy_VarCoords(hgt_nep_MAM(0,:), v_rc)
  ;   p_rc = (/p_rc*yearnum/)
  ;   u_rc = (/u_rc*yearnum/)
  ;   v_rc = (/v_rc*yearnum/)

  ;;;;去趋势
    t_nep_MAM =  (/dtrend_msg_n(ispan(1,yearnum,1),t_nep_MAM,False,True,0)/)
    q_nep_MAM =  (/dtrend_msg_n(ispan(1,yearnum,1),q_nep_MAM,False,True,0)/)
    hgt_nep_MAM =  (/dtrend_msg_n(ispan(1,yearnum,1),hgt_nep_MAM,False,True,0)/)
  


  ;;;correlation with sst
    t_nep_mam_change := regCoef_n(ssta_index, t_nep_MAM, 0, 0)
    copy_VarCoords(t_nep_MAM(0,:), t_nep_mam_change)

    ; t_nep_mam_change = (/t_nep_mam_change*ssta_index_rc/)
   ;;;;
   
    ;;;normalized before reression 
    ssta_index = (/dim_standardize_n(ssta_index, 1, 0)/)
    ;;;;

    ;;;correlation with sst
    q_nep_mam_change := regCoef_n(ssta_index, q_nep_MAM, 0, 0)
    copy_VarCoords(q_nep_MAM(0,:), q_nep_mam_change)

    ; q_nep_mam_change = (/q_nep_mam_change*ssta_index_rc/)
    ;;;;;


    ;;;correlation with sst
   hgt_nep_mam_change := regCoef_n(ssta_index, hgt_nep_MAM, 0, 0)
    copy_VarCoords(hgt_nep_MAM(0,:),hgt_nep_mam_change)

   ; hgt_nep_mam_change = (/hgt_nep_mam_change*ssta_index_rc/)
    ;;;;;


   printMinMax(t_nep_mam_change, 1)
   
   mse_change = hgt_nep_mam_change + t_nep_mam_change + q_nep_mam_change 




  print("*********开始画图********")

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
; end do 

; end do

  pres = True
  pres@gsnMaximize = True
  
  pres@gsnPanelBottom = 0.05

  ; gsn_panel(wks, plot, (/2,2/), pres)
  draw(plot(0))
  frame(wks)
  
; system("convert -geometry 1800x1800 -density 600x600 -trim  GPCP&olr回归到印度洋SST的回归系数.png  GPCP&olr回归到印度洋SST的回归系数_high.png")     

 
 print("done")

end