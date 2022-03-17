;;********************************************************
;;draw different season long-term change diabatic heat and Q1 HEAT rate profile averaged over SCS and 
;;philippines sea
;;diabatic and convection data from CFSR
;;q1 derived from NCEP_DOE calculate by eqution
;;                               lnshuheng 2020.1.1 
;**********************************************************

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1)/),float)
    do i=0,11   
    data_12(i,:,:) = data(i::12,:)
    end do 
return data_12
end 

begin
                ; send graphics to PNG file
  wks = gsn_open_wks("pdf","Q1profile-indian-strweakmon")                  ; send graphics to PNG file

  
; ==============================================================
; Mission 1. read data
; ==============================================================
;;;read Q1 monthly data from linshuheng
  fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/q1/q1"+"*") ;;;daily数据计算的q1
  ; fils = systemfunc("ls /home/yangsong3/data-observation/NCEP_DOE/monthly/q1/q1"+"*") ;;;daily数据计算的q1
  ; print(fils)
  f3 = addfiles(fils, "r")
  date = cd_calendar(f3[:]->time, 1)

  timeind = ind(date.le.201612.and.date.ge.197901)

  q1 = f3[:]->q1(timeind,:,:,:)

  q1 = (/q1*86400.0/)

  q1_clm = clmMonTLLL(q1)
  q1 = (/calcMonAnomTLLL(q1, q1_clm)/)

  ;;;区域平均
      g_latS=10.
      g_latN=30.
      g_lonW=70.
      g_lonE=90.

      ;  g_latS=0.
      ; g_latN=30.
      ; g_lonW=110.
      ; g_lonE=150.

      q1_indian = wgt_areaave_Wrap(q1(:,:,{g_latS:g_latN},{g_lonW:g_lonE}),1.0,1.0,0)
   ;;;;

  ; printVarSummary(q1)

  ; q1_annual = month_to_annual(q1, 1)
  ; q1_annual = month_to_season(q1_indian, "DJF")
  q1_12 = split_12(q1_indian)
  q1_annual = dim_avg_n_Wrap(q1_12(5:8,:,:), 0)
   printVarSummary(q1_annual)

   ;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    startyear=1979

    pre_jjas = jjas({startyear:2016}) 

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2016-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;

;;;strong and weak monsoon year 
    total_year = ispan(startyear,2016,1)
    std = dim_stddev_Wrap(DMI_ctrl_JJA)

    inds = ind(DMI_ctrl_JJA.gt.1.0*std)
    indw = ind(DMI_ctrl_JJA.lt.-1.0*std)
    indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
;;;

q1_annual_strmon = dim_avg_n_Wrap(q1_annual(inds,:), 0)
q1_annual_wekmon = dim_avg_n_Wrap(q1_annual(indw,:), 0)

; q1_annual_strmon = q1_annual_strmon/dim_max_n_Wrap(q1_annual_strmon, 0)
q1_annual_wekmon = q1_annual_wekmon/dim_max_n_Wrap(q1_annual_wekmon, 0)


;;;intepolate 
  linlog = 2   ; ln(p) interpolation

  pi = q1&lev
           
  xi = q1_annual_strmon

  ; po =(/3.545,7.389,13.97,23.94,37.23,53.11\    
  ; ,70.06,85.44,100.5,118.3,139.1,163.7,192.5,226.5,266.5\
  ; ,313.5,368.8,433.9,510.5,600.5,696.8,787.7,867.2,929.6\    
  ; ,970.6,992.6 /)

  po = (/3.64, 7.59, 14.35, 24.61, 38.26, 54.59, 72.0124,87.82, 103.31, 121.54, 142.99, \
  168.22, 197.90, 232.82, 273.91, 322.24, 379.10, 445.9925, 524.68, \
  609.77, 691.38, 763.40, 820.85,859.53, 887.02, 912.64, 936.19, \
  957.485, 976.325, 992.556/)

   heat_prfl=(/0.00000000 , 0.00000000 , 0.00000000 \
    , 0.00000000 , 0.00000000 , 0.00000000 , 0.00000000 , 0.00000000 ,\
      0.00000004 , 0.00000029 , 0.00000113 , 0.00000244 , 0.00000369 , 0.00000454 \
    , 0.00000674 , 0.00000929 , 0.00000970 , 0.00000970 , 0.00001157 ,\
      0.00000656 , 0.00000461 , 0.00000470 , 0.00000443 , 0.00000517 , 0.00000480\
    , 0.00000599 , 0.00000604 , 0.00000012 , 0.00000000 , 0.00000000/)
  heat_prfl = heat_prfl*86400.0

  print(heat_prfl)
  exit()
  xo = int2p_Wrap(pi,xi,po,linlog)

  xo = xo/dim_max_n_Wrap(xo, 0)

printVarSummary(xo)
print(xo)
; xo(0:7) = 0.0
; xo(23:25)=0.0
xo(0:7) = 0.0
xo(26:29)=0.0
asciiwrite("heat_mar_double_30.txt",xo)
exit()
; xo   =(/0., 0., 0., 0., 0., 0., 0., 0., -0.04108525, 0.1043775,\
;                              0.243452, 0.3307469, 0.3690013, 0.5299402, 0.748916,\
;                              0.9155024, 0.9713629, 0.9999998, 0.9381863, 0.6111271,\
;                              0.3761931, 0.3309714, 0.2640985, 0.1148081, 0.02364605, -0.01868893/)



; q1_annual_strmon = q1_annual_strmon/dim_max_n_Wrap(q1_annual_strmon, 0)


;;;;
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/2.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       :=  (/"blue","blue","red","red"/)          ; change line color
    resxy@xyMarkLineModes    :=  (/"MarkLines","MarkLines","MarkLines","MarkLines"/)
    resxy@xyDashPatterns := (/0,1,0,1/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           := (/16,16,7,7/)
    resxy@xyMarkerColors :=  (/"blue","blue","red","red"/)

    resxy@gsnLeftStringFontHeightF = 0.02
    ; resxy@gsnLeftString = /order(i)+" " +season(i) + " DPLH DIABH Q1 Profile"

  

    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
    resxy@tmXBMajorThicknessF=1.0
    resxy@tmYLMajorThicknessF=1.0

    
    resxy@tmXBMode  =  "Explicit"
    resxy@trXMaxF =  1.
    resxy@trXMinF = -1.
    resxy@tmXBValues = ispan(-10,10,2)*0.1
    resxy@tmXBLabels = ispan(-10,10,2)*0.1
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = ispan(-10,10,1)*0.1
    resxy@tiYAxisString = ""
   
    ; resxy@pmLegendDisplayMode = "Always"
    ; resxy@lgBoxLinesOn = False
    ; resxy@xyExplicitLegendLabels  =  (/"DPLH Clm","DPLH Clm+Change","DIABH Clm","DIABH Clm+Change"/)     ;加线的标注
    resxy@trYReverse        = True                       ; reverse Y-axis

    plot1 = gsn_csm_xy(wks, (/xo/),po ,resxy)
    ; plot1 = gsn_csm_xy(wks, (/q1_annual_strmon/),pi ,resxy)
    
draw(plot1)
frame(wks)

print("done stript")

   


end 