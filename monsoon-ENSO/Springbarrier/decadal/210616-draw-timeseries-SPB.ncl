;;*************************************************
;;draw SPB 10 years windows time series 
;;
;;*************************************************************
function cal_mean(length,nave,data)
begin
  EAESM_series = new((/length-nave+1/),float)
  do t = nave,length  
    EAESM  = dim_avg_n_Wrap(data(t-nave:t-1),0)
    EAESM_series(t-nave) =  EAESM
    delete(EAESM)
  end do 
return EAESM_series
end 


begin
;;;read SPB data
    diri="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/decadal/"
    fils1=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend.nc")
    ; fils1=systemfunc("ls "+diri+"SPB_series_1979-2019.nc")
    f1=addfile(fils1, "r")

    SPB_series = f1 ->SPB_series
;;;


;;;读入ERA wind data

   fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
   f1 = addfile(fils1, "r")
   time1 = f1->time
   date1=cd_calendar(time1, 1)
   timeind1=ind(date1.ge.197901.and.date1.le.201912)
 
  
;;;;U,V
    u_CT_200  = short2flt(f1->u(timeind1,{200},:,:))
    u_CT_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_CT_850 = (/rmMonAnnCycTLL(u_CT_850)/)
    u_CT_200 = (/rmMonAnnCycTLL(u_CT_200)/)
;;;

;;;; read  HadISST data
    startmon = 197901
    endmon = 201912

    diri = "/home/yangsong3/data-observation/linshh/data/sst/"
    fils1 = systemfunc("ls "+diri+"HadISST_sst.nc")
    f1 = addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta_annual = month_to_annual(ssta, 1)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
;;;

;;;;;EAESM指数 : 850hPA ZONAL wind diff 110-140 22.5-32.5 - 90-130,5-15
    DMI_ctrl = u_CT_850(:,1,1)

    DMI_ctrl  =  (/dim_avg_n_Wrap(u_CT_850(:,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_CT_850(:,{5:15},{90:130}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  ;;;jja 少一年，从25601开始
;;;


;;;;;Australian 指数 
    ASM = u_CT_850(:,1,1)

    ASM =  (/dim_avg_n_Wrap(u_CT_850(:,{-5:-15},{110:130}), (/1,2/))/)

    ASM_DJF = month_to_season(ASM, "DJF")  ;;;jja 少一年，从25601开始
;;;

;;;ISM index  850hPA ZONAL wind diff 5-15N 40-80E - 20-30N,70-90E
    ISM = u_CT_850(:,1,1)

    ISM = (/dim_avg_n_Wrap(u_CT_850(:,{5:15},{40:80}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_CT_850(:,{20:30},{70:90}), (/1,2/))/)

    ISM_JJA = month_to_season(ISM, "JJA")  ;;;jja 少一年，从25601开始
;;;;


;;;; AMO index 和 nino3.4 index 
    nrows = 2015-1948+1
    ncols = 13 
    AMO  = asciiread("/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/decadal/AMO_10smooth.txt", (/nrows,ncols/), "float")
    
    AMO_index = ndtooned(AMO(:,1:)) ;; TOTAL :1948-2015 

    AMO_annual = month_to_annual(AMO_index, 1)
;;;;

;;; IPO index 
    ;;;area mean 
      rad    = 4.0*atan(1.0)/180.0
      re     = 6371220.0
      rr     = re*rad

      dlon   = abs(lon(2)-lon(1))*rr
      dx     = dlon*cos(lat*rad)
      ;          
      jlat  = dimsizes( lat )
      dy     = new ( jlat, typeof(dx))      ; close enough
      
      dy(0)  = abs(lat(2)-lat(1))*rr
      dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5   
      dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr

      area   = dx*dy                                ; cell area function of latitude only
      copy_VarCoords(lat, area)
      printVarSummary(area)

     IPO_index = wgt_areaave_Wrap(ssta_annual(:,{-20:20},{40:100}), area({-20:20}), 1.0, 0) 

;;;;






;;;get 10years running windows EAESM strength 
    length = 2019-1979+1
    nave = 10 

    EAESM_series = cal_mean(length,nave,DMI_ctrl_JJA)
    ASM_series = cal_mean(length,nave,ASM_DJF)
    ISM_series = cal_mean(length,nave,ISM_JJA)

    ASM_series2 = runave_n(ASM_DJF, 10, 0, 0)

;;;;


;;;;
   EAESM_series = (/dim_standardize_Wrap(EAESM_series, 1)/)

   ; ASM_series = (/dim_standardize_Wrap(ASM_series, 1)/)

   ; ASM_series2 = (/dim_standardize_Wrap(ASM_series2, 1)/)

   ISM_series = (/dim_standardize_Wrap(ISM_series, 1)/)

   SPB_series = (/dim_standardize_Wrap(SPB_series, 1)/)

   ASM_DJF = (/dim_standardize_Wrap(ASM_DJF, 1)/)

   AMO_annual = (/dim_standardize_Wrap(AMO_annual, 1)/)
;;;


; ;;;corr 
;   corr = escorc(EAESM_series, SPB_series)
  ; corr2 = escorc(ISM_series, SPB_series)
  ; print(corr2)
  ; exit()
; ;;;



;;;draw plot 

    wks = gsn_open_wks("pdf", "SPB_intensity_timeseries")

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "SPB intensity 10-year windows"   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02
    res@tiYAxisString = "SPB intensity"           ; y axis title
    res@tiXAxisString = "Year"

    ; res@xyLineColors      = (/"black","red","blue","green"/)  ; line colors
    res@xyLineColor      = (/"black"/)  ; line colors
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    ; res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
    res@gsnYRefLine = 0.0

    res@tmXBMode="Explicit"

    res@trXMaxF = 2010
    res@trXMinF = 1949

    res@tmXTOn = False
    res@tmYROn  = False

    res@tmXBValues=ispan(1950,2010,5)
    res@tmXBLabels=ispan(1950,2010,5)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(1949,2010,1)
    res@tmYRMode="Automatic"

    res@vpHeightF = 0.3
    res@vpWidthF = 0.6

    plot = gsn_csm_xy(wks, ispan(1904, 2015, 1), SPB_series, res)

    resp = True
    resp@gsLineColor = "red"
    resp@gsLineDashPattern = 0
    resp@gsLineThicknessF = 4.0

    ; plot_line = gsn_add_polyline(wks, plot, ispan(1979,2010,1), ASM_series, resp)
    ; ; plot_line = gsn_add_polyline(wks, plot, ispan(1979,2010,1), ISM_series, resp)
    ; resp@gsLineColor = "blue"
    ; plot_line1 = gsn_add_polyline(wks, plot, ispan(1979,2010,1), ASM_series2, resp)
    ; ; plot_line1 = gsn_add_polyline(wks, plot, ispan(1949,2010,1), AMO_annual(2:dimsizes(AMO_annual)-6), resp)
  


    ; ;;添加legend 
    ; gres = True
    ; gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;       ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ; gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ; lineres = True
    ; ; lineres@lgLineColors = (/"black"，/) ; line colors
    ; lineres@lgLineColors       = res@xyLineColors ; colors chosen
    ; lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    ; lineres@LineLengthPercent = 9.
    ; lineres@lgDashIndexes   = res@xyDashPatterns                     ; expressed as %, 0->100, length of line

    ; textres = True

    ; textres@lgLabels = (/"MSE","CpT","Lq","Hgt"/)
    ; textres@lgLabelFontHeightF = 0.01
    ; ; textres@lgLabelOpacity = (/0/) 

    ; plot(0) = simple_legend(wks,plot(0),gres,lineres,textres)   


    ; plot(1) = gsn_csm_xy(wks, ispan(1950, 2019, 1), (/MSE_nep_MAM,MSE_nep_JJA,MSE_nep_SON,MSE_nep_DJF/), res)

draw(plot)
frame(wks)
end 