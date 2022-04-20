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

function slid_cor(nave,x,y,dimnum);滑动相关,一维的情况
begin  
if dimnum.eq.1 then
length = dimsizes(x)
R=new(length-nave+1,typeof(x))
do t=nave,length
  R(t-nave)=escorc(x(t-nave:t-1),y(t-nave:t-1))   
end do
end if
return R
end 

begin
;;;read SPB data
    diri="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/decadal/"
    ; fils1=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend.nc")
    fils1=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend-fang2-method.nc")
     ; fils1=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend-3456-method.nc")
    f1=addfile(fils1, "r")

    SPB_series = f1 ->SPB_series
;;;

; ;;;read 20 years SPB data
;     diri="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/"
;     ; fils2 = systemfunc("ls "+diri+"SPB_series_1949-2019-20windows.nc")
;     fils2=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend-20windows.nc")
;     f2 = addfile(fils2, "r")

;     SPB_series_20years = f2 ->SPB_series
; ;;;


; ;;;read 30 years SPB data
;     diri="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/"
;     ; fils3 = systemfunc("ls "+diri+"SPB_series_1949-2019-30windows.nc")
;     fils3 = systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend-30windows.nc")
;     ; fils1=systemfunc("ls "+diri+"SPB_series_1979-2019.nc")
;     f3 = addfile(fils3, "r")

;     SPB_series_30years = f3 ->SPB_series
; ;;;


;;;; read  HadISST data
    startmon = 189901
    endmon = 201912

    yearstart = 1899
    yearend = 2019

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
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
    nino34_annual = month_to_annual(nino34, 1)

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:dimsizes(nino34_djf)-1)
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1)

    nino34_jja = month_to_season(nino34, "JJA")
;;;


;;;get 10years running windows EAESM strength 
    length = 2019-1899+1
    nave = 10 

    ; EAESM_series = cal_mean(length,nave,DMI_ctrl_JJA)
    ; ASM_series = cal_mean(length,nave,ASM_DJF)
    ; ISM_series = cal_mean(length,nave,ISM_JJA)
    nino34_series = cal_mean(length,nave,nino34_annual)
;;;;

;;;; AMO index 和 nino3.4 index 
    nrows = 2015-1861+1
    ncols = 13 
    AMO  = asciiread("/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/decadal/amon.sm.long.data.txt", (/nrows,ncols/), "float")
    
    AMO_index = ndtooned(AMO(:,1:)) ;; TOTAL 1861-2015 

    AMO_annual = month_to_annual(AMO_index, 1)

    ind1 = 1899-1861
    ind2 = 2015-1861
;;;;


;;;;;; IPO index 
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

      IPO_index = wgt_areaave_Wrap(ssta(:,{-10:10},{170:270}), area({-10:10}), 1.0, 0) - \
               (/wgt_areaave_Wrap(ssta(:,{25:45},{140:215}), area({25:45}), 1.0, 0) + wgt_areaave_Wrap(ssta(:,{-50:-15},{150:200}), area({-50:-15}), 1.0, 0)/)*0.5

      IPO_index = runave_n(IPO_index, 121, 0, 0)

      IPO_annual = month_to_annual(IPO_index, 1)
;;;;
    ; SPB_series = (/dim_standardize_Wrap(SPB_series, 1)/)

    ; SPB_series_20years = (/dim_standardize_Wrap(SPB_series_20years, 1)/)

    ; SPB_series_30years = (/dim_standardize_Wrap(SPB_series_30years, 1)/)

    ; ASM_DJF = (/dim_standardize_Wrap(ASM_DJF, 1)/)
    AMO_annual = (/dim_standardize_Wrap(AMO_annual, 1)/)
    IPO_annual = (/dim_standardize_Wrap(IPO_annual, 1)/)
    nino34_series = (/dim_standardize_Wrap(nino34_series, 1)/)


;;;读入ERA wind data

    ; fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = 201212                 ;;; 1900-2019
    startmon = 190001
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  

 
;;;;U,V
    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_200  = f1->uwnd(timeind1,{200},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)

    DMI_ctrl = u_850(0:dimsizes(u_850&time)-1,1,1)
    DMI_ctrl  =  (/dim_avg_n_Wrap(u_850(0:dimsizes(u_850&time)-1,{22.5:32.5},{110:140}), (/1,2/)) - \
                  dim_avg_n_Wrap(u_850(0:dimsizes(u_850&time)-1,{5:15},{90:130}), (/1,2/))/)

    DMI_ctrl_JJA = month_to_season(DMI_ctrl, "JJA")  
;;;

;;滑动相关
    ;;nino34_djf:1899-2019  DMI_ctrl_JJA: 1900-2012
    ; corr = slid_cor(10,DMI_ctrl_JJA,nino34_djf(1:dimsizes(nino34_djf)-8),1) ;;; 1905-2007
    ; corr = slid_cor(10,DMI_ctrl_JJA,nino34_djf(0:dimsizes(nino34_djf)-9),1) ;;; 1905-2007
    corr = slid_cor(10,DMI_ctrl_JJA,nino34_jja(1:dimsizes(nino34_djf)-8),1) ;;; 1905-2007
;;;
    ; SPB_series = (/dim_standardize_n_Wrap(SPB_series, 1, 0)/)
    ; corr = (/dim_standardize_n_Wrap(corr, 1, 0)/)
    spb_corr_cor = escorc(SPB_series(16:dimsizes(SPB_series)-8), corr(15:))

    print(spb_corr_cor)
    p = rtest(spb_corr_cor, dimsizes(corr)-1, 0)
    print(p)
;;;;


;;;draw plot 
    wks = gsn_open_wks("pdf", "SPB_timing_timeseries&-jjaENSO-monsoon relation-fangmethod2")

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = "SPB timing 10-year windows"   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02
    res@tiYAxisString = "SPB intensity (month)"           ; y axis title
    res@tiXAxisString = "Year (the time mid window)"

    ; res@xyLineColors      = (/"black","red","blue","green"/)  ; line colors
    res@xyLineColor      = (/"black"/)  ; line colors
    res@xyMonoLineThickness = True
    res@xyLineThicknessF = (/4.0/)        ; line thicknesses
    ; res@xyDashPatterns    = (/0.0,0.0,0.0,0.0/)        ; line patterns
    res@gsnYRefLine = 0.0

    res@tmXBMode="Explicit"

    res@trXMaxF = yearend -5
    res@trXMinF = yearstart

    ; res@trYMaxF = 4.0
    ; res@trYMinF = -4.0

    res@tmXTOn = False
    res@tmYROn  = True


    res@tmXBValues=ispan(yearstart+1,yearend-5,10)
    res@tmXBLabels=ispan(yearstart+1,yearend-5,10)+""
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(yearstart+1,yearend-5,2)
    res@tmYRMode="Automatic"

    res@vpHeightF = 0.3
    res@vpWidthF = 0.6

    ; plot = gsn_csm_xy(wks, ispan(1949, 2010, 1), SPB_series, res)
    ; plot = gsn_csm_xy(wks, ispan(1904, 2015, 1), SPB_series, res)
    ; plot = gsn_csm_xy(wks, ispan(1909, 2010, 1), SPB_series_20years, res)

    ; resp = True
    ; resp@gsLineColor = "red"
    ; resp@gsLineDashPattern = 0
    ; resp@gsLineThicknessF = 4.0

    resp = res
    resp@xyLineColor = (/"red"/)
    resp@tiYAxisString = "Corr DJF(-1)nino34& JJA SEASM "           ; y axis title
    resp@tiXAxisString = ""
    resp@tiYAxisFontColor = "red"
    resp@tiYAxisFontHeightF = 0.015
    resp@tmYRLabelFontColor = "red"
    
    ; resp@tmBorderLineColor = "red"
    resp@tmYRMajorLineColor = "red"
    resp@tmYRMinorLineColor = "red"
    
    ; plot_line = gsn_add_polyline(wks, plot, ispan(1905,2008,1), corr, resp)


    plot = gsn_csm_x2y2(wks,ispan(1904, 2015, 1), ispan(1905,2008,1),SPB_series,corr,res,resp)
    ; plot_line = gsn_add_polyline(wks, plot, ispan(1979,2010,1), ISM_series, resp)
    ; plot_line = gsn_add_polyline(wks, plot, ispan(1899,2010,1), SPB_series_20years, resp)
    ; resp@gsLineColor = "red"
    ; plot_line2 = gsn_add_polyline(wks, plot, ispan(1899,2010,1), SPB_series_30years, resp)
    ; resp@gsLineColor = "blue"
    ; plot_line3 = gsn_add_polyline(wks, plot, ispan(1899,2019,1), IPO_annual, resp)

    ; ; resp@gsLineColor = "red"
    ; ; plot_line4 = gsn_add_polyline(wks, plot, ispan(1899,2015,1), AMO_annual(ind1:ind2), resp)

    ;  resp@gsLineColor = "red"
    ; plot_line4 = gsn_add_polyline(wks, plot, ispan(1904,2015,1), nino34_series, resp)



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