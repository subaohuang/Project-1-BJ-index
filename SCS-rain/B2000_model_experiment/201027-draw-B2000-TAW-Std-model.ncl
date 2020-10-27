;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 
begin

    ; season = "MAM"
    ; fils1 = "/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.TEMP.mon.0250-0299.nc"  
    ; fils1 = "/home/ys17-19/lsh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.TEMP.mon.0251-0280.nc"  
    ; fils2 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year_PACCLM/B2000_alt_north_year_PACCLM.TEMP.h0.0251-0280.nc"  
    ; fils1 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP.h0.0251-0280.nc"  
    


    fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/" + (/"B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h0.0251-0280.nc"/)

    fils2 = "/home/yangsong3/data-observation/linshh/CESM-data/B/" + (/"B2000_alt_north_year/B2000_alt_north_year.TEMP_taux_tauy.h0.0251-0280.nc", \
    "B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer.TEMP_taux_tauy.h0.0251-0280.nc"/)


    leftstring = (/"Atlantic Warming ", "Atlantic Warming and WEP fixed"/)



    wks = gsn_open_wks("pdf", "201027-TAW-std-change-WNPCLM")
    plot = new(1,graphic)

    mid  = new(1,graphic)
    f1=addfile(fils1, "r")
    time1 = f1->time 
    time1 = (/time1 -1/)
    date1=cd_calendar(time1, 1)

    timeind1=ind(date1.ge.25601.and.date1.le.28012)


    sst_CT  :=  f1->TEMP(timeind1,0,:,:)

    TLAT := f1->TLAT
    TLONG := f1->TLONG

    sst_CT!1 = "nlat"
    sst_CT!2 = "nlon"
    sst_CT@lat2d   = TLAT
    sst_CT@lon2d   = TLONG
    

   ;;;;差值TEMP 
    sst_CT_chazhi = PopLatLon(sst_CT,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_CT_chazhi)
   ;;;;
   
    sst_CT_chazhi =(/rmMonAnnCycTLL(sst_CT_chazhi)/)


    sst_ct_mean = dim_avg_n_Wrap(sst_CT_chazhi(:,{-5:5},{120:280}), 1)
    
    std_ct = dim_stddev_n_Wrap(sst_ct_mean, 0)

    title  = (/"Standard deviations of SSTA "/)

   data_all = new((/2,65/),float)
do i=0,1
 
    f2=addfile(fils2(i), "r")
    time2 = f2->time
    time2 = (/time2 -1/)
    date2=cd_calendar(time2, 1)
    timeind2=ind(date2.ge.25601.and.date2.le.28012)


    sst_indian  :=  f2->TEMP(timeind2,0,:,:)
    TLAT := f2->TLAT
    TLONG := f2->TLONG
   
    sst_indian!1 = "nlat"
    sst_indian!2 = "nlon"
    sst_indian@lat2d   = TLAT
    sst_indian@lon2d   = TLONG

    ;;;;差值TEMP 
    sst_indian_chazhi = PopLatLon(sst_indian,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_indian_chazhi)
    

    sst_indian_chazhi =(/rmMonAnnCycTLL(sst_indian_chazhi)/)

    sst_indian_mean = dim_avg_n_Wrap(sst_indian_chazhi(:,{-5:5},{120:280}), 1)

    std_indian = dim_stddev_n_Wrap(sst_indian_mean, 0)

    x = sst_indian_mean&lon


    data_all(i,:)  = std_indian
end do 


    
    
;;;;;draw plot

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiMainString  = title
    res@tiYAxisString = "SST anomaly s.d."           ; y axis title
   
    res@tiMainFontHeightF = 0.015
    res@xyLineColors      = (/"black","red","blue"/)  ; line colors
    res@xyLineThicknesses = (/4.0,4.0,4.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0,0.0/)        ; line patterns

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 280.0
    res@trXMinF = 120.0

    res@trYMinF = 0.0
    res@trYMaxF = 1.8

    

    res@tmXTOn = False
    res@tmYROn  = False
    ; res@tmYLOn  = False

    ; res@tmXBValues=(/0,90,180,270,357.5/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(0,360,10)
    res@vpWidthF = 0.6
    res@vpHeightF = 0.4

    res@tmXBValues=ispan(135,270,45)
    res@tmXBLabels=(/"135~S~o~N~E","180~S~o","135~S~o~N~W","90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(120,280,15)
    res@tmYRMode="Automatic"


    plot = gsn_csm_xy(wks,x,(/std_ct,data_all(0,:),data_all(1,:)/),res)          ; Draw an XY plot with 1 curve.

    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"black","red","blue"/) ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"CTRL","AW","WEPCLM"/)
    textres@lgLabelFontHeightF = 0.01
    mid = simple_legend(wks,plot,gres,lineres,textres)   

   draw(plot)
   frame(wks)
end
