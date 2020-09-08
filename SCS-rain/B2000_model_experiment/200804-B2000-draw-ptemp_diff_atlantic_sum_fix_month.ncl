;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
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

  season = (/"spring","summer","autumn","winter"/)
  ; fils1 = "/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.TEMP.mon.0250-0299.nc"  
  fils1 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP.h0.0251-0280.nc"
  ; fils2 = "/home/yangsong3/data-model/lsh/CESM/B/B2000_alt_north_winter/B2000_alt_north_winter.TEMP.h0.0251-0280.nc"

   
  fils2 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_"+season+"_fixmonth/B2000_alt_north_"+season+"_fixmonth.TEMP.h0.0251-0280.nc"

  diff_total = new((/4,4,384,320/),float)

   season = (/"MAM","JJA","SON","DJF"/)


do j = 0,3

    f1=addfile(fils1, "r")
    f2=addfile(fils2(j), "r")
    
    time1 = f1->time
    time2 = f2->time
    time1 = (/time1 -1/)
    time2 = (/time2 -1/)
    date1=cd_calendar(time1, 1)
    date2=cd_calendar(time2, 1)


    timeind1=ind(date1.ge.25601.and.date1.le.28012)
    timeind2=ind(date2.ge.25601.and.date2.le.28012)


    PSL_CT  =  f1->TEMP(timeind1,0,:,:)

  

    PSL_indian  =  f2->TEMP(timeind2,0,:,:)
    TLAT = f2->TLAT
    TLONG = f2->TLONG
   
    PSL_indian!1 = "nlat"
    PSL_indian!2 = "nlon"
    PSL_indian@lat2d   = TLAT
    PSL_indian@lon2d   = TLONG


    printVarSummary(PSL_indian)

    

  do i=0,3

    PSL_CT_clm = dim_avg_n_Wrap(month_to_season(PSL_CT, season(i)), 0)
    PSL_indian_clm = dim_avg_n_Wrap(month_to_season(PSL_indian, season(i)), 0)
   

    PSL_diff_indian = PSL_indian_clm
    PSL_diff_indian = (/PSL_indian_clm - PSL_CT_clm/)

    diff_total(j,i,:,:) = (/PSL_diff_indian/)

   end do 

end do 
;;;;;;;;;;;;;;;



;;;;计算累加

  diff_total_sum = dim_sum_n_Wrap(diff_total,0)
  diff_total_sum!0 ="time" 
  diff_total_sum!1 = "nlat"
  diff_total_sum!2 = "nlon"

  diff_total_sum@lat2d   = TLAT
  diff_total_sum@lon2d   = TLONG

   printVarSummary(diff_total_sum)
;;;;;;;;;;

    wks=gsn_open_wks("pdf", "diff-Potensial-temperature-restore-altantic-fourseason-sum-0280-fixmonth.pdf")

    plot = new(4,graphic)
    str3 = ispan(1, 4, 1) +""
    plot1 =new(4,graphic)
    iplot = new(4,graphic)
  

  ; do i = 1,3
  ; diff_total_sum(i,:,:) = (/diff_total_sum(i,:,:) - diff_total_sum(0,:,:)/)
  ; end do
  ; diff_total_sum(0,:,:) = 0.0

do i = 0,3

;;;;;画图
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    ; res@cnFillMode  = 
    ; res@cnRasterSmoothingOn=True
    res@trGridType = "TriangularMesh"

    ; res@cnFillPalette = "hotcolr_19lev"


    gsn_define_colormap(wks, "BlueWhiteOrangeRed")

    ; gsn_define_colormap(wks, "BlueRed")
    ; gsn_define_colormap(wks, "BlueDarkRed18")
 
 


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.1
    
    ; res@cnMinLevelValF       =   -1.5
    ; res@cnMaxLevelValF       =   1.5
    ; res@cnLevelSpacingF      =  0.1
    
    ; res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,0,0,16,18,20,21/)

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF = 180.0

    ; res@mpMinLatF             = -10.0
    ; res@mpMaxLatF             = 50.0
    ; res@mpMinLonF             = 60.0
    ; res@mpMaxLonF             = 160.0




    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3    
    res@gsnRightString = ""
    res@gsnLeftString = season(i)
 
    plot(i) = gsn_csm_contour_map(wks, diff_total_sum(i,:,:), res)
     

    ;;;画box
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   
    lat1 = 0
    lat2 = 15
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    
    plot1@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres) 
end do
    pres = True
    pres@gsnPanelLabelBar = True
    pres@lbOrientation = "Vertical"
    ; pres@gsnMaximize = True
    
    gsn_panel(wks, plot, (/2,2/), pres)
    print("done script")
  end
  


