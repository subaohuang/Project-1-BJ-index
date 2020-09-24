;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
function t_test(series1,series2,series_mam)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_pt = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,True)
      
      
     
      t = student_t(prob_pt(0,:,:), simple_size)
      std = dim_stddev_n(series1,0)
      interv = t*std/sqrt(simple_size)  ;;;置信区间     
    
      interv = where(interv.eq.interv@_FillValue, 0.0, interv)
       
       printVarSummary(interv)
    
      diff = avg1 
      diff = (/abs(avg2 - avg1 -series_mam)/)
      prob_p = diff
      printVarSummary(prob_p)
      printVarSummary(interv)
      prob_p = where(prob_p.ge.interv, 0.05, 2.0)
      prob_p = where(prob_p.eq.prob_p@_FillValue, 2.0, prob_p)
    
   

copy_VarCoords(avg1, prob_p)


return prob_p
end 



begin
  ; fils1 = "/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.TEMP.mon.0250-0299.nc"  
fils1 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h0.0251-0280.nc"
fils2 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year/B2000_alt_north_year.TEMP_taux_tauy.h0.0251-0280.nc"
 
 
  
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    
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

    
    season = (/"MAM","JJA","SON","DJF"/)
    wks=gsn_open_wks("pdf", "diff-Potensial-temperature-restore-altantic-year-fourseason-0280-miuns-spring.pdf")

    plot = new(4,graphic)
    str3 = ispan(1, 4, 1) +""
    plot1 =new(4,graphic)
    iplot = new(4,graphic)




  do i=0,3
    
    PSL_CT_clm = dim_avg_n_Wrap(month_to_season(PSL_CT, season(i)), 0)
    PSL_indian_clm = dim_avg_n_Wrap(month_to_season(PSL_indian, season(i)), 0)
   

    PSL_diff_indian = PSL_indian_clm
    PSL_diff_indian = (/PSL_indian_clm - PSL_CT_clm/)

    if season(i).eq."MAM"
        PSL_diff_spring = PSL_diff_indian
    end if

      PSL_diff_indian = (/PSL_diff_indian - PSL_diff_spring/)

    ;;;t_test 
    prob_slp_indian := t_test(PSL_CT,PSL_indian, PSL_diff_spring)
    prob_slp_indian!0 = "nlat"
    prob_slp_indian!1 = "nlon"
    prob_slp_indian@lat2d   = TLAT
    prob_slp_indian@lon2d   = TLONG

     

    
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

    res@cnFillPalette = "BlueDarkRed18"


    res@tmXBLabelFontHeightF=0.01
    res@tmYLLabelFontHeightF=0.01
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0
    res@gsnLeftStringFontHeightF = 0.02


    res@gsnRightString = ""
    res@gsnLeftString = season(i) + " Minus MAM"



    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.05
    
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

    
    

    plot(i) = gsn_csm_contour_map(wks, PSL_diff_indian, res)
     

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
   

     ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True
     res_pro@gsnLeftString = ""
    res_pro@gsnRightString = ""

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
   
    res_pro@gsnCenterString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  

 
    iplot(i)=gsn_csm_contour(wks,prob_slp_indian,res_pro)

      if i.ne.0 then 
    overlay(plot(i), iplot(i))
    end if

    ; ;;;打点
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@gsnLeftString = ""
    ; res_pro@gsnRightString = ""
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False
    ; res_pro@trGridType = "TriangularMesh"

    ; iplot=new(4,graphic)

 
    ; iplot(i)=gsn_csm_contour(wks,prob_slp_indian,res_pro)

    

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; opt@gsnShadeMid     = 17         ; stipple
    ; ; opt@gsnShadeFillDotSizeF =  0.0025
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

   
    ; if i.ne.0 then 
    ; iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    ; overlay(plot(i), iplot(i))
    ; end if
end do
    pres = True
     pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    ; pres@lbOrientation  = "Vertical"
    pres@lbOrientation  = "Horizontal"
    pres@pmLabelBarOrthogonalPosF = -0.02
    ; pres@pmLabelBarWidthF  = 0.03
    pres@pmLabelBarHeightF = 0.03
    pres@lbLabelStride = 2  
    
    gsn_panel(wks, plot, (/2,2/), pres)
    print("done script")
  end
  
  

  

