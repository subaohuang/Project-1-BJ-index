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


  ; fils1 = "/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.TEMP.mon.0250-0299.nc"  
  fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/C/lsh_C_CTRL/lsh_C_CTRL.TEMP_TAUX_TAUY.h0.0001-0030.nc"
  fils2 = "/home/yangsong3/data-observation/linshh/CESM-data/C/lsh_C_addaltwind_spring_new/lsh_C_addaltwind_spring_new.TEMP_TAUX_TAUY.h0.0001-0030.nc"
 
  
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    
    time1 = f1->time
    time2 = f2->time
    time1 = (/time1 -1/)
    time2 = (/time2 -1/)
    date1=cd_calendar(time1, 1)
    date2=cd_calendar(time2, 1)
    
    print(date1)

    timeind1=ind(date1.ge.101.and.date1.le.01012)
    timeind2=ind(date2.ge.101.and.date2.le.01012)


    PSL_CT  =  f1->TEMP(timeind1,0,:,:)

  

    PSL_indian  =  f2->TEMP(timeind2,0,:,:)
    TLAT = f2->TLAT
    TLONG = f2->TLONG
   
    PSL_indian!1 = "nlat"
    PSL_indian!2 = "nlon"
    PSL_indian@lat2d   = TLAT
    PSL_indian@lon2d   = TLONG

     ;;;;U,V
    u_CT  =  f1->TAUX(timeind1,:,:)
    v_CT  =  f1->TAUY(timeind1,:,:)
    
    printVarSummary(u_CT)

    u_id  =  f2->TAUX(timeind2,:,:)
    v_id  =  f2->TAUY(timeind2,:,:)
  
    
    season = (/"MAM","JJA","SON","DJF"/)
    wks=gsn_open_wks("pdf", "201016-diff-Potensial-temperature-by-addwind-POPmodel.pdf")

    plot = new(4,graphic)
    str3 = ispan(1, 4, 1) +""
    plot1 =new(4,graphic)
    iplot = new(4,graphic)

  do i=0,3

    PSL_CT_clm = dim_avg_n_Wrap(month_to_season(PSL_CT, season(i)), 0)
    PSL_indian_clm = dim_avg_n_Wrap(month_to_season(PSL_indian, season(i)), 0)
   

    PSL_diff_indian = PSL_indian_clm
    PSL_diff_indian = (/PSL_indian_clm - PSL_CT_clm/)

    

    ;;;t_test 
    prob_slp_indian = t_test(PSL_indian,PSL_CT)
   
    prob_slp_indian!0 = "nlat"
    prob_slp_indian!1 = "nlon"
    prob_slp_indian@lat2d   = TLAT
    prob_slp_indian@lon2d   = TLONG
    ;;;;


   
    u_CT_clm = dim_avg_n_Wrap(u_CT, 0)
    u_indian_clm = dim_avg_n_Wrap(u_id, 0)

    u_diff = (/u_CT_clm/)
    ; u_diff = (/u_indian_clm - u_CT_clm/)

    v_CT_clm = dim_avg_n_Wrap(v_CT, 0)
    v_indian_clm = dim_avg_n_Wrap(v_id, 0)

    v_diff = (/v_CT_clm/)
    ; v_diff = (/v_indian_clm - v_CT_clm/)

    
    u_diff!0 = "nlat"
    u_diff!1 = "nlon"
    u_diff@lat2d   = TLAT
    u_diff@lon2d   = TLONG

    
    v_diff!0 = "nlat"
    v_diff!1 = "nlon"
    v_diff@lat2d   = TLAT
    v_diff@lon2d   = TLONG

     

    
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


    ; gsn_define_colormap(wks, "BlueWhiteOrangeRed")

    ; gsn_define_colormap(wks, "BlRe")
    ; gsn_define_colormap(wks, "BlueDarkRed18")
 
 


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1
    
    ; res@cnMinLevelValF       =   -1.5
    ; res@cnMaxLevelValF       =   1.5
    ; res@cnLevelSpacingF      =  0.1
    
    ; res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,0,0,16,18,20,21/)

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    res@mpMinLatF             = -60.0
    res@mpMaxLatF             = 60.0
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

    res@tmXBLabelFontHeightF=0.01
    res@tmYLLabelFontHeightF=0.01
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0
    res@gsnLeftStringFontHeightF = 0.02

    
    res@gsnRightString = ""
    res@gsnLeftString = season(i) + " :SEN -CTRL"

      ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
 
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcMinMagnitudeF          =0.02
    res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 1.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    

 
    ; plot(i) = gsn_csm_contour_map(wks, PSL_diff_indian, res)
      plot(i) = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PSL_diff_indian, res)
    

    res2 = res
    res2@cnLevelSelectionMode = "ExplicitLevels"
    res2@cnFillOn = False
    res2@cnLinesOn = True
    res2@cnLevels = (/0.55/)
    res2@cnLineColors = "black"
    res2@cnMonoLineColor  = True
    res2@cnLineThicknessF = 2
    res2@cnMonoLineThickness  = True
    res2@gsnLeftString = ""
    res2@cnInfoLabelOn = False

  

    plot1(i) = gsn_csm_contour(wks, PSL_diff_indian,res2)
   
   
    overlay(plot(i), plot1(i))
 
     

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
    overlay(plot(i), iplot(i))
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
  

