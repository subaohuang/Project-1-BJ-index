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


  ; fils1 = "/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.TEMP.mon.0250-0299.nc"  
  fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/C/lsh_C_CTRL/lsh_C_CTRL.TEMP_TAUX_TAUY.h0.0001-0030.nc"
    ; fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h0.0251-0280.nc"
  fils2 = "/home/yangsong3/data-observation/linshh/CESM-data/C/lsh_C_addaltwind_spring/lsh_C_addaltwind_spring.TEMP_TAUX_TAUY.h0.0001-0030.nc"
 
  
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    
    time1 = f1->time
    time2 = f2->time
    time1 = (/time1 -1/)
    time2 = (/time2 -1/)
    date1=cd_calendar(time1, 1)
    date2=cd_calendar(time2, 1)
    
  

    timeind1=ind(date1.ge.501.and.date1.le.03012)

    ; timeind1=ind(date1.ge.25101.and.date1.le.28012)
    timeind2=ind(date2.ge.501.and.date2.le.03012)


    PSL_CT_r  =  f1->TEMP(timeind1,0,:,:)

  
    PSL_indian_r  =  f2->TEMP(timeind2,0,:,:)
    TLAT = f2->TLAT
    TLONG = f2->TLONG
   
    PSL_indian_r!1 = "nlat"
    PSL_indian_r!2 = "nlon"
    PSL_indian_r@lat2d   = TLAT
    PSL_indian_r@lon2d   = TLONG
    
    PSL_CT_r!1 = "nlat"
    PSL_CT_r!2 = "nlon"
    PSL_CT_r@lat2d   = TLAT
    PSL_CT_r@lon2d   = TLONG



   ;;;;差值
    PSL_CT = PopLatLon(PSL_CT_r,"gx1v6","1x1d","bilin","da","100716")

    PSL_indian = PopLatLon(PSL_indian_r,"gx1v6","1x1d","bilin","da","100716")
    

     ;;;;U,V
    u_CT_r  =  f1->TAUX(timeind1,:,:)
    v_CT_r  =  f1->TAUY(timeind1,:,:)


    u_CT_r!1 = "nlat"
    u_CT_r!2 = "nlon"
    u_CT_r@lat2d   = TLAT
    u_CT_r@lon2d   = TLONG


    v_CT_r!1 = "nlat"
    v_CT_r!2 = "nlon"
    v_CT_r@lat2d   = TLAT
    v_CT_r@lon2d   = TLONG

    
    u_id_r  =  f2->TAUX(timeind2,:,:)
    v_id_r  =  f2->TAUY(timeind2,:,:)

    u_id_r!1 = "nlat"
    u_id_r!2 = "nlon"
    u_id_r@lat2d   = TLAT
    u_id_r@lon2d   = TLONG


    v_id_r!1 = "nlat"
    v_id_r!2 = "nlon"
    v_id_r@lat2d   = TLAT
    v_id_r@lon2d   = TLONG
  
    ;;;;差值
    u_CT = PopLatLon(u_CT_r,"gx1v6","1x1d","bilin","da","100716")

    v_CT = PopLatLon(v_CT_r,"gx1v6","1x1d","bilin","da","100716")
    
    
    u_id = PopLatLon(u_id_r,"gx1v6","1x1d","bilin","da","100716")

    v_id = PopLatLon(v_id_r,"gx1v6","1x1d","bilin","da","100716")


    season = (/"MAM","JJA","SON","DJF"/)
    wks=gsn_open_wks("pdf", "201016-B2000-CLM-TAU-vector-by-addwind-POPmodel.pdf")

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
   
    ; prob_slp_indian!0 = "nlat"
    ; prob_slp_indian!1 = "nlon"
    ; prob_slp_indian@lat2d   = TLAT
    ; prob_slp_indian@lon2d   = TLONG
    ;;;;


   
    u_CT_clm = dim_avg_n_Wrap(u_CT, 0)
    u_indian_clm = dim_avg_n_Wrap(u_id, 0)

    u_diff = u_CT_clm
    u_diff = (/u_indian_clm - u_CT_clm/)

    ; u_diff = where(ismissing(u_diff), 0.0, u_diff)
    ; u_diff = where(u_diff.le.2.0, 0.0, u_diff)


    v_CT_clm = dim_avg_n_Wrap(v_CT, 0)
    v_indian_clm = dim_avg_n_Wrap(v_id, 0)

    v_diff = v_CT_clm
    v_diff = (/v_indian_clm - v_CT_clm/)

    ; v_diff = where(ismissing(v_diff), 0.0, v_diff)
    ; v_diff = where(v_diff.le.2.0, 0.0, v_diff)   
    
    ; u_diff!0 = "nlat"
    ; u_diff!1 = "nlon"
    ; u_diff@lat2d   = TLAT
    ; u_diff@lon2d   = TLONG

    
    ; v_diff!0 = "nlat"
    ; v_diff!1 = "nlon"
    ; v_diff@lat2d   = TLAT
    ; v_diff@lon2d   = TLONG

    
   




    
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
    ; res@cnMinLevelValF       =   -0.05
    ; res@cnMaxLevelValF       =   0.05
    ; res@cnLevelSpacingF      =  0.01
    
    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1
    
    ; res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,0,0,16,18,20,21/)

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    ; res@mpMinLatF             = -90.0
    ; res@mpMaxLatF             = 90.0
    ; res@mpMinLonF             = 0.0
    ; res@mpMaxLonF             = 360.0
    ; res@mpCenterLonF = 180.0

    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 290.0




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
    ; res@vcMinMagnitudeF          =0.2
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 0.05
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    

 
    ; plot(i) = gsn_csm_contour_map(wks, PSL_diff_indian, res)
     ; plot(i) = gsn_csm_contour_map(wks, v_diff, res)
    ; v_diff = 0.0
    ; plot(i) =  gsn_csm_vector_map(wks, u_diff_chazhi, v_diff_chazhi, res)
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
  

