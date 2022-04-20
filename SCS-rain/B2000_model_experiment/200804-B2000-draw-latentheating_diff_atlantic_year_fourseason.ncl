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


  fils1 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.cam.h1.0251-0280_OLR_heatflux.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year/B2000_alt_north_year.cam.h1.0251-0280_OLR_heatflux.nc"
 
  
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



   
;;;;;FSNS,FLNS,LHFLX,SHFLX
  LHFLX =  f1->LHFLX(timeind1,:,:)
  SHFLX =  f1->SHFLX(timeind1,:,:)
  FLNS = f1->FLNS(timeind1,:,:)
  FSNS = f1->FSNS(timeind1,:,:)

  PRE_CT = FSNS
  ; PRE_CT = (/FSNS - LHFLX - SHFLX - FLNS /)
  ; PRE_CT = (/FSNS  - SHFLX - FLNS /)
 

  LHFLX_id =  f2->LHFLX(timeind2,:,:)
  SHFLX_id =  f2->SHFLX(timeind2,:,:)
  FLNS_id = f2->FLNS(timeind2,:,:)
  FSNS_id = f2->FSNS(timeind2,:,:)

  PRE_id = FSNS_id
   
    
    season = (/"MAM","JJA","SON","DJF"/)
    wks=gsn_open_wks("pdf", "diff-solar-heating-flux-restore-altantic-year-fourseason-0280.pdf")

    plot = new(4,graphic)
    str3 = ispan(1, 4, 1) +""
    plot1 =new(4,graphic)
    iplot = new(4,graphic)


  do i=0,3


    PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season(i)), 0)
    PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season(i)), 0)

    PRE_diff = PRE_CT_clm
    PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

    ;;;t_test 
     prob_pre_indian = t_test(month_to_season(PRE_CT, season(i)),month_to_season(PRE_id, season(i)))
     copy_VarCoords(PRE_indian_clm, prob_pre_indian)
    
    
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


        ; gsn_define_colormap(wks, "MPL_BrBG")
        ; res@cnFillPalette = "precip_diff_12lev"
        ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
        ; res@cnLevels = (/-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0/)
        ; res@cnFillColors= (/4,5,6,8,11,12,1/) - 2
  
           
        res@tmXBLabelFontHeightF=0.02
        res@tmYLLabelFontHeightF=0.02
        res@tmXBMajorThicknessF=1.0
        res@tmYLMajorThicknessF=1.0
        res@tmBorderThicknessF=1.0
        res@gsnLeftStringFontHeightF = 0.02


        res@gsnRightString = ""
        res@gsnLeftString = season(i) + " :PRE SEN - CTRL"




        res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
        res@cnMinLevelValF       =   -8.0
        res@cnMaxLevelValF       =   8.0
        res@cnLevelSpacingF      =  1.0

        ; res@cnMinLevelValF       =   -1.5
        ; res@cnMaxLevelValF       =   1.5
        ; res@cnLevelSpacingF      =  0.5


        ; res@cnFillPalette = "precip_diff_12lev"
        ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
        ; res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
        ; ; gsn_define_colormap(wks, "precip_diff_12lev")
        ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
        ; ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
 

        ; res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,0,0,16,18,20,21/)

        res@tmXTOn=False
        res@tmYROn=False

        res@mpLimitMode           = "LatLon"
         
        ; res@mpMinLatF             = -60.0
        ; res@mpMaxLatF             = 60.0
        ; res@mpMinLonF             = 0.0
        ; res@mpMaxLonF             = 360.0
        ; res@mpCenterLonF = 180.0

        res@mpMinLatF             = -30.0
        res@mpMaxLatF             = 30.0
        ; res@mpMaxLatF             = 50.0
        res@mpMinLonF             = 30.0
        res@mpMaxLonF             = 170.0




        res@mpFillOn = False
        res@mpGeophysicalLineColor      = "Grey18"
        res@mpGeophysicalLineThicknessF = 1.5



        ; res@mpShapeMode = "FreeAspect"
        ; res@vpWidthF = 0.7
        ; res@vpHeightF = 0.3



        res@gsnRightString = ""
    

        plot(i) = gsn_csm_contour_map(wks, PRE_diff, res)
     

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
        res_pro@cnLevels = (/0,0.05/)
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


        iplot(i)=gsn_csm_contour(wks,prob_pre_indian,res_pro)

      overlay(plot(i), iplot(i))
end do 
    pres = True
    pres@gsnPanelLabelBar = True
    
    gsn_panel(wks, plot, (/2,2/), pres)
    print("done script")
  end
  