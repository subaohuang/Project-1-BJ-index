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

    season = "MAM"
    fils1 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_spring_WNPCLM_spring/B2000_alt_north_spring_WNPCLM.cam.h0.0251-0280.nc"  
    fils2 = "/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_spring/B2000_alt_north_spring.cam.h0.0251-0280.nc"
    

  
    f1=addfile(fils1, "r")
  
    f2=addfile(fils2, "r")
    
    time1 = f1->time
    time2 = f2->time
    time1 = (/time1 -1/)
    time2 = (/time2 -1/)
    date1=cd_calendar(time1, 1)
    date2=cd_calendar(time2, 1)


    timeind1=ind(date1.ge.25101.and.date1.le.28012)
    timeind2=ind(date2.ge.25101.and.date2.le.28012)
    


    PRECL_CT  =  f1->PRECL(timeind1,:,:)
    PRECC_CT  =  f1->PRECC(timeind1,:,:)

    PRE_CT = PRECC_CT
    PRE_CT = (/PRECC_CT+PRECL_CT/)
    ; PRE_CT = (/PRECL_CT/)  ;;;大尺度降水
    PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

    PRECL_id  =  f2->PRECL(timeind2,:,:)
    PRECC_id  =  f2->PRECC(timeind2,:,:)

    PRE_id = PRECC_id
    PRE_id = (/PRECC_id+PRECL_id/)
    ; PRE_id = (/PRECL_id/)
    PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day

    PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season), 0)
    PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season), 0)

    PRE_diff = PRE_CT_clm
    PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

    ;;;t_test 
    prob_p = t_test( month_to_season(PRE_CT, season) ,month_to_season(PRE_id, season))
    PRE_diff@_FillValue = 10^9
     

    wks=gsn_open_wks("pdf", "diff-precipitation-restore-altantic-minus-WNPCLM-"+season+"-0280")
    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    ; res@cnFillMode  = 
    ; res@cnRasterSmoothingOn=True
    res@trGridType = "TriangularMesh"

    ; res@cnFillPalette = "hotcolr_19lev"


    gsn_define_colormap(wks, "MPL_BrBG")
 


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -0.5
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.1
    
    ; res@cnMinLevelValF       =   -1.0
    ; res@cnMaxLevelValF       =   1.0
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

    
    
    res@gsnRightString = ""
    res@gsnLeftString = "Indian SST forcing"
 
    plot = new(1,graphic)
    plot(0) = gsn_csm_contour_map(wks, PRE_diff, res)
     

    ;;;画box
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   
    lat1 = 0
    lat2 = 15
    lon1 = 110
    lon2 = 155

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    str3 = ispan(1, 3, 1) +""
    plot1 =new(1,graphic)
    ;;;
    
    plot1@$str3(0)$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres) 
   

    ;  ;;;打点
    ; res_pro=True
    ; res_pro@gsnAddCyclic = False
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnFillOn = True

    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLevels = (/0,0.1/)
    ; res_pro@lbLabelBarOn = False
    ; res_pro@cnLinesOn=False
    ; ; res_pro@cnMonoLineColor =False
    ; ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    ; res_pro@cnMonoFillPattern = True
    ; res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.002
  
    ; res_pro@cnMonoFillColor = False
    ; res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False
   
    ; res_pro@gsnCenterString = ""
    ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; ; res_pro@cnLineThicknessF = 3.
    ; ; res_pro@cnFillOn        = True
    ; iplot = new(1,graphic)

    ; iplot(0)=gsn_csm_contour(wks,prob_slp_indian,res_pro)
    

    
    ; overlay(plot(0), iplot(0))
    
    
   getvalues plot@contour
   "cnFillColors" :cnFillColors
   end getvalues
   print(cnFillColors)

   draw(plot)
   frame(wks)
    
    print("done script")
  end
  

