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
season_r = (/"spring","summer","autumn","winter"/)
season = (/"MAM","JJA","SON","DJF"/)
  fils1 = systemfunc("ls /home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h0.00"+"*")

  
  wks=gsn_open_wks("pdf", "CESM-F2000-OLR-indian-warming-foul-season")
  
  plot = new(4,graphic)
  iplot = new(4,graphic)
  str3 = ispan(1,4,1) + ""

f1=addfiles(fils1, "r")


do i = 0,3

  fils2 = systemfunc("ls /home/ys17-19/lsh/CESM-data/F/F_2000_addindwarimg_" + season_r(i) +"_notest/"+"F_2000_ind_"+season_r(i)+"_notest.cam.h0.00"+"*")

  f2=addfiles(fils2, "r")


  date=cd_calendar(f1[:]->time, 1)
  timeind=ind(date.ge.601.and.date.le.4012)

  PRECL_CT  =  f1[:]->FLUT(timeind,:,:)


  PRE_CT = PRECL_CT
 
  PRECL_id  =  f2[:]->FLUT(timeind,:,:)
  delete(f2)

  PRE_id = PRECL_id

  
  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season(i)), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season(i)), 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)



  ;;;t_test 
    prob_p = t_test( month_to_season(PRE_CT, season(i)) ,month_to_season(PRE_id, season(i)))
  
    PRE_diff@_FillValue = 10^9




;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; ; res@cnFillPalette = "MPL_BrBG"

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res@cnLevels = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)*10.0
    ; ; gsn_define_colormap(wks, "precip_diff_12lev")
    ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
    ; ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
 

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    ; res@cnMinLevelValF       =   -8
    ; res@cnMaxLevelValF       =   8
    ; res@cnLevelSpacingF      =  1.0
    
    res@tmXBLabelFontHeightF=0.025
    res@tmYLLabelFontHeightF=0.025
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5.0
    res@mpMaxLatF             = 30.0
    ; res@mpMinLonF             = 40.0
    res@mpMinLonF             = 100.0
    res@mpMaxLonF             = 170.0
    ; res@mpMinLatF             = -60.0
    ; res@mpMaxLatF             = 60.0
    ; res@mpMinLonF             = 0.0
    ; res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
    ;;
     ;;;LineArrow
    ; res@gsnScalarContour        = True
    res@gsnScalarContour        = False
    res@gsnLeftString = season(i)+ " OLR diff("+season(i)+"_ind - CTRL)"
    res@gsnLeftStringFontHeightF  = 0.02


    
    
    res@gsnRightString = ""
   
    plot(i) = gsn_csm_contour_map(wks, PRE_diff,res)

  

   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
 
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
    plot@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

    ; ;;;打点
    ; res_pro=True
    ; res_pro@gsnAddCyclic = False
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnFillOn = True

    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLevels = (/0,0.01/)
    ; res_pro@lbLabelBarOn = False
    ; res_pro@cnLinesOn=False
    ; ; res_pro@cnMonoLineColor =False
    ; ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    ; res_pro@cnMonoFillPattern = True
    ; res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.003
  
    ; res_pro@cnMonoFillColor = False
    ; res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False
   
    ; res_pro@gsnCenterString = ""
    ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; ; res_pro@cnLineThicknessF = 3.
    ; ; res_pro@cnFillOn        = True
    ; iplot(i)=gsn_csm_contour(wks,prob_p,res_pro)
      ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.01     ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(i) = gsn_csm_contour(wks,prob_p,res_pro)

     
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.01,opt)
    overlay(plot(i), iplot(i))
    print("finish subset figure")
end do
  print("panel")
  pres = True
  pres@gsnPanelLabelBar = True
  pres@lbLabelFontHeightF = 0.015

  gsn_panel(wks, plot, (/2,2/), pres)
  print("done script")
end
  

