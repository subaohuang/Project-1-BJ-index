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
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012_OLR_heatflux.nc"
  
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addind_wep_" + season_r +"/"+"F_2000_addind_wep_"+season_r+".cam.h1.0101-4012_OLR_heatflux.nc" 
  
  wks=gsn_open_wks("pdf", "CESM-F2000-WEP-ind-Warimg-fourseaon")
  
  plot = new(4,graphic)
  iplot = new(4,graphic)
  str3 = ispan(1,4,1) + ""


do i = 0,3
  f1=addfile(fils1, "r")
  f2=addfile(fils2(i), "r")

  date=cd_calendar(f1->time, 1)

  timeind=ind(date.ge.601.and.date.le.4012)
  ; print(date(timeind))
  PRE_CT  =  f1->FLUT(timeind,:,:)   ;;;OLR data
  
  PRE_id  =  f2->FLUT(timeind,:,:)

  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season(i)), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season(i)), 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)


  ;;;t_test 
    prob_p = t_test( month_to_season(PRE_CT, season(i)) ,month_to_season(PRE_id, season(i)))
  ;;;

  
  ;;;
   PRE_diff@_FillValue = 10^9
   ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
 
;;;


;;;;;;



    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    ; res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "BlueRed"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -10
    res@cnMaxLevelValF       =   10
    res@cnLevelSpacingF      =  2
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -5.0
    ; res@mpMaxLatF             = 30.0
    ; res@mpMinLonF             = 100.0
    ; res@mpMaxLonF             = 170.0
    res@mpMinLatF             = -60.0
    res@mpMaxLatF             = 60.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3
    ;;
  
    res@gsnLeftString = season(i)

  
    res@gsnRightString = ""
    plot(i) = gsn_csm_contour_map(wks, PRE_diff, res)

  

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

    ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.001
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
   
    res_pro@gsnCenterString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  
    iplot(i)=gsn_csm_contour(wks,prob_p,res_pro)
    overlay(plot(i), iplot(i))
end do

  pres = True
  pres@gsnPanelLabelBar = True

  gsn_panel(wks, plot, (/2,2/), pres)
  print("done script")
end
  

