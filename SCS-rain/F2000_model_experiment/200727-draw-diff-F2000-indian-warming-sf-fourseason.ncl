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
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h1.0101-4012_vp_sf.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindwarimg_" + season_r +"_notest/"+"F_2000_ind_"+season_r+"_notest.cam.h0.0101-4012_vp_sf.nc" 
  
  wks=gsn_open_wks("pdf", "CESM-F2000-sf-indian-warming-foul-season")
  
  plot = new(4,graphic)
  iplot = new(4,graphic)
  str3 = ispan(1,4,1) + ""


do i = 0,3
  f1=addfile(fils1, "r")
  f2=addfile(fils2(i), "r")
  print("读入file"+fils2(i))

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.4012)

  PRECL_CT  =  f1->sf(timeind,{850},:,:)

  PRE_CT = PRECL_CT

  PRECL_id  =  f2->sf(timeind,{850},:,:)


  PRE_id = PRECL_id

   
  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season(i)), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season(i)), 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

  PRE_diff = (/PRE_diff/10^6/)

;   ;;;;U,V
;   u_CT  =  f1->U(timeind,{850},:,:)
;   v_CT  =  f1->V(timeind,{850},:,:)

;   u_id  =  f2->U(timeind,{850},:,:)
;   v_id  =  f2->V(timeind,{850},:,:)


;   u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, season(i)), 0)
;   u_indian_clm = dim_avg_n_Wrap(month_to_season(u_id, season(i)), 0)

;   u_diff = u_CT_clm
;   u_diff = (/u_indian_clm - u_CT_clm/)

;   v_CT_clm = dim_avg_n_Wrap(month_to_season(v_CT, season(i)), 0)
;   v_indian_clm = dim_avg_n_Wrap(month_to_season(v_id, season(i)), 0)

;   v_diff = v_CT_clm
;   v_diff = (/v_indian_clm - v_CT_clm/)

;   ;;;t_test 
    prob_p = t_test( month_to_season(PRE_CT, season(i)) ,month_to_season(PRE_id, season(i)))
;     prob_u = t_test( month_to_season(u_CT, season(i)) ,month_to_season(u_id, season(i)))
;     prob_v = t_test( month_to_season(v_CT, season(i)) ,month_to_season(v_id, season(i)))
;   ;;;
;    PRE_diff@_FillValue = 10^9
;    u_diff@_FillValue = 10^9
;    v_diff@_FillValue = 10^9
;    ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
;     ; u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
;     ; v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
; ;;;


;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = False
    res@cnLinesOn =True
    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 3.0
    res@cnLineThicknesses = 2.0
    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=True
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "precip_diff_12lev"
      res@cnFillPalette = "BlRe"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -1.6
    res@cnMaxLevelValF       =   1.6
    res@cnLevelSpacingF      =  0.4
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -10.0
    ; res@mpMaxLatF             = 30.0
    ; ; res@mpMinLonF             = 40.0
    ; res@mpMinLonF             = 70.0
    ; res@mpMaxLonF             = 170.0
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 170.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
    ;;
    ;  ;;;LineArrow
    ; ; res@gsnScalarContour        = True
    ; res@gsnScalarContour        = False
    ; res@gsnLeftString = season(i)

    
    ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; ; res@vcGlyphStyle="CurlyVector"

    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; ; res@vcMinMagnitudeF          =0.02
    ; ; res@vcMaxMagnitudeF = 7.0
    ; res@vcRefMagnitudeF          = 2.0
    ; res@vcRefLengthF             = 0.05     ; define length of v
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ; res@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
    ; res@vcVectorDrawOrder="PostDraw"

    
    
    res@gsnRightString = ""
    ; plot(i) = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)
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
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.01     ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False


    ; iplot(i) = gsn_csm_contour(wks,prob_p,res_pro)

     
    

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; opt@gsnShadeMid     = 0         ; stipple
    ; opt@gsnshade

    ; ; opt@gsnShadeFillDotSizeF =  0.0025
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    ; iplot(i) = gsn_contour_shade(iplot(i) , 0,0.01,opt)
    ; overlay(plot(i), iplot(i))
    ; print("finish subset figure")

    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True
    res_pro@cnFillDrawOrder = "PreDraw"
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
    res_pro@lbLabelBarOn = False

    ; res_pro@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res_pro@cnMinLevelValF       =   0
    ; res_pro@cnMaxLevelValF       =   0.05
    ; res_pro@cnLevelSpacingF      =  0.05

    res_pro@cnLinesOn=True
    res_pro@cnMonoLineColor =False
    res_pro@cnLineColors = ((/"Transparent","gray","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 0
    ; res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","gray","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
   
     
    iplot(i) = gsn_csm_contour(wks,prob_p,res_pro)

      overlay(plot(i), iplot(i))


end do
  print("panel")
  pres = True
  pres@gsnPanelLabelBar = False

  gsn_panel(wks, plot, (/2,2/), pres)
  print("done script")
end
  

