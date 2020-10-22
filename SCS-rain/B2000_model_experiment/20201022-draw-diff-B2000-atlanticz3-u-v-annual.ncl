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
  fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.cam.h1.0251-0280_chazhi.nc"
  fils2 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_summer_fixmonth/B2000_alt_north_summer_fixmonth.cam.h1.0251-0280_chazhi.nc" 
  
  wks=gsn_open_wks("pdf", "CESM-B2000-Z3-alt-jjawarming-JJA")
  

  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")
  print("读入file"+fils2)

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.25601.and.date.le.28012)


  PRECC_CT  =  f1->Z3(timeind,{200},:,:)

  PRE_CT = PRECC_CT

  PRECC_id  =  f2->Z3(timeind,{200},:,:)
  
 
  PRE_id = PRECC_id

   
  ; PRE_CT_clm = dim_avg_n_Wrap(month_to_annual(PRE_CT, 1), 0)
  ; PRE_indian_clm =dim_avg_n_Wrap(month_to_annual(PRE_id, 1), 0)
  
  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, "JJA"), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, "JJA"), 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)



  ;;;;U,V
  u_CT  =  f1->U(timeind,{200},:,:)
  v_CT  =  f1->V(timeind,{200},:,:)

  u_id  =  f2->U(timeind,{200},:,:)
  v_id  =  f2->V(timeind,{200},:,:)


  ; u_CT_clm = dim_avg_n_Wrap(month_to_annual(u_CT, 1), 0)
  ; u_indian_clm = dim_avg_n_Wrap(month_to_annual(u_id, 1), 0)

  u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, "JJA"), 0)
  u_indian_clm = dim_avg_n_Wrap(month_to_season(u_id, "JJA"), 0)


  u_diff = u_CT_clm
  u_diff = (/u_indian_clm - u_CT_clm/)

  ; v_CT_clm = dim_avg_n_Wrap(month_to_annual(v_CT, 1), 0)
  ; v_indian_clm = dim_avg_n_Wrap(month_to_annual(v_id, 1), 0)

  v_CT_clm = dim_avg_n_Wrap(month_to_season(v_CT, "JJA"), 0)
  v_indian_clm = dim_avg_n_Wrap(month_to_season(v_id, "JJA"), 0)


  v_diff = v_CT_clm
  v_diff = (/v_indian_clm - v_CT_clm/)

  ; ;;;t_test 
  ;   prob_p = t_test( PRE_CT ,PRE_id)
  ;   prob_u = t_test(u_CT,u_id)
  ;   prob_v = t_test( v_CT ,v_id)
  ; ;;;

    prob_p = t_test( month_to_season(PRE_CT, "JJA")  ,month_to_season(PRE_id, "JJA") )
    prob_u = t_test( month_to_season(u_CT, "JJA")  ,month_to_season(u_id, "JJA") )
    prob_v = t_test( month_to_season(v_CT, "JJA")  ,month_to_season(v_id, "JJA") )

   PRE_diff@_FillValue = 10^9
   u_diff@_FillValue = 10^9
   v_diff@_FillValue = 10^9
   
   ; if 1.eq."JJA" then
     ; PRE_diff = where(prob_p.le.0.4, PRE_diff, PRE_diff@_FillValue)
   ; else 
    ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
   ; end if

;     u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
;     v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
; ;;;


;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =True

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0




    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "precip_diff_12lev"
      res@cnFillPalette = "BlueWhiteOrangeRed"
      ; res@cnFillPalette = "BlueDarkRed18"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -30
    res@cnMaxLevelValF       =   30
    res@cnLevelSpacingF      =  5
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -10.0
    ; res@mpMaxLatF             = 30.0
    ; ; res@mpMinLonF             = 40.0
    ; res@mpMinLonF             = 70.0
    ; res@mpMaxLonF             = 170.0
    res@mpMinLatF             = -90.0
    res@mpMaxLatF             = 90.0
    res@mpMinLonF             = -180.0
    res@mpMaxLonF             = 180.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True

    res@gsnLeftString = "GCM JJA alt warming induced 200hpa geopotensial Height and U,V diff (SEN - CTRL)"
  
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.5
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 5.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    
    
    res@gsnRightString = ""
    plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)
    ; plot(i) = gsn_csm_contour_map(wks, PRE_diff,res)

  

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
      str3 = "1"
      plot@$str3$ = gsn_add_polyline(wks, plot, pgx_pool, pgy_pool, pgres)

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


      res_pro@cnMonoFillPattern = False
      res_pro@cnFillPatterns = (/-1,17,-1/) 
      ; res_pro@cnFillPattern = 17
      res_pro@cnFillDotSizeF = 0.003

      res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
      res_pro@cnLineLabelsOn=False
      res_pro@tmYRMode = "Automatic"
      res_pro@cnInfoLabelOn=False

      res_pro@gsnCenterString = ""
      ; res_pro@cnLineLabelBackgroundColor="Transparent"
      ; res_pro@cnLineThicknessF = 3.
      ; res_pro@cnFillOn        = True
      iplot=gsn_csm_contour(wks,prob_p,res_pro)
      ; overlay(plot(i), iplot(i))
      print("finish subset figure")

      draw(plot)
      frame(wks)
end 

  

