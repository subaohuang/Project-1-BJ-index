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
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h1.0101-4012_OLR_heatflux.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addatlanticwarimg_" + season_r +"_notest/"+"F_2000_atl_"+season_r+"_notest.cam.h1.0101-4012_OLR_heatflux.nc" 
  
  wks=gsn_open_wks("pdf", "CESM-F2000-所有通量-alt-warming-foul-season-global")
  
  plot = new(4,graphic)
  iplot = new(4,graphic)
  str3 = ispan(1,4,1) + ""


do i = 0,3
  ; f1=addfile(fils2(0), "r")
  f1=addfile(fils1, "r")
  f2=addfile(fils2(i), "r")
  print("读入file"+fils2(i))

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.4012)

;;;;;FSNS,FLNS,LHFLX,SHFLX
  LHFLX =  f1->LHFLX(timeind,:,:)
  SHFLX =  f1->SHFLX(timeind,:,:)
  FLNS = f1->FLNS(timeind,:,:)
  FSNS = f1->FSNS(timeind,:,:)

  PRE_CT = SHFLX
  PRE_CT = (/FSNS - LHFLX - SHFLX - FLNS /)
 

  LHFLX_id =  f2->LHFLX(timeind,:,:)
  SHFLX_id =  f2->SHFLX(timeind,:,:)
  FLNS_id = f2->FLNS(timeind,:,:)
  FSNS_id = f2->FSNS(timeind,:,:)


  PRE_id =  SHFLX_id
  PRE_id = (/ FSNS_id -  LHFLX_id - SHFLX_id - FLNS_id/)

  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season(i)), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season(i)), 0)
  




  PRE_diff = PRE_indian_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

    if i.eq.0 
    PRE_diff_mam = PRE_diff
  end if
  
    PRE_diff = (/PRE_diff - PRE_diff_mam/)


  ; ;;;;U,V
  ; u_CT  =  f1->U(timeind,{850},:,:)

  ; u_id  =  f2->U(timeind,{850},:,:)



  ; u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, season(i)), 0)
  ; u_indian_clm = dim_avg_n_Wrap(month_to_season(u_id, season(i)), 0)

  ; u_diff = u_CT_clm
  ; u_diff = (/u_indian_clm - u_CT_clm/)

 

  ;;;t_test 
    prob_p = t_test( month_to_season(PRE_CT, season(i)) ,month_to_season(PRE_id, season(i)))
    ; prob_u = t_test( month_to_season(u_CT, season(i)) ,month_to_season(u_id, season(i)))
  
  ;;;
   PRE_diff@_FillValue = 10^9
   ; u_diff@_FillValue = 10^9

   
   ; if season(i).eq."JJA" then
     ; PRE_diff = where(prob_p.le.0.4, PRE_diff, PRE_diff@_FillValue)
   ; else 
    ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
   ; end if

    ; u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
    ; v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
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
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "precip_diff_12lev"
      res@cnFillPalette = "BlueWhiteOrangeRed"
      ; res@cnFillPalette = "matlab_jet"
      ; res@cnFillPalette = "BlueDarkRed18"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -20.0
    res@cnMaxLevelValF       =   20.0
    res@cnLevelSpacingF      =  5.0

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -60.0
    ; res@mpMaxLatF             = 60.0
    ; ; res@mpMinLonF             = 40.0
    ; res@mpMinLonF             = 0.0
    ; res@mpMaxLonF             = 360.0


    res@mpMinLatF             = -40.0
    res@mpMaxLatF             = 40.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 150.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2.0
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3
    ; ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
  
    res@gsnLeftString = season(i)

    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.02           ; thin vectors
    ; res@vcMinMagnitudeF          =0.02
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 1.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    
    
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
    iplot(i)=gsn_csm_contour(wks,prob_p,res_pro)
    ; overlay(plot(i), iplot(i))
    print("finish subset figure")
end do
  print("panel")
  pres = True
  pres@gsnPanelLabelBar = True

  gsn_panel(wks, plot, (/2,2/), pres)
  print("done script")
end
  

