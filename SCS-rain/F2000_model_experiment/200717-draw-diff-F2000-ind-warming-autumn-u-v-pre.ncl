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

  season  = "SON"

  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012_chazhi.nc"  
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindwarimg_autumn/F_2000_ind_autumn.cam.h0.0101-4012_chazhi.nc"  


  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.4012)

  PRECL_CT  =  f1->PRECL(timeind,:,:)
  PRECC_CT  =  f1->PRECC(timeind,:,:)

  PRE_CT = PRECC_CT
  PRE_CT = (/PRECC_CT+PRECL_CT/)
  PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

  PRECL_id  =  f2->PRECL(timeind,:,:)
  PRECC_id  =  f2->PRECC(timeind,:,:)

  PRE_id = PRECC_id
  PRE_id = (/PRECC_id+PRECL_id/)
  PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day
   

  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season), 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

  ;;;;U,V
  u_CT  =  f1->U(timeind,{850},:,:)
  v_CT  =  f1->V(timeind,{850},:,:)

  u_id  =  f2->U(timeind,{850},:,:)
  v_id  =  f2->V(timeind,{850},:,:)


  u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, season), 0)
  u_indian_clm = dim_avg_n_Wrap(month_to_season(u_id, season), 0)

  u_diff = u_CT_clm
  u_diff = (/u_indian_clm - u_CT_clm/)

  v_CT_clm = dim_avg_n_Wrap(month_to_season(v_CT, season), 0)
  v_indian_clm = dim_avg_n_Wrap(month_to_season(v_id, season), 0)

  v_diff = v_CT_clm
  v_diff = (/v_indian_clm - v_CT_clm/)

  ;;;t_test 
    prob_p = t_test( month_to_season(PRE_CT, season) ,month_to_season(PRE_id, season))
    prob_u = t_test( month_to_season(u_CT, season) ,month_to_season(u_id, season))
    prob_v = t_test( month_to_season(v_CT, season) ,month_to_season(v_id, season))
  ;;;
   PRE_diff@_FillValue = 10^9
   u_diff@_FillValue = 10^9
   v_diff@_FillValue = 10^9
   ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
    ; u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
    ; v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)


   wks=gsn_open_wks("pdf", "CESM-F2000-u,v,p_pattern-indian-"+season)

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =True
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "precip_diff_12lev"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -2
    res@cnMaxLevelValF       =   2
    res@cnLevelSpacingF      =  0.2
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 160.0

    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3
    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    res@gsnLeftString = "Atlantic SST forcing"

    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcMinMagnitudeF          =0.02
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 2.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    
    
    res@gsnRightString = ""
    plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)


   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
    plot@$str3$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)

   
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
    res_pro@cnFillDotSizeF = 0.004
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
   
    res_pro@gsnCenterString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(1,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_p,res_pro)

    i = 0
    overlay(plot(i), iplot(i))

    draw(plot)
    frame(wks)
    print("done script")
  end
  

