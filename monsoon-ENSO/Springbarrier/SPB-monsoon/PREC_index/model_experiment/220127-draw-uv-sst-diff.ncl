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
      ; simple_size =  avg1
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
  fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/Cam/B2000_f19g16_CTRL.cam.h1.0250-0300_chazhi.nc"
 
  year = "0264"
    
  fils2 = systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/Cam/B2000_f19g16_CTRL.cam.h0."+year+"*")
  f2 = addfiles(fils2, "r")
  
 
  wks=gsn_open_wks("pdf", year+"_DIFF-B2000-u,v,p")
  ; wks=gsn_open_wks("pdf", "DIFF-CESM-WNPIDLHeating-u,v,p")

  plot = new(4,graphic)
  iplot = new(4,graphic)
  str3 = ispan(1,4,1) + ""

  startmon = 25001
  endmon = 30012

  f1=addfile(fils1, "r")


  time1= f1->time  
  time1 = (/time1-1/)  
  date=cd_calendar(time1, 1)
  
  timeind=ind(date.ge.startmon.and.date.le.endmon)
  
  ; time2= f2->time 
  ; time2 = (/time2-1/)
  ; date2=cd_calendar(time2, 1)

  ; timeind2=ind(date2.ge.startmon.and.date2.le.endmon)

  PRECL_CT  =  f1->PRECL(timeind,:,:)
  PRECC_CT  =  f1->PRECC(timeind,:,:)

  PRE_CT = PRECC_CT
  PRE_CT = (/PRECC_CT+PRECL_CT/)
  ; PRE_CT = (/PRECL_CT/)  ;;;大尺度降水
  PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

  PRECL_id  =  f2[:]->PRECL
  PRECC_id  =  f2[:]->PRECC

  PRE_id = PRECC_id
  PRE_id = (/PRECC_id+PRECL_id/)
  ; PRE_id = (/PRECL_id/)
  PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day
   

  ;;;;U,V
  u_CT  =  f1->U(timeind,{850},:,:)
  v_CT  =  f1->V(timeind,{850},:,:)

  u_id  =  f2[:]->U(:,{850},:,:)
  v_id  =  f2[:]->V(:,{850},:,:)


   PRE_CT_clm = clmMonTLL(PRE_CT)
   u_CT_clm = clmMonTLL(u_CT)
   v_CT_clm = clmMonTLL(v_CT)



    PRE_diff = PRE_CT_clm
    PRE_diff = (/PRE_id - PRE_CT_clm/)  ;;;12,lat,lon

    u_diff = u_CT_clm
    u_diff = (/u_id - u_CT_clm/)

    v_diff = v_CT_clm
    v_diff = (/v_id - v_CT_clm/)
;;;  

;;; jjas mean 

    PRE_diff_jjas = dim_avg_n_Wrap(PRE_diff(5:8,:,:), 0)
    u_diff_jjas = dim_avg_n_Wrap(u_diff(5:8,:,:), 0)
    v_diff_jjas = dim_avg_n_Wrap(v_diff(5:8,:,:), 0)
;;;

  ; ;;;t_test 
  ;   prob_p = t_test( month_to_season(PRE_CT, season(i)) ,month_to_season(PRE_id, season(i)))
  ;   prob_u = t_test( month_to_season(u_CT, season(i)) ,month_to_season(u_id, season(i)))
  ;   prob_v = t_test( month_to_season(v_CT, season(i)) ,month_to_season(v_id, season(i)))
  ; ;;;

;   ;;;t_test 
;   prob_p = t_test( month_to_season(PRE_CT, season(i)), month_to_season(PRE_id, season(i)))
;   prob_u = t_test( month_to_season(u_CT, season(i)), month_to_season(u_id, season(i)))
;   prob_v = t_test( month_to_season(v_CT, season(i)), month_to_season(v_id, season(i)))
;   PRE_diff@_FillValue = 10^9
;   u_diff@_FillValue = 10^9
;   v_diff@_FillValue = 10^9
;   PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
;   u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
;   ; v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
; ;;;


;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn= True
    res@pmLabelBarOrthogonalPosF = 0.3

    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "MPL_BrBG"


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -6
    res@cnMaxLevelValF       =   6
    res@cnLevelSpacingF      =  1.0

    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0

    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -40.0
    ; res@mpMaxLatF             = 40.0
    ; ; res@mpMinLonF             = 40.0
    ; res@mpMinLonF             = 30.0
    ; res@mpMaxLonF             = 170.0
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 300.0
   
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
    res@gsnLeftString = "UV850&Pre DIFF"+year
    res@gsnLeftStringFontHeightF = 0.02

    res@gsnLeftStringOrthogonalPosF  = -0.2

    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.012           ; thin vectors
    res@vcMinMagnitudeF          =0.6
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 6.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"


    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "6 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.3              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013

    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False

    
    i=0
    
    res@gsnRightString = ""
    ; plot = gsn_csm_vec/tor_scalar_map(wks,u_diff,v_diff, PRE_diff,res)
    plot(i) = gsn_csm_vector_scalar_map(wks,u_diff_jjas,v_diff_jjas, PRE_diff_jjas,res)

    ; plot(i) = gsn_csm_contour_map(wks,PRE_diff,res)


     pgres=True
     pgres@gsLineColor="black"
     pgres@gsLineThicknessF=2.0
     pgres@gsLineDashPattern=0

     
     ;;;画box
     lat1 = 0
     lat2 = 20
     lon1 = 125
     lon2 = 160

     pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
     pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
     ;;;
      plot@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
      ; plot@$str3$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

    ;;;打点
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


    ; iplot(i) = gsn_csm_contour(wks,prob_p,res_pro)

     
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    ; iplot(i) = gsn_contour_shade(iplot(i) , 0,0.01,opt)
    ; overlay(plot(i), iplot(i))
    print("finish subset figure")
draw(plot(0))
frame(wks)
  ; print("panel")
  ; pres = True
  ; pres@gsnPanelLabelBar = True

  ; gsn_panel(wks, plot, (/2,2/), pres)
  ; print("done script")
end
  

