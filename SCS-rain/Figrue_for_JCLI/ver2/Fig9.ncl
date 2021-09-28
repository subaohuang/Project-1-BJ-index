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
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012_chazhi.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindwarimg_" + season_r +"_notest/"+"F_2000_ind_"+season_r+"_notest.cam.h1.0101-4012_chazhi.nc" 
   ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_ind_" + season_r +"_notest_wholeind/"+"F_2000_ind_"+season_r+"_notest_wholeind.cam.h1.0101-4012_chazhi.nc" 

  fils3 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012_OLR_heatflux.nc"
  fils4 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindwarimg_" + season_r +"_notest/"+"F_2000_ind_"+season_r+"_notest.cam.h1.0101-4012_OLR_heatflux.nc" 
    ; fils4 = "/home/ys17-19/lsh/CESM-data/F/F_2000_ind_" + season_r +"_notest_wholeind/"+"F_2000_ind_"+season_r+"_notest_wholeind.cam.h1.0101-4012_OLR_heatflux.nc" 
  
  wks=gsn_open_wks("pdf", "Fig9")
  
  plot = new(8,graphic)
  iplot = new(8,graphic)
  str3 = ispan(1,8,1) + ""

  order  = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/)
 


do i = 0,3
  f1=addfile(fils1, "r")
  f2=addfile(fils2(i), "r")
  print("读入file"+fils2(i))

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.4012)

  PRECL_CT  =  f1->PRECL(timeind,:,:)
  PRECC_CT  =  f1->PRECC(timeind,:,:)

  PRE_CT = PRECC_CT
  PRE_CT = (/PRECC_CT+PRECL_CT/)
  ; PRE_CT = (/PRECL_CC/)  ;;;大尺度降水
  PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

  PRECL_id  =  f2->PRECL(timeind,:,:)
  PRECC_id  =  f2->PRECC(timeind,:,:)

  PRE_id = PRECC_id
  PRE_id = (/PRECC_id+PRECL_id/)
  ; PRE_id = (/PRECC_id/)
  PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day
   
  PRE_CT_clm = dim_avg_n_Wrap(month_to_season(PRE_CT, season(i)), 0)
  PRE_indian_clm =dim_avg_n_Wrap(month_to_season(PRE_id, season(i)), 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

  ;;;;U,V
  u_CT  =  f1->U(timeind,{850},:,:)
  v_CT  =  f1->V(timeind,{850},:,:)

  u_id  =  f2->U(timeind,{850},:,:)
  v_id  =  f2->V(timeind,{850},:,:)


  u_CT_clm = dim_avg_n_Wrap(month_to_season(u_CT, season(i)), 0)
  u_indian_clm = dim_avg_n_Wrap(month_to_season(u_id, season(i)), 0)

  u_diff = u_CT_clm
  u_diff = (/u_indian_clm - u_CT_clm/)

  v_CT_clm = dim_avg_n_Wrap(month_to_season(v_CT, season(i)), 0)
  v_indian_clm = dim_avg_n_Wrap(month_to_season(v_id, season(i)), 0)

  v_diff = v_CT_clm
  v_diff = (/v_indian_clm - v_CT_clm/)

  ;;;t_test 
  prob_p = t_test( month_to_season(PRE_CT, season(i)) ,month_to_season(PRE_id, season(i)))
  prob_u = t_test( month_to_season(u_CT, season(i)) ,month_to_season(u_id, season(i)))
  prob_v = t_test( month_to_season(v_CT, season(i)) ,month_to_season(v_id, season(i)))
  ;;;
  PRE_diff@_FillValue = 10^9
  u_diff@_FillValue = 10^9
  v_diff@_FillValue = 10^9
  ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
  

  if i.eq.2 then 
    prob_u({10:30},{100:130}) = where(prob_u({10:30},{100:130}).le.1.0, 0.1, prob_u({10:30},{100:130})) 
  end if

  u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
  v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
;;;
    
;;;OLR
   f3 = addfile(fils3, "r")
   f4 = addfile(fils4(i), "r")

   OLR_CT = f3->FLUT(timeind,:,:)

   OLR_indian = f4->FLUT(timeind,:,:)

   OLR_CT_clm = dim_avg_n_Wrap(month_to_season(OLR_CT, season(i)), 0)
   OLR_indian_clm = dim_avg_n_Wrap(month_to_season(OLR_indian, season(i)), 0)

   OLR_diff = OLR_CT_clm
   OLR_diff = (/OLR_indian_clm - OLR_CT_clm/)

    region_mean = avg(OLR_diff({0:15},{125:160}))
    ;;;;;

    print("*********"+season(i)+"**********OLR_diff="+region_mean)
    
   
   prob_OLR = t_test( month_to_season(OLR_CT, season(i)) ,month_to_season(OLR_indian, season(i)))
   
   if i.eq.0 then 
   prob_OLR({0:15},{125:160}) = where(prob_OLR({0:15},{125:160}).le.0.01, 0.02, prob_OLR({0:15},{125:160}))
   end if 

   if i.eq.3 then 
   prob_OLR({0:15},{125:160}) = where(OLR_diff({0:15},{125:160}).ge.2.0, 0.005, prob_OLR({0:15},{125:160}))
   end if 


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

    if i.eq.3 then
    res@lbLabelBarOn= True
    else
    res@lbLabelBarOn= False
    end if
    res@pmLabelBarOrthogonalPosF = 0.2

    res@lbLabelStride=1
    res@lbLabelFontHeightF = 0.015
    res@lbAutoManage = False

    res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "MPL_BrBG"

    res@cnFillPalette = "precip_diff_12lev"
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ; res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
    ; res@cnLevels := (/-2.5,-2.,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2

    

    ; res@cnLevels := (/-1.5,-1.0,-0.25,0.25,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/2,3,4,8,9,10,11,12,13/) - 2 

    ; res@cnLevels := (/-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0/)
    ; res@cnFillColors= (/2,3,4,6,8,9,10,11,12,13/) - 2 
    
    res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
    res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5/)



    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF = 1.5
    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -5.0
    ; res@mpMaxLatF             = 30.0
    ; ; res@mpMinLonF             = 40.0
    ; res@mpMinLonF             = 100.0
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
     ;;;LineArrow
    res@gsnScalarContour        = True
    ; res@gsnScalarContour        = False

    res@gsnLeftString = order(i*2)+" "+ season(i)+ " Precip & 850UV diff(TIO_"+season(i)+" - CTRL)"
    res@gsnLeftStringOrthogonalPosF  = -0.12
    res@gsnLeftStringFontHeightF  = 0.02

    
    ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.3
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 2.0      
    res@vcRefLengthF             = 0.05     ; define length of v
   
    res@vcLineArrowThicknessF=1.2
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
    res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.006
    res@vcVectorDrawOrder="PostDraw"

    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "2 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013

    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False

  
    
    res@gsnRightString = ""
    plot(i*2) = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)
    

;;;;draw OLR result
 ; delete(res@cnLevelSelectionMode)
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillColors)
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "MPL_bwr"

  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  res@lbLabelStride=1

  res@cnLevels := (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)*20.0

   ; res@cnLevels := (/-20.0,-15.,-10.,-5.,5.0,10.,15.,20./)
  ; res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  ; res@cnMinLevelValF      = -2.5       ; set min contour level
  ; res@cnMaxLevelValF      = 2.5      ; set max contour level
  ; res@cnLevelSpacingF     = 0.5       ; set contour spacing

  res@gsnLeftString = order(i*2+1)+" "+ season(i)+ " OLR diff(TIO_"+season(i)+" - CTRL)"
    res@gsnLeftStringOrthogonalPosF  = 0.06


  plot(i*2+1) = gsn_csm_contour_map(wks, OLR_diff, res)

   
  ; 
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern= 1
 
   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   plot@$str3(i*2)$ = gsn_add_polyline(wks, plot(i*2), pgx_pool, pgy_pool, pgres)
   plot@$str3(i*2+1)$ = gsn_add_polyline(wks, plot(i*2+1), pgx_pool, pgy_pool, pgres)

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


    iplot(i*2) = gsn_csm_contour(wks,prob_p,res_pro)

    iplot(i*2+1) = gsn_csm_contour(wks,prob_OLR,res_pro)

     
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    iplot(i*2) = gsn_contour_shade(iplot(i*2) , 0,0.01,opt)
    iplot(i*2+1) = gsn_contour_shade(iplot(i*2+1) , 0,0.01,opt)
    overlay(plot(i*2), iplot(i*2))
    overlay(plot(i*2+1), iplot(i*2+1))
    print("finish subset figure")
end do
  print("panel")
  pres = True
  pres@gsnPanelLabelBar = False
  
   pres@gsnPanelBottom = 0.05

  gsn_panel(wks, plot, (/4,2/), pres)
  print("done script")
end
  

