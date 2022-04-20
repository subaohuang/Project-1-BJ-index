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

  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012_chazhi.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/shanshan/F_2000/F_2000_IOBM_fivetimes/F_2000_IOBM_fivetimes.cam.h1.0101-4012_chazhi.nc"
  
  wks=gsn_open_wks("pdf", "U,V,Pre-IOBM_fivetimes")
  
  plot = new(4,graphic)
  iplot = new(4,graphic)
  str3 = ispan(1,4,1) + ""

  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")

  i = 0

  date=cd_calendar(f1->time, 0)  ;;; dimsizes(data)×6
  Yr = date(:,0)
  mon = date(:,1)

  ; print(Yr)
 
  timeind=ind(Yr.ge.6.and.Yr.le.30.and.(mon.ge.1.and.mon.le.4))

  PRECL_CT  =  f1->PRECL(timeind,:,:)
  PRECC_CT  =  f1->PRECC(timeind,:,:)

  PRE_CT = PRECC_CT
  PRE_CT = (/PRECC_CT+PRECL_CT/)
  ; PRE_CT = (/PRECL_CT/)  ;;;大尺度降水
  PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

  PRECL_id  =  f2->PRECL(timeind,:,:)
  PRECC_id  =  f2->PRECC(timeind,:,:)

  PRE_id = PRECC_id
  PRE_id = (/PRECC_id+PRECL_id/)
  ; PRE_id = (/PRECL_id/)
  PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day
   

  ;;;;U,V
  u_CT  =  f1->U(timeind,{850},:,:)
  v_CT  =  f1->V(timeind,{850},:,:)

  u_id  =  f2->U(timeind,{850},:,:)
  v_id  =  f2->V(timeind,{850},:,:)



  u_CT_clm = dim_avg_n_Wrap(u_CT, 0)
  u_indian_clm = dim_avg_n_Wrap(u_id, 0)

  u_diff = u_CT_clm
  u_diff = (/u_indian_clm - u_CT_clm/)

  v_CT_clm = dim_avg_n_Wrap(v_CT, 0)
  v_indian_clm = dim_avg_n_Wrap(v_id, 0)

  v_diff = v_CT_clm
  v_diff = (/v_indian_clm - v_CT_clm/)


  PRE_CT_clm = dim_avg_n_Wrap(PRE_CT, 0)
  PRE_indian_clm =dim_avg_n_Wrap(PRE_id, 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

  ; PRE_diff = (/PRE_diff*30/)


  ;;;t_test 
    prob_p = t_test( PRE_CT,PRE_id)
    prob_u = t_test( u_CT , u_id)
    prob_v = t_test( v_CT ,v_id)
  ;;;
   PRE_diff@_FillValue = 10^9
   u_diff@_FillValue = 10^9
   v_diff@_FillValue = 10^9
   ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
    ; u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
    ; v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
;;;
printMinMax(PRE_diff, 1)

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    ; res@cnLineColor = "White"
    ; res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    ; res@cnRasterSmoothingOn=True

    res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette="BlueRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -2
    res@cnMaxLevelValF       =   2
    res@cnLevelSpacingF      =  0.4
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = 25.0
    ; res@mpMaxLatF             = 50.0
    ; res@mpMinLonF             = 75.0
    ; res@mpMaxLonF             = 110.

    res@mpMinLatF             = -20.0
    res@mpMaxLatF             = 50.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 150.
    
    
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)

    ;;;;画省界

    res@mpOutlineOn   = True
    res@mpGeophysicalLineColor = "black"

    ; res@mpFillOn      = False
    res@mpDataBaseVersion = "MediumRes"
    res@mpDataSetName="Earth..4"
    res@mpOutlineSpecifiers=(/"China","Xinjiang Uygur"/)  ;NCL自带的地图轮廓，比较粗糙，边界划分失误严重
    ; res@mpProvincialLineColor="red" 
    ; res@mpProvincialLineThicknessF =4 



     ;;;LineArrow
    ; res@gsnScalarContour        = True
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcMinMagnitudeF          =0.02
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 3.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    
    
    res@gsnLeftString = "Precipitation DIFF (IOBM_fivetimes - Control)"
    res@gsnRightString = "mm/day"
    ; plot(i) = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)
    plot(i) = gsn_csm_contour_map(wks, PRE_diff, res)

  

   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
 
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
    ; plot@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

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
    res_pro@cnFillDotSizeF = 0.003
  
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
draw(plot(i))
frame(wks)
end
  

