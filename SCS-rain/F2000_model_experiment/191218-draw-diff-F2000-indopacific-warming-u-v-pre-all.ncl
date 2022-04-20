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

function read_exp(f2,timeind,PRE_CT,u_CT,v_CT)
begin

  PRECL_id  =  f2->PRECL(timeind,:,:)
  PRECC_id  =  f2->PRECC(timeind,:,:)

  PRE_id = PRECC_id
  PRE_id = (/PRECC_id+PRECL_id/)
  PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day

  PRE_CT_clm = dim_avg_n_Wrap(PRE_CT, 0)
  PRE_indian_clm = dim_avg_n_Wrap(PRE_id, 0)

  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)


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
  ;;;t_test 
  prob_p =  t_test(PRE_CT,PRE_id)
  prob_u =  t_test(u_CT,u_id)
  prob_v =  t_test(v_CT,v_id)
  ;;;;
  my_list = NewList("lifo")
  ListPush(my_list, prob_v)
  ListPush(my_list, prob_u)
  ListPush(my_list, prob_p)
  ListPush(my_list, v_diff)
  ListPush(my_list, u_diff)
  ListPush(my_list, PRE_diff)
  ; (/PRE_diff,u_diff,v_diff,prob_p,prob_u,prob_v/)

return my_list
end 

begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"

  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012_chazhi.nc"  
  fils3 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addpacifictrend/F_2000_pacific.cam.h0.0101-4012_chazhi.nc"  
  fils4 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindopacific/F_2000_indopacific.cam.h0.0101-4012_chazhi.nc"  


  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")
  f3=addfile(fils3, "r")
  f4=addfile(fils4, "r")




  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.4012)

  PRECL_CT  =  f1->PRECL(timeind,:,:)
  PRECC_CT  =  f1->PRECC(timeind,:,:)

  PRE_CT = PRECC_CT
  PRE_CT = (/PRECC_CT+PRECL_CT/)
  PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

  ;;;;U,V
  u_CT  =  f1->U(timeind,{850},:,:)
  v_CT  =  f1->V(timeind,{850},:,:)
  ;;;;
;;;


;;;;get_diff 
    
   exp_indian =  read_exp(f2,timeind,PRE_CT,u_CT,v_CT)
   exp_pacific =  read_exp(f3,timeind,PRE_CT,u_CT,v_CT)
   exp_indopacific = read_exp(f4,timeind,PRE_CT,u_CT,v_CT)
   


   PRE_diff_indian = exp_indian[0]
   PRE_diff_pacific = exp_pacific[0]
   PRE_diff_indopacific = exp_indopacific[0]

   ; PRE_diff_indopacific = (/PRE_diff_indian + PRE_diff_pacific/)



   u_diff_indian = exp_indian[1]
   u_diff_pacific = exp_pacific[1]
   u_diff_indopacific = exp_indopacific[1]

   v_diff_indian = exp_indian[2]
   v_diff_pacific = exp_pacific[2]
   v_diff_indopacific = exp_indopacific[2]
;;;;;


  ;;;t_test 
    prob_p_indian = exp_indian[3]
    prob_u_indian = exp_indian[4]
    prob_v_indian = exp_indian[5]

    prob_p_pacific = exp_pacific[3]
    prob_u_pacific = exp_pacific[4]
    prob_v_pacific = exp_pacific[5]
    
    prob_p_indopacific = exp_indopacific[3]
    prob_u_indopacific = exp_indopacific[4]
    prob_v_indopacific = exp_indopacific[5]
  ;;;




  ;;;
   ; PRE_diff@_FillValue = 10^9
   ; u_diff@_FillValue = 10^9
   ; v_diff@_FillValue = 10^9
   ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
    ; u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
    ; v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)


   wks=gsn_open_wks("pdf", "CESM-F2000-u,v,p_pattern-in-do-all")

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
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
    res@cnLevelSpacingF      =  0.25
    

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
    res@gsnLeftString = "Indian SST forcing"

    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.015           ; thin vectors
    res@vcMinMagnitudeF          =0.1
    res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 2.0
    res@vcRefLengthF             = 0.04     ; define length of v
    res@vcRefAnnoOrthogonalPosF    = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    
    
    res@gsnRightString = ""
    plot = new(3,graphic)
    plot(0) = gsn_csm_vector_scalar_map(wks,u_diff_indian,v_diff_indian, PRE_diff_indian,res)
        res@gsnLeftString = "Pacific SST forcing"
    plot(1) = gsn_csm_vector_scalar_map(wks,u_diff_pacific,v_diff_pacific, PRE_diff_pacific,res)
        res@gsnLeftString = "Indo-Pacific SST forcing"
    plot(2) = gsn_csm_vector_scalar_map(wks,u_diff_indopacific,v_diff_indopacific, PRE_diff_indopacific,res)


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
     str2 = ispan(1, 3, 1)+""
     do i=0,2
      plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
     end do
   
    ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
   
    res_pro@gsnLeftString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(3,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_p_indian,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob_p_pacific,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob_p_indopacific,res_pro)

    do i = 0,2
    overlay(plot(i), iplot(i))
    end do
    pres = True
    pres@gsnPanelLabelBar = True
    gsn_panel(wks, plot, (/3,1/), pres)

    
    print("done script")
  end
  

