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
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"
  ; fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/CTRL.cam.h0.0101-4012_chazhi.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012_chazhi.nc"  
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addpacifictrend/F_2000_pacific.cam.h0.0101-4012_chazhi.nc"  
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindopacific/F_2000_indopacific.cam.h0.0101-4012_chazhi.nc"
  

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.601.and.date.le.4012)
    

    PSL_CT  =  f1->PSL(timeind,:,:)
    PSL_CT = (/PSL_CT/100/)
 

    u_CT  =  f1->U(timeind,{850},:,:)
    v_CT  =  f1->V(timeind,{850},:,:)

    
    PSL_indian  =  f2->PSL(timeind,:,:)
    PSL_indian = (/PSL_indian/100/)

    u_indian  =  f2->U(timeind,{850},:,:)
    v_indian  =  f2->V(timeind,{850},:,:)
    
    PSL_CT_clm = dim_avg_n_Wrap(PSL_CT, 0)
    PSL_indian_clm = dim_avg_n_Wrap(PSL_indian, 0)
    
    u_CT_clm = dim_avg_n_Wrap(u_CT, 0)
    u_indian_clm = dim_avg_n_Wrap(u_indian, 0)

    v_CT_clm = dim_avg_n_Wrap(v_CT, 0)
    v_indian_clm = dim_avg_n_Wrap(v_indian, 0)


    PSL_diff = PSL_CT_clm
    PSL_diff = (/PSL_indian_clm - PSL_CT_clm/)

    PSL_glbal = avg(dim_avg_n_Wrap(PSL_diff, (/1/)))  ;;;;为什么要去掉全球平均后才是想要的结果？可能CTRL有问题？
    

    ; PSL_diff = (/PSL_diff - PSL_glbal/)

    u_diff = u_CT_clm
    u_diff = (/u_indian_clm - u_CT_clm/)
    
    v_diff = v_CT_clm
    v_diff = (/v_indian_clm - v_CT_clm/)

    
    ;;;显著性检验
    prob_psl = t_test(PSL_indian,PSL_CT)
    prob_u = t_test(u_indian,u_CT)
    prob_v = t_test(v_indian,v_CT)

    
    u_diff = where(prob_u.le.0.05.or.prob_v.le.0.05, u_diff, u_diff@_FillValue)
     
    
   
     wks=gsn_open_wks("pdf", "CESM-F2000-slp_u_v_pattern-indopacific-SCS")
    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnMinLevelValF       =   -1
    ; res@cnMaxLevelValF       =   1
    ; res@cnLevelSpacingF      =  0.1
    
    res@cnMinLevelValF       =   -0.8
    res@cnMaxLevelValF       =   0.8
    res@cnLevelSpacingF      =  0.2

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
   
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 45
    res@mpMinLonF             =90
    res@mpMaxLonF             = 160
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


;
    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.025           ; thin vectors
    res@vcMinMagnitudeF          =0.04
    res@vcRefMagnitudeF          = 3.0
    res@vcMaxMagnitudeF = 6.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=2.0
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"
    
    
    res@gsnRightString = ""
    plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PSL_diff,res)

     ;;添加区域
   ;;warm pool
   pgres = True
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
   
    res_pro@gsnCenterString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(1,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_psl,res_pro)

    i = 0
    overlay(plot(i), iplot(i))

    draw(plot)
    frame(wks)
    print("done script")
  end
  

