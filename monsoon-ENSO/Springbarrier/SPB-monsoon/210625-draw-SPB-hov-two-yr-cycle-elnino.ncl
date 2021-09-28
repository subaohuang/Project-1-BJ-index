begin 
  ; f1 = addfile("WNP_notbo_acr_hov_before_monsoon.nc", "r")
  ; f1 = addfile("ISM_acr_hov_before_monsoon.nc", "r")
  f1 = addfile("acr_hov_before_enino.nc", "r")

  spb_strmon_before = f1 ->spb_strmon 
  spb_weakmon_before = f1->spb_weakmon
  spb_normon_before = f1->spb_normon

  ; f2 = addfile("WNP_notbo_acr_hov_after_monsoon.nc", "r")
  ; f2 = addfile("ISM_acr_hov_after_monsoon.nc", "r")
  f2 = addfile("acr_hov_after_enino.nc", "r")

  spb_strmon_after = f2 ->spb_strmon 
  spb_weakmon_after = f2->spb_weakmon
  spb_normon_after = f2->spb_normon


  maind_strmon_before = spb_strmon_before@maind
  maind_weakmon_before = spb_weakmon_before@maind
  maind_normon_before = spb_normon_before@maind
  
  maind_strmon_after = spb_strmon_after@maind
  maind_weakmon_after = spb_weakmon_after@maind
  maind_normon_after = spb_normon_after@maind
  


  spb_strmon = new((/24,13/),float)
  spb_weakmon = new((/24,13/),float)
  spb_normon = new((/24,13/),float)



  spb_strmon(0:11,:) = (/spb_strmon_before/)
  spb_strmon(12:,:) = (/spb_strmon_after /)


  spb_weakmon(0:11,:) = (/spb_weakmon_before /)
  spb_weakmon(12:,:) = (/spb_weakmon_after /)


  spb_normon(0:11,:) = (/spb_normon_before /)
  spb_normon(12:,:) = (/spb_normon_after /)

  spb_diff = spb_strmon - spb_weakmon

  
  maind_strmon = new((/24/),float)
  maind_weakmon = new((/24/),float)
  maind_normon = new((/24/),float)

  maind_strmon(0:11) = (/maind_strmon_before/)
  maind_strmon(12:) = (/maind_strmon_after /)


  maind_weakmon(0:11) = (/maind_weakmon_before /)
  maind_weakmon(12:) = (/maind_weakmon_after /)


  maind_normon(0:11) = (/maind_normon_before /)
  maind_normon(12:) = (/maind_normon_after /)

  

;;;ACF center diff 
   acr_gra = center_finite_diff_n(spb_strmon, 2, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   ; SPB_b = dim_sum_n_Wrap(acr_gra_max(0:11), 0)
   SPB_str_a = dim_sum_n_Wrap(acr_gra_max(11:), 0)

   acr_gra = center_finite_diff_n(spb_weakmon, 2, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   SPB_weak_a = dim_sum_n_Wrap(acr_gra_max(11:), 0)   

   acr_gra = center_finite_diff_n(spb_normon, 2, False, 1, 1)
   acr_gra = -1*acr_gra  ;;;change to positive
   acr_gra_max = dim_max_n(acr_gra, 1) 
 
   SPB_nor_a = dim_sum_n_Wrap(acr_gra_max(11:), 0) 

   print(SPB_str_a)
   print(SPB_weak_a)
   print(SPB_nor_a)
exit()


  plot = new(3,graphic)

;;;;
;;;draw root mean square R(lagt) =  sqrt(((1/N)*sum(ACF(lagt)²))

   ;;分成 ：before and during monsoon   ： after monsoon
   N = 12
   
   ;;;before 
   acf_square_strmon_b = spb_strmon(0:11,:)^2
   acf_square_weakmon_b = spb_weakmon(0:11,:)^2
  
    ;;;after 
   acf_square_strmon_a = spb_strmon(11:,:)^2
   acf_square_weakmon_a = spb_weakmon(11:,:)^2
  


   R_strmon_a =  sqrt(dim_sum_n(acf_square_strmon_a(:,:),0)/N)

   R_weakmon_a =  sqrt(dim_sum_n(acf_square_weakmon_a(:,:),0)/N)

   R_strmon_b =  sqrt(dim_sum_n(acf_square_strmon_b(:,:),0)/N)

   R_weakmon_b =  sqrt(dim_sum_n(acf_square_weakmon_b(:,:),0)/N)

wks = gsn_open_wks("pdf", "twoyear-cycle-acr_hov-Enino")

wks1 = gsn_open_wks("pdf", "root_mean_square-Enino")


;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=True
    res@gsnAddCyclic = False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn= True
    ; res@pmLabelBarOrthogonalPosF = 0.2
    ; res@lbLabelStride=1
    res@lbOrientation="Vertical"
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmBorderThicknessF=1.0
    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF = 1.5

    res@tmXBValues = ispan(0,12,1)
    res@tmXBLabels = ispan(0,12,1) +""

    res@tmYLValues = ispan(0,23,2)
    res@tmYLMinorValues = ispan(0,23,1)


    res@tmYLLabels = (/"Dec(-1)","Feb(0)","Apr(0)","Jun(0)","Aug(0)","Oct(0)","Dec(0)","Feb(+1)","Apr(+1)","Jun(+1)","Aug(+1)","Oct(+1)"/) 

    ; res@cnFillPalette = "BlueDarkOrange18"
    gsn_define_colormap(wks,"BlueDarkRed18")


    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels

    res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 
  ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

  res@gsnLeftString = "EL Nino"
   plot(0) = gsn_csm_contour(wks,spb_strmon,res)
   res@gsnLeftString = "La Nina"
   plot(1) = gsn_csm_contour(wks,spb_weakmon,res)
   res@gsnLeftString = "EL Nino - La Nina"
   ; plot(2) = gsn_csm_contour(wks,spb_normon,res)
   plot(2) = gsn_csm_contour(wks,spb_diff,res)

   resp = True
   resp@gsMarkerColor = "black"
   resp@gsMarkerSizeF = 0.05
   resp@gsMarkerIndex = 1

   plot_gs = gsn_add_polymarker(wks, plot(0), maind_strmon, ispan(0,23,1), resp)
   plot_gs1 = gsn_add_polymarker(wks, plot(1), maind_weakmon, ispan(0,23,1), resp)
   ; plot_gs2 = gsn_add_polymarker(wks, plot(2), maind_normon, ispan(0,23,1), resp)



   pres = True

   gsn_panel(wks, plot, (/2,2/), pres)

   ;;;;;; draw root_mean_square

      resxy = True
      resxy@gsnDraw = False
      resxy@gsnFrame = False
      resxy@trXMinF = 0
      resxy@trXMaxF = 12

      resxy@trYMaxF = 1.0
      resxy@trYMinF = 0.0
      
      resxy@xyMonoLineColor = False
      resxy@xyLineColors = (/"red","blue"/)
      resxy@xyDashPatterns= (/0,1/)

      resxy@gsnLeftString = "Rb"
      resxy@vpHeightF = 0.4
      resxy@vpWidthF = 0.6
      plotxy = new(2,graphic)

      plotxy(0) = gsn_csm_xy(wks1, ispan(0,12,1), (/R_strmon_b,R_weakmon_b/) , resxy) 
      resxy@gsnLeftString = "Ra"
      plotxy(1) = gsn_csm_xy(wks1, ispan(0,12,1), (/R_strmon_a,R_weakmon_a/) , resxy) 
     
    
    gsn_panel(wks1, plotxy, (/2,1/), pres)
    
   end 
