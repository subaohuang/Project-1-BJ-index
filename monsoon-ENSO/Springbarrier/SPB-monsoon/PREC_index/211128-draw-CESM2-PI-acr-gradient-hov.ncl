;;*************************************************
;; main program 
;;*************************************************************
begin
   

wksstring = (/"warm-ENSO-SPB-grad","cold-ENSO-SPB-grad","norm-ENSO-SPB-grad"/)

do i=0,2
;;; 
 
  fi_out := addfile("Picontrol_acr_hov_after_monsoon_"+wksstring(i)+".nc", "r")
  spb_strmon = fi_out->spb_strmon 
  spb_weakmon = fi_out->spb_weakmon
  spb_normon = fi_out->spb_normon

; exit()
 print("strongmonsoon="+spb_strmon@SPB)
 print("weakmonsoon="+spb_weakmon@SPB)
 print("normalmonsoon="+spb_normon@SPB)
; exit()
  plot = new(3,graphic)

;;;;
;;;draw root mean square R(lagt) =  sqrt(((1/N)*sum(ACF(lagt)²))

   ;;分成3-8月 ：before and during monsoon   9-2月： after monsoon
   
   ; N = 12

   acf_square_strmon = spb_strmon^2
   acf_square_weakmon = spb_weakmon^2
   acf_square_normon = spb_normon^2
   N = dimsizes(acf_square_strmon(4:,1))
   R_strmon =  sqrt(dim_sum_n(acf_square_strmon(4:,:),0)/N) ;;;9月后

   R_weakmon =  sqrt(dim_sum_n(acf_square_weakmon(4:,:),0)/N)

   R_normon =  sqrt(dim_sum_n(acf_square_normon(4:,:),0)/N)

    
   P1_3_str =  sum(R_strmon(1:3))/(3-1+1)
   P4_10_str =  sum(R_strmon(4:10))/(10-4+1)
   P1_10_str =  sum(R_strmon(1:10))/(10-1+1)

   P1_3_weak =  sum(R_weakmon(1:3))/(3-1+1)
   P4_10_weak =  sum(R_weakmon(4:10))/(10-4+1)
   P1_10_weak =  sum(R_weakmon(1:10))/(10-1+1)
   
   P1_3_nor =  sum(R_normon(1:3))/(3-1+1)
   P4_10_nor =  sum(R_normon(4:10))/(10-4+1)
   P1_10_nor =  sum(R_normon(1:10))/(10-1+1)


print("P1_3_str="+P1_3_str)
print("P4_10_str="+P4_10_str)
print("P1_10_str="+P1_10_str)

print("P1_3_weak="+P1_3_weak)
print("P4_10_weak="+P4_10_weak)
print("P1_10_weak="+P1_10_weak)

print("P1_3_nor="+P1_3_nor)
print("P4_10_nor="+P4_10_nor)
print("P1_10_nor="+P1_10_nor)


; wks = gsn_open_wks("pdf", "WNP-after-acr_hov-strong-weak-monsoon")
; wks1 = gsn_open_wks("pdf", "WNP-after-root_mean_square-newmonsoon.pdf")

wks = gsn_open_wks("pdf", "CESM2-PIcontrol-ISM-"+wksstring(i))

;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=False
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

    res@tmYLValues = ispan(0,11,1)
    res@tmYLLabels = (/"Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May"/) 

    ; res@cnFillPalette = "BlueDarkOrange18"
    ; gsn_define_colormap(wks,"BlueDarkRed18")
    gsn_define_colormap(wks,"amwg256")
    gsn_define_colormap(wks,"BlAqGrYeOrReVi200")



    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels

    ; res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    ; res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 

    res@cnLevels = (/-0.25,-0.2,-0.15,-0.1,-0.05,0.0/)
    ; res@cnFillColors = (/2,4,5,8,9,10,11,13/) 
    res@cnFillColors = (/200,185,176,161,145,0,0/) 
    ; res@cnFillColors = (/17,33,49,65,81,0,0,129/) 



  ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

   ; res@gsnLeftString = "warm year"
    res@gsnLeftString = "Strong ISMR"
   res@gsnRightString = ""+spb_strmon@SPB*-1
   plot(0) = gsn_csm_contour(wks,spb_strmon,res)
   ; res@gsnLeftString = "cold year"
    res@gsnLeftString = "Weak ISMR"
   res@gsnRightString = ""+spb_weakmon@SPB*-1
   plot(1) = gsn_csm_contour(wks,spb_weakmon,res)
   ; res@gsnLeftString = "normal year"
   res@gsnLeftString = "Normal ISMR"
   res@gsnRightString = ""+spb_normon@SPB*-1
   plot(2) = gsn_csm_contour(wks,spb_normon,res)

   ; resp = True
   ; resp@gsMarkerColor = "black"
   ; resp@gsMarkerSizeF = 0.1
   ; resp@gsMarkerIndex = 1

   ; plot_gs = gsn_add_polymarker(wks, plot(0), spb_strmon@maind, ispan(0,11,1), resp)
   ; plot_gs1 = gsn_add_polymarker(wks, plot(1), spb_weakmon@maind, ispan(0,11,1), resp)
   ; plot_gs2 = gsn_add_polymarker(wks, plot(2), spb_normon@maind, ispan(0,11,1), resp)
   x1 = ispan(0,9,1)
   y1 = x1 ;;; line for mar

   x2 = ispan(3,12,1)
   y2 = ispan(2,11,1) ;;; line for aug
   
   x3 = ispan(0,2,1)
   y3 = ispan(0,2,1) ;;; line for aug


   resp = True
   resp@gsLineColor = "green"
   resp@gsLineThicknessF = 6.0
   
   


   plot_gs = gsn_add_polyline(wks, plot(0), x1(::-1), y1, resp)
   plot_gs1 = gsn_add_polyline(wks, plot(1),  x1(::-1), y1, resp)
   plot_gs2 = gsn_add_polyline(wks, plot(2),  x1(::-1), y1, resp)
    resp@gsLineColor = "blue"
   plot_gs_aug = gsn_add_polyline(wks, plot(0), x2(::-1), y2, resp)
   plot_gs_aug1 = gsn_add_polyline(wks, plot(1),  x2(::-1), y2, resp)
   plot_gs_aug2 = gsn_add_polyline(wks, plot(2),  x2(::-1), y2, resp)


   plot_gs_augs = gsn_add_polyline(wks, plot(0), x3, y3(::-1), resp)
   plot_gs_augs1 = gsn_add_polyline(wks, plot(1),  x3, y3(::-1), resp)
   plot_gs_augs2 = gsn_add_polyline(wks, plot(2),  x3, y3(::-1), resp)

   rest = True
   rest@txFontColor = "black"
   rest@txFontHeightF = 0.02
   
    plot_gt = gsn_add_text(wks, plot(0),"MAR", 5 ,5, rest)
    plot_gt1 = gsn_add_text(wks, plot(1),"MAR", 5 ,5, rest)
    plot_gt2 = gsn_add_text(wks, plot(2),"MAR", 5 ,5, rest)
  
    
    plot_gx = gsn_add_text(wks, plot(0),"AUG", 8.2 ,6.5, rest)
    plot_gx1 = gsn_add_text(wks, plot(1),"AUG", 8.2 ,6.5, rest)
    plot_gx2 = gsn_add_text(wks, plot(2),"AUG", 8.2 ,6.5, rest)
    

    plot_gxs = gsn_add_text(wks, plot(0),"AUG", 2 ,1, rest)
    plot_gxs1 = gsn_add_text(wks, plot(1),"AUG", 2 ,1, rest)
    plot_gxs2 = gsn_add_text(wks, plot(2),"AUG", 2 ,1, rest)

   pres = True

   gsn_panel(wks, plot, (/2,2/), pres)
   ; draw(plot)
   ; frame(wks)
   ;;;;;;; draw root_mean_square

    ;   resxy = True
    ;   resxy@gsnDraw = False
    ;   resxy@gsnFrame = False
    ;   resxy@trXMinF = 0
    ;   resxy@trXMaxF = 12

    ;   resxy@trYMaxF = 1.0
    ;   resxy@trYMinF = 0.0
      
    ;   resxy@xyMonoLineColor = False
    ;   resxy@xyLineColors = (/"red","blue","black"/)
    ;   resxy@xyDashPatterns= (/0,0,0/)
    ;   resxy@xyMonoLineThickness = True
    ;   resxy@xyLineThicknessF = 1.5

    ;   resxy@gsnLeftString = "R "
    ;   resxy@vpHeightF = 0.4
    ;   resxy@vpWidthF = 0.6
    ;   plotxy = new(1,graphic)
    ;   plotxy(0) = gsn_csm_xy(wks1, ispan(0,12,1), (/R_strmon,R_weakmon,R_normon/) , resxy) 
    ;   ;;添加legend 
    ;   gres = True
    ;   gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ;     ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    ;   gres@XPosPercent = 60    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    ;   lineres = True
    ;   ; lineres@lgLineColors = (/"black"，/) ; line colors
    ;   lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
    ;   lineres@lgLineThicknesses = resxy@xyLineThicknessF                       ; line thicknesses
    ;   lineres@LineLengthPercent = 9.
    ;   lineres@lgDashIndexes   = resxy@xyDashPatterns                     ; expressed as %, 0->100, length of line

    ;   textres = True
    ;   textres@lgLabels = (/"Str monsoon","Weak monsoon","Normal monsoon"/)
    ;   ; textres@lgLineLabelFontHeights = (/0.05,
    ;   ; textres@lgLabelOpacity = (/0/) 

    ;   plotxy(0) = simple_legend(wks1,plotxy(0),gres,lineres,textres)
     
    ; draw(plotxy)
    frame(wks)
   end do 
   end 
