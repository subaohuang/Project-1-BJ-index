;;*************************************************
;; calculate SPB intensity using the method in Jin et al 2021 by shuheng Lin 2021.06.16
;; we need to calculate ACF (autocorrlation function) firstly 
;; then calculate the gradient of ACF and then find the maximum ACF decline location 
;; more detailed can be found in Jin et al 2019,GRL
;;*************************************************************
function spb_cal(data,data2)  ;;  data:two dimension 12×yearnum array i.e.12×10
;;data1 : the monsoon anomous year  data2: next year
begin 
  dim = dimsizes(data)
  yearnum = dim(1)
  ;;; cal ACF, want to get a two dimensional array 12×13,initial month start June 
     acr = new((/12,13/),float)
     do i=0,11
         do j=0,12
            k=i+j 
            if k.le.11 then 
              acr(i,j) = escorc(data(i,0:yearnum-1),data(i+j,0:yearnum-1))
            else 
              acr(i,j) = escorc(data(i,0:yearnum-1),data2(i+j-12,0:yearnum-1))
            end if 
         end do 
     end do
   ;;;
 ;;;ACF center diff 
   acr_gra = center_finite_diff_n(acr, 1, False, 1, 1)
   ; acr_gra = acr_gra  ;;;change to positive
   ; acr_gra_max = dim_max_n(acr_gra, 1) 
 
   acr_gra!0 = "initial_month"
   acr_gra!1 = "lag_month"

   ; maind = dim_maxind(acr_gra(lag_month|:,initial_month|:),0)
   ; SPB = dim_sum_n_Wrap(acr_gra_max, 0)
   ; acr@maind = maind
   ; acr@SPB = SPB
  ;;
  acr_gra_mar2aug = new(12,float)

  do i=0,11
    if i.le.2 then 
    acr_gra_mar2aug(i) = sum(acr_gra(i,6-i+3:12))
    end if

    if i.le.9.and.i.gt.2 then 
    acr_gra_mar2aug(i) = sum(acr_gra(i,6-i+3:6-i+3+5))
    end if
    if i.gt.9 then 
    acr_gra_mar2aug(i) = sum(acr_gra(i,0:5+9-i))
    end if
  end do 
   SPB = dim_avg_n_Wrap(acr_gra_mar2aug, 0)
   acr_gra@SPB = SPB
  return acr_gra
  end 

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
;;;;
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99712
    startyear = 1
    endyear = 997
;;;读入海温数据
    diri="/home/yangsong3/CMIP6/1.0x1.0/CESM2/piControl/Omon/"
    fils1=systemfunc("ls "+diri+"tos_Omon_CESM2_piControl_r1i1p1f1_gn_0{0..9}"+"*")
  
    f1=addfiles(fils1, "r")
    date:=cd_calendar(f1[:]->time, 1)

    lat = f1[0]->lat
    lon = f1[0]->lon

    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    sst=lonFlip(f1[:]->tos(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    delete(ssta)
    nino34_djf = month_to_season(nino34, "DJF")

    printVarSummary(nino34_djf)

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)


    std_nino = dim_stddev(nino34_djf)
    sst_total_year = ispan(startyear, endyear+2, 1) ;;;0-199

    ; nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    ; lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))

    print("************** finish SST data*****************")
;;;



;;;PRE DATA 
    diri="/home/yangsong3/CMIP6/1.0x1.0/CESM2/piControl/"
    fils2=systemfunc("ls "+diri+"pr_Amon_CESM2_piControl_r1i1p1f1_gn_0{0..9}"+"*")

    f2=addfiles(fils2, "r")
    date:=cd_calendar(f2[:]->time, 1)
    lat = f2[0]->lat
    lon = f2[0]->lon
    
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    pre = lonFlip(f2[:]->pr(timeind,:,:))
  
    prea = pre
    
    prea = (/rmMonAnnCycTLL(pre)/)
    delete(pre)

    prea = (/dtrend_msg_n(ispan(1,dimsizes(prea&time),1),prea,False,True,0)/)

    prea_indian = dim_avg_n_Wrap(prea(:,{10:30},{70:90}), (/1,2/))
    delete(prea)

    prea_indian_12 = split_12(prea_indian)

    pre_jjasa = month_to_season(prea_indian, "MAM")
    delete(prea_indian)

    pre_jjasa = (/dim_avg_n_Wrap(prea_indian_12(5:8,:), 0)/)
    print("************** finish PRE data*****************")
;;;

   
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;

;;;strong and weak monsoon year 
   total_year = ispan(startyear,endyear,1)   ;;1-399


   std = dim_stddev_Wrap(DMI_ctrl_JJA)

   inds = ind(DMI_ctrl_JJA.ge.0.43*std)
   indw = ind(DMI_ctrl_JJA.le.-0.43*std)

   ; print(DMI_ctrl_clm)

   indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)
 
   monstryear = total_year(inds)
   monwekyear = total_year(indw)
   noryear = total_year(indnor)

   ; print(monstryear)
   ; print(monwekyear)
   ; print(noryear)
   ; ; exit()
;;;;;  


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal


 nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-3)  ;;;少取一年

 wet_warm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 wet_cold_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 dry_warm_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 dry_cold_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 dry_norm_year =  total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 
 nor_warm_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 nor_cold_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 nor_norm_year =  total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 print(wet_warm_year)
 print(wet_cold_year)
 print(wet_norm_year)
 
 print(dry_warm_year)
 print(dry_cold_year)
 print(dry_norm_year)

 print(nor_warm_year)
 print(nor_cold_year)
 print(nor_norm_year)

;;;split to 12 month nino34_12 is a array size (12,yearnum)
  yearnum = dimsizes(nino34)/12 -1
  nino34_12 = new((/12,yearnum/),float)

  do i=0,11   
    if i.le.6 then 
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1-12:12)
     else
     nino34_12(i,:) = nino34(i+5:dimsizes(nino34)-1:12)
    end if
  end do 
;;;;


; case1 = [/wet_warm_year,wet_cold_year,wet_norm_year/]

; case2 = [/dry_warm_year,dry_cold_year,dry_norm_year/]

; case3 = [/nor_warm_year,nor_cold_year,nor_norm_year/]

case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

case = [/case1,case2,case3/]


wksstring = (/"warm-ENSO-SPB-grad","cold-ENSO-SPB-grad","norm-ENSO-SPB-grad"/)

do i=0,2
;;; 
   case_need:= case[i]
   
   
   sst_year = ispan(startyear,endyear,1) ;; nino34_12 start from startyear-1 because nino34 start from startyear-2

   ind_s1 :=  get1Dindex(sst_year,case_need[0]);;;
   ind_w1 :=  get1Dindex(sst_year,case_need[1])
   ind_nor1 :=  get1Dindex(sst_year,case_need[2])
;;;;

;;; strong and weak monsoon year ssta 

  nino34_12_strmon := nino34_12(:,ind_s1)
  nino34_12_weakmon := nino34_12(:,ind_w1)
  nino34_12_normon := nino34_12(:,ind_nor1)

  nino34_12_strmon2 := nino34_12(:,ind_s1+1)
  nino34_12_weakmon2 := nino34_12(:,ind_w1+1)
  nino34_12_normon2 := nino34_12(:,ind_nor1+1)

  spb_strmon = spb_cal(nino34_12_strmon,nino34_12_strmon2)
  spb_weakmon = spb_cal(nino34_12_weakmon,nino34_12_weakmon2)
  spb_normon = spb_cal(nino34_12_normon,nino34_12_normon2)
  
   
  system("rm Picontrol_acr_hov_after_monsoon_"+wksstring(i)+".nc")
  fi_out = addfile("Picontrol_acr_hov_after_monsoon_"+wksstring(i)+".nc", "c")
  fi_out->spb_strmon = spb_strmon 
  fi_out->spb_weakmon = spb_weakmon
  fi_out->spb_normon = spb_normon

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

    res@tmYLValues = ispan(0,11,1)
    res@tmYLLabels = (/"Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May"/) 

    ; res@cnFillPalette = "BlueDarkOrange18"
    gsn_define_colormap(wks,"BlueDarkRed18")


    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels

    ; res@cnLevels = (/-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
    ; res@cnFillColors = (/3,6,0,12,13,14,15,16,17,18,19/) 

    res@cnLevels = (/-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05/)
    res@cnFillColors = (/2,4,5,8,9,10,11,13/) 



  ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
    res@lbLabelStride=1

   ; res@gsnLeftString = "warm year"
    res@gsnLeftString = "Strong ISMR"
   res@gsnRightString = ""+spb_strmon@SPB
   plot(0) = gsn_csm_contour(wks,spb_strmon,res)
   ; res@gsnLeftString = "cold year"
    res@gsnLeftString = "Weak ISMR"
   res@gsnRightString = ""+spb_weakmon@SPB
   plot(1) = gsn_csm_contour(wks,spb_weakmon,res)
   ; res@gsnLeftString = "normal year"
   res@gsnLeftString = "Normal ISMR"
   res@gsnRightString = ""+spb_normon@SPB
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
   y2 = ispan(2,11,1) ;;; line for mar
   

   resp = True
   resp@gsLineColor = "black"
   resp@gsLineThicknessF = 3.0
   
   


   plot_gs = gsn_add_polyline(wks, plot(0), x1(::-1), y1, resp)
   plot_gs1 = gsn_add_polyline(wks, plot(1),  x1(::-1), y1, resp)
   plot_gs2 = gsn_add_polyline(wks, plot(2),  x1(::-1), y1, resp)

   plot_gs_aug = gsn_add_polyline(wks, plot(0), x2(::-1), y2, resp)
   plot_gs_aug1 = gsn_add_polyline(wks, plot(1),  x2(::-1), y2, resp)
   plot_gs_aug2 = gsn_add_polyline(wks, plot(2),  x2(::-1), y2, resp)


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
    ; frame(wks)
   end do 
   end 
