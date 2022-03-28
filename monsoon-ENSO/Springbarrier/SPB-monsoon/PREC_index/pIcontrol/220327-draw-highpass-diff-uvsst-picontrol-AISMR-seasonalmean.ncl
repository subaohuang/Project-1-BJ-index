;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var = dim_variance_n(pv_ano(east_time,:,:), 0)
    var  =where(var.eq.0, var@_FillValue, var)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
    var@dim = n
 return var
 end
function composite_mean(ssta_month_12,u_850_month_12,v_850_month_12,ind1,ind2,ind3)
begin 

    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/))

    prob_ssta_str = one_test(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)
    prob_u_str = one_test(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    prob_v_str = one_test(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,u_month_str,v_month_str,prob_ssta_str,prob_u_str,prob_v_str/]
end 
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99712
    startyear = 1
    endyear = 997

    lat1 = -30
    lat2 = 30
    lon1 = -180
    lon2 = 180

; ; ; ;;;读入海温数据
;     f1 = addfile("ssta_picontrol_highpass_bw.nc", "r")
;     ssta = f1->ssta
; ;;;;;;

; ; ; ;;;;
;     f2 = addfile("uv_picontrol_highpass_bw.nc", "r")
;     u_850 = f2->u_850
;     v_850 = f2->v_850
; ;;;;
; ; ;;;读入海温数据
    f1 = addfile("ssta_picontrol_highpass.nc", "r")
    ssta = f1->ssta
;;;;;;

; ; ;;;;
    f2 = addfile("uv_picontrol_highpass.nc", "r")
    u_850 = f2->u_850
    v_850 = f2->v_850
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    year_out = addfile("piControl_monsoon_and_ensoyear_highpass.nc", "r") 
    ;;;; 03-997

    nino_year = year_out->nino_year 
    lanina_year = year_out->lanina_year 
    normal_year = year_out->normal_year 

    monstryear = year_out->monstryear 
    monwekyear = year_out->monwekyear
    noryear = year_out->noryear
 
    wet_warm_year =  year_out->wet_warm_year
    wet_cold_year = year_out->wet_cold_year
    wet_norm_year = year_out->wet_norm_year

    dry_warm_year = year_out->dry_warm_year
    dry_cold_year = year_out->dry_cold_year
    dry_norm_year = year_out->dry_norm_year

    nor_warm_year = year_out->nor_warm_year
    nor_cold_year = year_out->nor_cold_year
    nor_norm_year = year_out->nor_norm_year


    yearnum = endyear+2 - startyear +1 -1
    dim_sst = dimsizes(ssta)
    dim_u = dimsizes(u_850)


;;;;;
   sst_year = ispan(startyear,endyear+2,1) 
;;;;



    case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

    case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

    case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]



    case = [/case1,case2,case3/]


  wksstring = (/"ElNino","LaNina","Normal"/)


do k=0,2 

   case_need:= case[k]

   ind_ts1 :=  get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])
  

   u_year = ispan(startyear,endyear,1) 

   ind_us1_or :=  get1Dindex(u_year,case_need[0])
   ind_uw1_or :=  get1Dindex(u_year,case_need[1])
   ind_unor1 :=  get1Dindex(u_year,case_need[2])
   

   ind_uw1 :=  ind_uw1_or(0:dimsizes(ind_uw1_or)-2)  ;;去掉最后一年

   ind_us1 :=  ind_us1_or(0:dimsizes(ind_us1_or)-2)  ;;去掉最后一年

   
  plot = new(6,graphic)
  plot1 = new(6,graphic)
  plot2 = new(6,graphic)
  iplot = plot1
  iplot1 = plot
  iplot2 = plot
  ; wks = gsn_open_wks("pdf", "seasonalmean-WNP-UVSST_strong-monsoon-after")
  ; wks1 = gsn_open_wks("pdf", "seasonalmean-WNP-UVSST_weak-monsoon-after")

  wks = gsn_open_wks("pdf", "highpass-diff-UVSST-PIcontrol-seasonalmean-str-"+" "+wksstring(k))
  wks1 = gsn_open_wks("pdf", "highpass-diff-UVSST-PIcontrol-seasonalmean-weak-"+" "+wksstring(k))
  ; wks2 = gsn_open_wks("pdf", "UVPRE-PIcontrol-seasonalmean-normal-"+" "+wksstring(k))

 ; month = (/"JJA","SON","DJF(+1)","MAM(+1)"/)  
  month = (/"MAM","JJA","SON","DJF(+1)","MAM(+1)","JJA(+1)"/) 
 
 
;;;split to 12 month 
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 -1
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

  dim2  = dimsizes(u_850)
  yearnum := dim2(0)/12 -1
  u_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)
  v_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)


 do i=0,11   
     if i.le.9 then
       ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1-12:12,:,:)
       u_850_month_12(i,:,:,:) = u_850(i+2:dimsizes(u_850&time)-1-12:12,:,:)
       v_850_month_12(i,:,:,:) = v_850(i+2:dimsizes(v_850&time)-1-12:12,:,:)
     else
        ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1:12,:,:)
        u_850_month_12(i,:,:,:) = u_850(i+2:dimsizes(u_850&time)-1:12,:,:)
        v_850_month_12(i,:,:,:) = v_850(i+2:dimsizes(v_850&time)-1:12,:,:)
     end if
 end do 

do i=0,5
;;; strong and weak monsoon year ssta , u ,v 
   if i.le.3 then 
    print("i="+i)
    com_month_str  = composite_mean(ssta_month_12(i*3:i*3+2,:,:,:),u_850_month_12(i*3:i*3+2,:,:,:),v_850_month_12(i*3:i*3+2,:,:,:),ind_ts1,ind_us1,ind_us1)
    
    ssta_month_str = com_month_str[0]
    u_month_str = com_month_str[1]
    v_month_str = com_month_str[2]

    prob_ssta_str = com_month_str[3]
    prob_u_str = com_month_str[4]
    prob_v_str = com_month_str[5]

    com_month_weak  = composite_mean(ssta_month_12(i*3:i*3+2,:,:,:),u_850_month_12(i*3:i*3+2,:,:,:),v_850_month_12(i*3:i*3+2,:,:,:),ind_tw1,ind_uw1,ind_uw1)
    
    ssta_month_weak = com_month_weak[0]
    u_month_weak = com_month_weak[1]
    v_month_weak = com_month_weak[2]

    prob_ssta_weak = com_month_weak[3]
    prob_u_weak = com_month_weak[4]
    prob_v_weak = com_month_weak[5]

    com_month_nor  = composite_mean(ssta_month_12(i*3:i*3+2,:,:,:),u_850_month_12(i*3:i*3+2,:,:,:),v_850_month_12(i*3:i*3+2,:,:,:),ind_tnor1,ind_unor1,ind_unor1)
    
    ssta_month_nor = com_month_nor[0]
    u_month_nor = com_month_nor[1]
    v_month_nor = com_month_nor[2]

    prob_ssta_nor = com_month_nor[3]
    prob_u_nor = com_month_nor[4]
    prob_v_nor = com_month_nor[5]


    else 
    com_month_str  = composite_mean(ssta_month_12(i*3-12:i*3+2-12,:,:,:),u_850_month_12(i*3-12:i*3+2-12,:,:,:),v_850_month_12(i*3-12:i*3+2-12,:,:,:),ind_ts1+1,ind_us1+1,ind_us1+1)
    
    ssta_month_str = com_month_str[0]
    u_month_str = com_month_str[1]
    v_month_str = com_month_str[2]

    prob_ssta_str = com_month_str[3]
    prob_u_str = com_month_str[4]
    prob_v_str = com_month_str[5]

    com_month_weak  = composite_mean(ssta_month_12(i*3-12:i*3+2-12,:,:,:),u_850_month_12(i*3-12:i*3+2-12,:,:,:),v_850_month_12(i*3-12:i*3+2-12,:,:,:),ind_tw1+1,ind_uw1+1,ind_uw1+1)
    
    ssta_month_weak = com_month_weak[0]
    u_month_weak = com_month_weak[1]
    v_month_weak = com_month_weak[2]

    prob_ssta_weak = com_month_weak[3]
    prob_u_weak = com_month_weak[4]
    prob_v_weak = com_month_weak[5]

    com_month_nor  = composite_mean(ssta_month_12(i*3-12:i*3+2-12,:,:,:),u_850_month_12(i*3-12:i*3+2-12,:,:,:),v_850_month_12(i*3-12:i*3+2-12,:,:,:),ind_tnor1+1,ind_unor1+1,ind_unor1+1)
    
    ssta_month_nor = com_month_nor[0]
    u_month_nor = com_month_nor[1]
    v_month_nor = com_month_nor[2]

    prob_ssta_nor = com_month_nor[3]
    prob_u_nor = com_month_nor[4]
    prob_v_nor = com_month_nor[5]
    
    end if

    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
 

;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
; 
  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "precip_diff_12lev"
  ; res@cnFillPalette = "BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 300

  res@mpFillOn =False

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.01
  res@tmYLLabelFontHeightF=0.01
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.5
  res@cnMaxLevelValF       =   0.5
  res@cnLevelSpacingF      =  0.1

   ; res@cnFillColors= (/2,3,4,5,6,8,8,9,10,11,13,14/) - 2 

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; ; res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False
  
  res@gsnLeftString = month(i) + " SST&850UV " 
  res@gsnLeftStringFontHeightF = 0.015

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

   res@gsnScalarContour        = True

    ;;;;;;;LineArrow
    res@vcMinDistanceF           = 0.02           ; thin vectors
    res@vcMinMagnitudeF          =0.4
    res@vcRefMagnitudeF          = 4.0
    res@vcRefLengthF             = 0.1     ; define length of v
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF = 1.0
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.006

    if i.eq.5 then 
    res@vcRefAnnoOn = True
    res@vcRefAnnoOrthogonalPosF = 0.2
    else 
    res@vcRefAnnoOn = False
    end if
  
    res@gsnRightString = wksstring(k) + " event (Str Mon - Without Mon)"

    u_month_str = (/u_month_str - u_month_nor/)
    v_month_str = (/v_month_str - v_month_nor/)
    ssta_month_str = (/ssta_month_str - ssta_month_nor/)

    
    u_month_weak = (/u_month_weak - u_month_nor/)
    v_month_weak = (/v_month_weak - v_month_nor/)
    ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)


    var1 = (/prob_ssta_str/)
    var2 = (/prob_ssta_nor/)
    var3 = (/prob_ssta_weak/)

    s1 = prob_ssta_str@dim
    s2 = prob_ssta_nor@dim
    s3 = prob_ssta_weak@dim


    prob_ssta_str = (/ttest(ssta_month_str, var1,s1 , ssta_month_nor, var2, s2, False, False)/)
    prob_ssta_weak = (/ttest(ssta_month_weak, var3, s3, ssta_month_nor, var2, s2, False, False)/)

    copy_VarCoords(ssta_month_str,prob_ssta_str)
    copy_VarCoords(ssta_month_str,prob_ssta_weak)


    ; res@gsnRightString = "Strong ISM"
    plot(i) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, ssta_month_str, res)
       res@gsnRightString = wksstring(k) + " event (Weak Mon - Without Mon)"
    ; res@gsnRightString = "Weak ISM"

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -0.5
    ; res@cnMaxLevelValF       =   0.5
    ; res@cnLevelSpacingF      =  0.1

   plot1(i) = gsn_csm_vector_scalar_map(wks1,u_month_weak, v_month_weak, ssta_month_weak, res)
    ; res@gsnRightString = wksstring(k) + " event (Weak Mon - Without Mon)"

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -1.0
    ; res@cnMaxLevelValF       =   1.0
    ; res@cnLevelSpacingF      =  0.2

   ; plot2(i) = gsn_csm_vector_scalar_map(wks2,u_month_nor, v_month_nor, ssta_month_nor, res)

  
  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.0
  pgres@gsLineDashPattern=0
  ;;;画box
  lat1 = -5
  lat2 = 5
  lon1 = 190
  lon2 = 240

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;

   plot@$month(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   plot1@$month(i)$ = gsn_add_polyline(wks1, plot1(i), pgx_pool, pgy_pool, pgres)
   ; plot2@$month(i)$ = gsn_add_polyline(wks2, plot2(i), pgx_pool, pgy_pool, pgres)

   ; str1 = unique_string("dum1")  ; "unique_string" will return a unique    

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
    res_pro@cnFillDotSizeF = 0.001
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

  
    iplot(i) = gsn_csm_contour(wks,prob_ssta_str,res_pro)
    iplot1(i) = gsn_csm_contour(wks1,prob_ssta_weak,res_pro)
    ; iplot2(i) = gsn_csm_contour(wks2,prob_ssta_nor,res_pro)

    overlay(plot(i), iplot(i))
    overlay(plot1(i), iplot1(i))
    ; overlay(plot2(i), iplot2(i))
    delete(prob_ssta_str)
    delete(prob_ssta_weak)
end do
pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelStride=1
gsn_panel(wks, plot(0:3), (/2,2/), pres)
gsn_panel(wks1, plot1(0:3), (/2,2/), pres)
; gsn_panel(wks2, plot2, (/3,2/), pres)

end do

end 

