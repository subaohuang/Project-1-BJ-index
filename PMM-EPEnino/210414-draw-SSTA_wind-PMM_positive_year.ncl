load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;;**************************************************************
;;;main program
;;***************************************************************
;;;*********************************************************
function one_test(pv_ano,pv_clm,east_time)
begin
    pv_east=dim_avg_n_Wrap( pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east- pv_clm)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
 end

begin
 ;;;time span
    startmon =27601
    endmon = 32512

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;
 
;;;;model name 
   modelname = "B2000_f19g16_CTRL"
;;;;


;;;;读入控制实验结果
    diri="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/"
    fils1=systemfunc("ls "+diri+"B2000_f19g16_CP_CTRL_TEMP.h1.0251-0350.nc")
    f1=addfile(fils1, "r")
    time = f1->time
    time = (/time - 1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)
   
    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})
;;;;


;;;u,v 

    diri2="/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/"
    fils2=systemfunc("ls "+diri+"B2000_f19g16_CP_CTRL.cam.h1.0251-0340_chazhi.nc")
    f2=addfile(fils2, "r")
    time := f2->time
    time = (/time - 1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    uwnd = f2->U(timeind,0,:,:)
 
    vwnd = f2->V(timeind,0,:,:)

    printMinMax(uwnd, 1)

    uwnd = (/rmMonAnnCycTLL(uwnd)/)
    vwnd = (/rmMonAnnCycTLL(vwnd)/)

;;;;read PMM index  
  ; pmmindex=asciiread("/home/ys17-19/lsh/Project/PMM-EPEnino/pmm_index_CTRL_276-325.txt", -1, "float")
  pmmindex=asciiread("/home/ys17-19/lsh/Project/PMM-EPEnino/pmm_index_CTRL.txt", -1, "float")


  printVarSummary(pmmindex)
  pmmindex = pmmindex
;;;

  pmmindex!0= "time"

;;前期春季的PMM index 3,4 5
  pmmindex_mam = month_to_season(pmmindex, "MAM")
;;;

  pmm_posi_ind = ind(pmmindex_mam.ge.1.0)
  
;;; sst tendency 
  sst = (/rmMonAnnCycTLL(sst)/)
  printVarSummary(sst)

  cti = dim_avg_n_Wrap(sst(:,{-5:5},{180:270}), (/1,2/))
  cti_3d = conform(sst, cti, 0)

  reg_cti_sst = regCoef_n(cti, sst, 0, 0)
  reg_cti_sst_3d = conform(sst, reg_cti_sst, (/1,2/))
  sst_anom_rmcti = sst - cti_3d * reg_cti_sst_3d
  copy_VarCoords(sst, sst_anom_rmcti)

  sst_tendency = center_finite_diff_n(sst, 1, False, 0, 0)
  copy_VarCoords(sst, sst_tendency)


  sst_tendency_nino1_2 = dim_avg_n_Wrap(sst_tendency(:,{-10:0},{270:280}), (/1,2/))

  sst_tendency_nino1_2_jja = month_to_season(sst_tendency_nino1_2, "JJA")

  sst_tendency_std = dim_stddev_n(sst_tendency_nino1_2_jja, 0)

  sst_tendency_jja_ind = ind(sst_tendency_nino1_2_jja.ge.sst_tendency_std)
  print(sst_tendency_jja_ind)
  print(pmm_posi_ind)
  
  ind_overlap = sst_tendency_jja_ind
  ind_overlap = 0
  do i=0,dimsizes(sst_tendency_jja_ind)-1
     if any(sst_tendency_jja_ind(i).eq.pmm_posi_ind) then 
       ind_overlap(i) = sst_tendency_jja_ind(i)
    end if
  end do  
  ind_overlap_final = ind_overlap(get1Dindex_Collapse(ind_overlap,0))
print(ind_overlap_final)

;;;

sst_reg = regCoef_n(pmmindex_mam, month_to_season(sst, "MAM"), 0, 0)

; sst_reg = regCoef_n(pmmindex, sst_anom_rmcti, 0, 0)

copy_VarCoords(sst(1,:,:), sst_reg)

;;;
plot = new(8,graphic)
Leftstring = (/"(a) June","(b) July","(c) August","(d) September","(e) Octorber","(f) November","(g) December(+1)","(h) January(+1)"/)
; wks = gsn_open_wks("pdf", "JJA and SON SSTA tendency PMM positive year& JJA SST tendency largest")
wks = gsn_open_wks("pdf", "JJA and SON SSTA PMM positive year& JJA SST tendency largest-wind")
; wks = gsn_open_wks("pdf", "JJA and SON SSTA PMM positive year")
; wks = gsn_open_wks("pdf", "JJA and SON SSTA tendency PMM positive year")
iplot = plot
plot_vc = plot

do i =0,7
;;;SST JJA and SON tendency  next year JF
; sst_tendency_june_pmm = sst_reg
sst_tendency = sst
   if i.le.6 then 
   sst_tendency_june = sst_tendency(i+5::12,:,:)

   uwnd_june = uwnd(i+5::12,:,:)
   vwnd_june = vwnd(i+5::12,:,:)

;;;sst_tendency pmm year 
   ; sst_tendency_june_pmm = dim_avg_n_Wrap(sst_tendency_june(pmm_posi_ind,:,:),0)
   ; uwnd_pmm = dim_avg_n_Wrap(uwnd(pmm_posi_ind,:,:),0)
   ;  vwnd_pmm = dim_avg_n_Wrap(vwnd(pmm_posi_ind,:,:),0)

   sst_tendency_june_pmm = dim_avg_n_Wrap(sst_tendency_june(ind_overlap_final,:,:),0)
   uwnd_pmm = dim_avg_n_Wrap(uwnd(ind_overlap_final,:,:),0)
    vwnd_pmm = dim_avg_n_Wrap(vwnd(ind_overlap_final,:,:),0)
   

  else 
    sst_tendency_june = sst_tendency(0::12,:,:)

    uwnd_june = uwnd(0::12,:,:)
    vwnd_june = vwnd(0::12,:,:)

    ; sst_tendency_june_pmm = dim_avg_n_Wrap(sst_tendency_june(pmm_posi_ind(0:dimsizes(pmm_posi_ind)-2)+1,:,:),0)
    ; uwnd_pmm = dim_avg_n_Wrap(uwnd(pmm_posi_ind(0:dimsizes(pmm_posi_ind)-2)+1,:,:),0)
    ; vwnd_pmm = dim_avg_n_Wrap(vwnd(pmm_posi_ind(0:dimsizes(pmm_posi_ind)-2)+1,:,:),0)

    sst_tendency_june_pmm = dim_avg_n_Wrap(sst_tendency_june(ind_overlap_final+1,:,:),0)
    uwnd_pmm = dim_avg_n_Wrap(uwnd(ind_overlap_final+1,:,:),0)
    vwnd_pmm = dim_avg_n_Wrap(vwnd(ind_overlap_final+1,:,:),0)

  end if

; ;;;
   ; prob = one_test(sst_tendency_june, dim_avg_n_Wrap(sst_tendency_june, 0), pmm_posi_ind )
   ; prob_u = one_test(uwnd_june, dim_avg_n_Wrap(uwnd_june, 0), pmm_posi_ind )
   ; prob_v = one_test(vwnd_june, dim_avg_n_Wrap(vwnd_june, 0), pmm_posi_ind )

  prob = one_test(sst_tendency_june, dim_avg_n_Wrap(sst_tendency_june, 0), ind_overlap_final )
   prob_u = one_test(uwnd_june, dim_avg_n_Wrap(uwnd_june, 0), ind_overlap_final )
   prob_v = one_test(vwnd_june, dim_avg_n_Wrap(vwnd_june, 0), ind_overlap_final )



   uwnd_pmm = where(prob_u.le.0.1.or.prob_v.le.0.1, uwnd_pmm, uwnd_pmm@_FillValue)
;;;draw plot 
;;**************************************************************
;画图
;;***************************************************************
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnLeftStringFontHeightF=20.0
    res@cnFillPalette="BlueWhiteOrangeRed"
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@gsnSpreadColorStart = 130
    ; res@gsnSpreadColorEnd = 253


    res@cnFillOn=True
    res@cnLinesOn=False
    res@cnLineLabelInterval = 2
    res@gsnAddCyclic = False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
    res@lbLabelStride=2
    res@lbOrientation="Vertical"
    res@mpLimitMode           = "LatLon"
    res@mpFillDrawOrder = "PostDraw"
    res@mpMinLatF             = -40
    res@mpMaxLatF             = 40
    res@mpMinLonF             =125
    res@mpMaxLonF             = 290

    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF=0.015
    res@tmYLLabelFontHeightF=0.015
    ; res@tmXBMajorThicknessF=3.0
    ; res@tmYLMajorThicknessF=3.0
    ; res@tmBorderThicknessF=4.0
    res@gsnRightString = ""

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -0.45
    ; res@cnMaxLevelValF       =   0.45
    ; res@cnLevelSpacingF      =  0.05
    res@cnMinLevelValF       =   -2.0
    res@cnMaxLevelValF       =   2.0
    res@cnLevelSpacingF      =  0.5

    ; res@cnMinLevelValF       =   -0.4
    ; res@cnMaxLevelValF       =   0.4
    ; res@cnLevelSpacingF      =  0.1
    
    res@gsnScalarContour        = True

    ;;;;;;;LineArrow
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.02           ; thin vectors
    ; res@vcMinMagnitudeF          =0.04
    res@vcRefMagnitudeF          = 1.5
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"

    res@gsnLeftString = Leftstring(i)
    plot(i) = gsn_csm_vector_scalar_map(wks, uwnd_pmm,vwnd_pmm, sst_tendency_june_pmm, res)
   

 
  ; ;;;添加indo-pacific warm poo
  ;   ; ;;添加区域
  ;  ;;warm pool
  ;  pgres@gsLineThicknessF=2.0
  ;  pgres@gsLineDashPattern=0
  ;  str3=unique_string("dum3")
  ;  str4=unique_string("dum4")
   
  ;  ;;;画box
  ;  lat1 = 0
  ;  lat2 = 15
  ;  lon1 = 125
  ;  lon2 = 160

  ;  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  ;  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;  ;;;
  ; plot@$str2(i*3)$ = gsn_add_polyline(wks, plot(i*3), pgx_pool, pgy_pool, pgres)
  ; plot@$str2(i*3+1)$ = gsn_add_polyline(wks, plot(i*3+1), pgx_pool, pgy_pool, pgres)
  ; plot@$str2(i*3+2)$ = gsn_add_polyline(wks, plot(i*3+2), pgx_pool, pgy_pool, pgres)
  
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
    
      ; iplot(i)=gsn_csm_contour(wks,prob,res_pro)
      ; overlay(plot(i), iplot(i))

end do 

  
  pres = True
  pres@gsnPanelLabelBar  = True
  pres@tiMainString = "+PMM>1.0 SSTA tendency"

  ; pres@gsnPanelFigureStrings = (/"a","b","c"/)
  gsn_panel(wks, plot, (/4,2/), pres)
 

end 
  




