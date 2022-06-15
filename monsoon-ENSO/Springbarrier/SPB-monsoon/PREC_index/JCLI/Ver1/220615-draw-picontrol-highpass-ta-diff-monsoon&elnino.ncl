;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function var_get(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var = dim_variance_n(pv_ano(east_time,:,:), 0)
    var  =where(var.eq.0, var@_FillValue, var)
    ; var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    ; var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    ; var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    ; t_pv_ano_east = (pv_east)/var_pv_ano_east
    ; confi_pv_east=student_t(t_pv_ano_east , n-1)
    ; copy_VarCoords(pv_east, confi_pv_east)
    var@dim = n
 return var
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
function composite_mean(ssta_month_12,p_month_12,u_850_month_12,v_850_month_12,ind1,ind2,ind3)
begin 

    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0,1/))
    p_month_str = dim_avg_n_Wrap(p_month_12(:,ind1,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/))
    
    var_ssta_str = var_get(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)  ;;; get variance
    var_p_str = var_get(dim_avg_n_Wrap(p_month_12(:,:,:,:),0), ind1)
    var_u_str = var_get(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    var_v_str = var_get(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,p_month_str,u_month_str,v_month_str,var_ssta_str,var_p_str,var_u_str,var_v_str/]
end 


function bw_bandpass(var)
 begin
    ca    = 9.0*12       ; band start (longer period)
    cb    = 4.0        ; band end

    fca   = 1.0/ca      ; 'left'  frequency
    fcb   = 1.0/cb      ; 'right' frequency

    opt   = False       ; use default options (time series of filtered
                        ; values will be returned)

    var = (/bw_bandpass_filter (var,fca,fcb,opt,0)/)  
 return var
 end 
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999

    lat1 = -50
    lat2 = 50
    lon1 = -180
    lon2 = 180

; ;;;
  diri_in = "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"
;;;


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    year_out = addfile( "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"+\
                       "piControl_monsoon_and_ensoyear_highpass_bw.nc", "r")
      ; year_out = addfile( "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"+\
      ;                  "piControl_monsoon_and_ensoyear.nc", "r")

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
;;;;




;read data========================================================
    ; f1 = addfile("/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/ta_CESM2_1-1200_200hpa.nc", "r")
    ; ; f1 = addfile(diri_in+"vp_ud_vd_CESM2_1-1200.nc", "r")
    ; time = f1->time 
    ; date=cd_calendar(f1->time, 1)
    ; timeind:=ind(date.le.endmon.and.date.ge.startmon)
    ; vp = f1->ta(timeind,{20000.0},{-40:40},:)
    
    fils1 = systemfunc("ls " + "/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/ta_Amon_CESM2_piControl_r1i1p1f1_gn_"+"*"+".nc")
    f1 = addfiles(fils1, "r")
    ; f1 = addfile(diri_in+"vp_ud_vd_CESM2_1-1200.nc", "r")
    time = f1[:]->time 
    date=cd_calendar(f1[:]->time, 1)
    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    vp = dim_avg_n_Wrap(f1[:]->ta(timeind,{70000:20000},:,:),1)
    printVarSummary(vp)


    system("rm "+"/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/ta_CESM2_1-1200_200-700mean.nc")
    ncout = addfile("/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/ta_CESM2_1-1200_200-700mean.nc", "w")
  
     ncout->ta = vp
 
   
    f2 = addfile("/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/ua_CESM2_1-1200_200hpa.nc", "r")
    f3 = addfile("/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/va_CESM2_1-1200_200hpa.nc", "r")

    ud = f2->ua(timeind,:,{-40:40},:)
    vd = f3->va(timeind,:,{-40:40},:)
;;;;



;;; band pass 
     ud = (/bw_bandpass(ud)/)
     vd = (/bw_bandpass(vd)/)
     vp = (/bw_bandpass(vp)/)
;;;

print("========= finish read data ==============")


;;;;
  levs = (/20000.0/)

  wksstring = (/"El Nino","La Nina","Normal"/)
 
  u_year = ispan(startyear,endyear,1) 
;;;;

;;

  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

print("=========================Draw plot==================================")


do k=0,2

  case_need:= case[k]
  ind_us1 :=  get1Dindex(u_year,case_need[0])
  ind_uw1 :=  get1Dindex(u_year,case_need[1])
  ind_unor1 :=  get1Dindex(u_year,case_need[2])


  wks = gsn_open_wks("pdf", "airtemp-highpass-piControl-diff-9month-Str-normon-"+" "+wksstring(k))
  wks1 = gsn_open_wks("pdf", "airtemp-highpass-piControl-diff-9month-Weak-normon-"+" "+wksstring(k))  


lv = 0
;========================================================
    lev = levs(lv)
    u_850 = ud(:,{lev},:,:)
    v_850 = vd(:,{lev},:,:)
    ssta  = vp(:,:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    v_850 = (/rmMonAnnCycTLL(v_850)/)
    ssta= (/rmMonAnnCycTLL(ssta)/)

    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    yearnum = endyear - startyear +1 -1
    dim_u = dimsizes(u_850)

    ;;;split to 12 month    3月开始，第一年是1900年
    dim1  = dimsizes(ssta)
    yearnum = dim1(0)/12 - 1
    ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

    dim2  = dimsizes(u_850)
    yearnum := dim2(0)/12 - 1 
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



;;;******************************************************************
;;;;;开始画图
;;;******************************************************************



  month = (/"MAM","JJAS","OND","JF(+1)"/) 

  order = (/"(a)","(b)","(c)","(d)",\
           "(e)","(f)","(g)","(h)"/)

  month_str = ispan(1,8,1) + ""

  plot = new(dimsizes(month)*2,graphic)
  plot1 = new(dimsizes(month)*2,graphic)
  plot2 = new(dimsizes(month)*2,graphic)
  iplot = plot1
  iplot1 = plot
  iplot2 = plot




do i=0,dimsizes(month)-1
;;; strong and weak monsoon year p , u ,v 
   if i.eq.0 then 
    index1 = i*3
    index2 = i*3+2
    h = 0
   end if

   if i.eq.1 then 
    index1 = i*3
    index2 = i*3+3
    h = 0
   end if

   if i.eq.2 then 
    index1 = i*3+1
    index2 = i*3+3
    h = 0
   end if

   if i.eq.3 then 
    index1 = i*3+1
    index2 = i*3+2
    h = 0
   end if


    com_month_str  = composite_mean(ssta_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_us1+h,ind_us1+h,ind_us1+h)
    
    ssta_month_str = com_month_str[0]
    p_month_str = com_month_str[1]
    u_month_str = com_month_str[2]
    v_month_str = com_month_str[3]
    
    var_ssta_str = com_month_str[4]
    var_p_str = com_month_str[5]
    var_u_str = com_month_str[6]
    var_v_str = com_month_str[7]

    com_month_weak  = composite_mean(ssta_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                     v_850_month_12(index1:index2,:,:,:),ind_uw1+h,ind_uw1+h,ind_uw1+h)
    

    ssta_month_weak = com_month_weak[0]
    p_month_weak = com_month_weak[1]
    u_month_weak = com_month_weak[2]
    v_month_weak = com_month_weak[3]
    
    var_ssta_weak = com_month_weak[4]
    var_p_weak = com_month_weak[5]
    var_u_weak = com_month_weak[6]
    var_v_weak = com_month_weak[7]

    com_month_nor  = composite_mean(ssta_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_unor1+h,ind_unor1+h,ind_unor1+h)
    
    ssta_month_nor = com_month_nor[0]
    p_month_nor = com_month_nor[1]
    u_month_nor = com_month_nor[2]
    v_month_nor = com_month_nor[3]
    
    var_ssta_nor = com_month_nor[4]
    var_p_nor = com_month_nor[5]
    var_u_nor = com_month_nor[6]
    var_v_nor = com_month_nor[7]
    ; p_month_str = where(var_p_str.le.0.1, p_month_str, p_month_str@_FillValue)
    ; p_month_weak = where(var_p_weak.le.0.1, p_month_weak, p_month_weak@_FillValue)
;;;;;;


;;;;get data 
;;;;;;

  u_month_str = (/u_month_str - u_month_nor/)
  v_month_str = (/v_month_str - v_month_nor/)
  p_month_str = (/p_month_str - p_month_nor/)
  ssta_month_str = (/ssta_month_str - ssta_month_nor/)

  u_month_weak = (/u_month_weak - u_month_nor/)
  v_month_weak = (/v_month_weak - v_month_nor/)
  p_month_weak = (/p_month_weak - p_month_nor/)
  ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)
  


  ;;;; p
  var1 := (/var_p_str/)
  var2 := (/var_p_nor/)
  var3 := (/var_p_weak/)

  s1 = var_p_str@dim
  s2 = var_p_nor@dim
  s3 = var_p_weak@dim

  printVarSummary(var1)
  printVarSummary(p_month_str)


  prob_p_str = (/ttest(p_month_str, var1, s1, p_month_nor, var2, s2, False, False)/)
  prob_p_weak = (/ttest(p_month_weak, var3, s3, p_month_nor, var2, s2, False, False)/)

  copy_VarCoords(p_month_str,prob_p_str)
  copy_VarCoords(p_month_str,prob_p_weak)

  var1_u := (/var_u_str/)
  var2_u := (/var_u_nor/)
  var3_u := (/var_u_weak/)

  s1_u = var_u_str@dim
  s2_u = var_u_nor@dim
  s3_u = var_u_weak@dim

  prob_u_str = (/ttest(u_month_str, var1_u,s1_u , u_month_nor, var2_u, s2_u, False, False)/)
  prob_u_weak = (/ttest(u_month_weak, var3_u, s3_u, u_month_nor, var2_u, s2_u, False, False)/)

  copy_VarCoords(u_month_str,prob_u_str)
  copy_VarCoords(u_month_str,prob_u_weak)

  var1_v := (/var_v_str/)
  var2_v := (/var_v_nor/)
  var3_v := (/var_v_weak/)

  s1_v = var_v_str@dim
  s2_v = var_v_nor@dim
  s3_v = var_v_weak@dim

  prob_v_str = (/ttest(v_month_str, var1_v,s1_v , v_month_nor, var2_v, s2_v, False, False)/)
  prob_v_weak = (/ttest(v_month_weak, var3_v, s3_v, v_month_nor, var2_v, s2_v, False, False)/)

  copy_VarCoords(v_month_str,prob_v_str)
  copy_VarCoords(v_month_str,prob_v_weak)
  ;;;


  ;;; ssta 
  var1 := (/var_ssta_str/)
  var2 := (/var_ssta_nor/)
  var3 := (/var_ssta_weak/)

  s1 := var_ssta_str@dim
  s2 := var_ssta_nor@dim
  s3 := var_ssta_weak@dim

  prob_ssta_str = (/ttest(ssta_month_str, var1,s1 , ssta_month_nor, var2, s2, False, False)/)
  prob_ssta_weak = (/ttest(ssta_month_weak, var3, s3, ssta_month_nor, var2, s2, False, False)/)

  copy_VarCoords(ssta_month_str,prob_ssta_str)
  copy_VarCoords(ssta_month_str,prob_ssta_weak)
  ;;;;

  ; if i.eq.1 then 
  ; prob_u_str({0:30},{90:150}) = 0.05
  ; end if
  u_month_str = where(prob_u_str.le.0.1.or.prob_v_str.le.0.1, u_month_str, u_month_str@_FillValue)
  u_month_weak = where(prob_u_weak.le.0.1.or.prob_v_weak.le.0.1, u_month_weak, u_month_weak@_FillValue)

  ; p_month_str = where(prob_p_str.le.0.1, p_month_str, p_month_str@_FillValue)
  ; p_month_weak = where(prob_p_weak.le.0.1, p_month_weak, p_month_weak@_FillValue)



;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  ; 
  res@cnFillPalette = "BlueWhiteOrangeRed"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False

  res@lbLabelBarOn=False


  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -40
  res@mpMaxLatF             =40
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 290

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@mpFillOn = False

  res@tmXTOn=True
  res@tmYROn=True
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  
  res@tmXBLabelsOn = False

  res@tmXBLabelFontHeightF=0.02
  res@tmYLLabelFontHeightF=0.02
  res@tmXBMajorThicknessF=1.5
  res@tmYLMajorThicknessF=1.5
  res@tmXBMinorThicknessF = 1.5
  res@tmYLMinorThicknessF = 1.5
  res@tmBorderThicknessF=1.5

  res@cnLevelSelectionMode = "ExplicitLevels"
  ; res@cnLevels := (/-5,-4,-3,-2,-1,1,2,3,4,5/)*0.1

  if  k.eq.2 then 
  res@cnLevels := ispan(-10,10,1)*0.05
  else
   res@cnLevels := ispan(-10,10,1)*0.1
   end if 
  res@gsnScalarContour        = False

  ;;;;;;;LineArrow
  ; res@vcMinDistanceF           = 0.015           ; thin vectors

  if k.eq.2 then 
  res@vcMinDistanceF           = 0.014          ; thin vectors
    res@vcRefMagnitudeF          = 2.5
  res@vcRefLengthF             = 0.1    ; define length of v
  else
   res@vcMinDistanceF           = 0.016          ; thin vectors
     res@vcRefMagnitudeF          = 5.0
  res@vcRefLengthF             = 0.1    ; define length of v
  end if 
  ; res@vcMinDistanceF           = 0.016          ; thin vectors
  res@vcMinMagnitudeF          =0.5
  ; res@vcRefMagnitudeF          = 3.0

  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF = 1.0
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  if i.eq.3 then 
  res@tmXBLabelsOn = True
  res@vcRefAnnoOn = True
  res@vcRefAnnoOrthogonalPosF = 0.3
  res@vcRefAnnoFontHeightF = 0.01
  res@vcRefAnnoPerimOn = False
  res@vcRefAnnoString1 = "5 m/s"
  else 
  res@vcRefAnnoOn = False
  end if


  res@gsnLeftString = order(lv*4+i)+" "+month(i)
  res@gsnLeftStringFontHeightF = 0.02
  res@gsnRightStringFontHeightF = 0.02
  res@gsnRightString = "VP" 

  res@tmYLLabelsOn  = True

  ; if i.eq.3 then 
   ; res@lbLabelBarOn=True
   ; res@lbOrientation="Horizontal"
   ; res@pmLabelBarWidthF = 0.6
   ; res@pmLabelBarHeightF =0.05
   ; res@pmLabelBarOrthogonalPosF = 0.35
   ; res@lbLabelStride = 1
   ; res@lbLabelAutoStride = False
   ; res@lbAutoManage = False
  ; end if 
   

  res@gsnScalarContour        = True
 

  
  plot(lv*4 +i) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, ssta_month_str, res)
  plot1(lv*4+i) = gsn_csm_vector_scalar_map(wks1,u_month_weak, v_month_weak, ssta_month_weak, res)

 

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

   plot@$month_str(lv*4 +i)$ = gsn_add_polyline(wks, plot(lv*4+i) , pgx_pool, pgy_pool, pgres)
   plot1@$month_str(lv*4 +i)$ = gsn_add_polyline(wks1, plot1(lv*4+i) , pgx_pool, pgy_pool, pgres)
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

    
    iplot(lv*4+i)  = gsn_csm_contour(wks,prob_ssta_str,res_pro)
    iplot1(lv*4+i) = gsn_csm_contour(wks1,prob_ssta_weak,res_pro)

   
    ; iplot2(i) = gsn_csm_contour(wks2,prob_p_nor,res_pro)

    overlay(plot(lv*4+i), iplot(lv*4+i))
    overlay(plot1(lv*4+i), iplot1(lv*4+i))

    ; overlay(plot2(i), iplot2(i))
    delete(prob_p_str)
    delete(prob_p_weak)
    delete(prob_u_str)
    delete(prob_u_weak)
    delete(prob_ssta_str)
    delete(prob_ssta_weak)
end do

pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.02
pres@gsnPanelMainString =  wksstring(k) +" events (Str Mon - Normal Mon)"

; gsn_panel(wks, plot, (/4,2/), pres)
; pres@gsnPanelMainString =  wksstring(k) +" events (Weak Mon - Normal Mon)"
; gsn_panel(wks1, plot1, (/4,2/), pres)

gsn_panel(wks, (/plot/), (/2,2/), pres)
; gsn_panel(wks, (/plot(0),plot(4),plot(1),plot(5)/), (/2,2/), pres)
pres@gsnPanelMainString =  wksstring(k) +" events (Weak Mon - Normal Mon)"
gsn_panel(wks1, (/plot1/), (/2,2/), pres)

end do 

end  