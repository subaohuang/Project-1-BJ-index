;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function var_get(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(:,east_time,:),1)
    n=dimsizes(east_time)
    var = dim_variance_n(pv_ano(:,east_time,:),1)
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
function composite_mean(ssta_month_12,uo_month_12,taux_month_12,ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
  ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,{lat1:lat2},:),(/1,2/))
  uo_month_str = dim_avg_n_Wrap(uo_month_12(:,ind1,{lat1:lat2},:),(/1,2/))
  taux_month_str = dim_avg_n_Wrap(taux_month_12(:,ind1,{lat1:lat2},:),(/1,2/))

  var_ssta_str = var_get(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),2), ind1)  ;;; get variance
  var_uo_str = var_get(dim_avg_n_Wrap(uo_month_12(:,:,:,:),2), ind1)
  var_taux_str = var_get(dim_avg_n_Wrap(taux_month_12(:,:,:,:),2), ind1)
return [/ssta_month_str,uo_month_str,taux_month_str,var_ssta_str,var_uo_str,var_taux_str/]
end 


function get_ttest(ssh_month_str,ssh_month_nor,ssh_month_weak, var_ssh_str,var_ssh_nor,var_ssh_weak)
begin
  var1 := (/var_ssh_str/)
  var2 := (/var_ssh_nor/)
  var3 := (/var_ssh_weak/)

  s1 = var_ssh_str@dim
  s2 = var_ssh_nor@dim
  s3 = var_ssh_weak@dim

  printVarSummary(var1)

  prob_ssh_str = (/ttest(ssh_month_str, var1, s1, ssh_month_nor, var2, s2, False, False)/)
  prob_ssh_weak = (/ttest(ssh_month_weak, var3, s3, ssh_month_nor, var2, s2, False, False)/)

  copy_VarCoords(ssh_month_str,prob_ssh_str)
  copy_VarCoords(ssh_month_str,prob_ssh_weak)
return [/prob_ssh_str,prob_ssh_weak/]
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
    ; endmon = 201912
    endmon = 99912
    startyear = 1
    endyear = 999
;;;;
   diri_in = "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    year_out = addfile(diri_in+"piControl_monsoon_and_ensoyear_highpass_bw.nc", "r")
    ; year_out = addfile(diri_in+"piControl_monsoon_and_ensoyear_highpass.nc", "r")

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
;;;;;


;;;read data 
    fils1 = diri_in + "hb/uo_his_mixed_anom.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)              
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)    

    uo = f1->uo_his_mixed_anom(timeind1,:,:)

    fils2 = diri_in + "hb/tauuo_picontrol.nc"
    f2 = addfile(fils2, "r")
    taux = f2->tauuo_ano(timeind1,:,:)

    fils3 = diri_in + "hb/zos_picontrol.nc.nc"
    f3 = addfile(fils3, "r")
    ssh = f3->zos_ano(timeind1,:,:)
;;;;;


;;;;anomaly
    ssh_CLM = clmMonTLL(ssh)
    ssha = calcMonAnomTLL(ssh, ssh_CLM)
    ssha = (/dtrend_msg_n(ispan(1,dimsizes(ssha&time),1),ssha,False,True,0)/)

    ssha = (/ssha*100.0/) ;;convert to cm
    
    ; uo_CLM = clmMonTLL(uo)
    ; uoa = calcMonAnomTLL(uo, uo_CLM)
    uoa = uo
    uoa = (/dtrend_msg_n(ispan(1,dimsizes(uoa&time),1),uoa,False,True,0)/)
    uoa!0 = "time"

    taux_CLM = clmMonTLL(taux)
    tauxa = calcMonAnomTLL(taux, taux_CLM)
    tauxa = (/dtrend_msg_n(ispan(1,dimsizes(tauxa&time),1),tauxa,False,True,0)/)
    tauxa!0 = "time"
;;;;


;;;;
   ssha = (/bw_bandpass(ssha)/)
   tauxa = (/bw_bandpass(tauxa)/)
   uoa = (/bw_bandpass(uoa)/)




;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;

;;-=================================================

;;;split to 12 month    3月开始，第一年是1871年
  dim1  = dimsizes(ssha)
  yearnum = dim1(0)/12 - 1
  ssh_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

  dim2  = dimsizes(uoa)
  yearnum := dim2(0)/12 - 1 
  uo_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)

  taux_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)
  
printVarSummary(ssha)
printVarSummary(uoa)
printVarSummary(tauxa)

 do i=0,11   
     if i.le.9 then
       ssh_month_12(i,:,:,:) = ssha(i+2:dimsizes(ssha&time)-1-12:12,:,:)
       uo_month_12(i,:,:,:) = uoa(i+2:dimsizes(uoa&time)-1-12:12,:,:)
       taux_month_12(i,:,:,:) = tauxa(i+2:dimsizes(tauxa&time)-1-12:12,:,:)
     else
       ssh_month_12(i,:,:,:) = ssha(i+2:dimsizes(ssha&time)-1:12,:,:)
       uo_month_12(i,:,:,:) = uoa(i+2:dimsizes(uoa&time)-1:12,:,:)
       taux_month_12(i,:,:,:) = tauxa(i+2:dimsizes(tauxa&time)-1:12,:,:)
     end if
 end do 
;;;;

;;======================================================
  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  stringlesft = (/"El Nino"," La Nina"," Normal"/)
  plot = new(6,graphic)
  iplot = plot
  midplot = plot

  stss = ispan(1,3,1)+""

  ntilde  = "n~H-13V2F35~D~FV-2H3~"

do k=0,2 

 wks = gsn_open_wks("pdf", "pIcontrol-"+stringlesft(k)+"-TAUXSSHUO-minusnormal-diffmonsoon")
 ; wks = gsn_open_wks("pdf", "SODASST-"+stringlesft(k)+"-TAUXSSHUO-minusnormal-diffmonsoon")
  print("***********start*********************"+k)

  case_need:= case[k]

  ind_ts1 :=  get1Dindex(sst_year,case_need[0])
  ind_tw1 :=  get1Dindex(sst_year,case_need[1])
  ind_tnor1 :=  get1Dindex(sst_year,case_need[2])

  ; ind_tnor1_r := ind_tnor1(0:dimsizes(ind_tnor1)-2)

  ; ind_tw1 :=  ind_tw1_or(0:dimsizes(ind_tw1_or)-2)  ;;去掉最后一年


  com_month_str  = composite_mean(ssh_month_12,uo_month_12, taux_month_12, ind_ts1)

  ssh_month_str = com_month_str[0]
  uo_month_str = com_month_str[1]
  taux_month_str = com_month_str[2]

  var_ssh_str = com_month_str[3]
  var_uo_str = com_month_str[4]
  var_taux_str = com_month_str[5]
  
  com_month_weak  = composite_mean(ssh_month_12,uo_month_12, taux_month_12, ind_tw1)

  ssh_month_weak = com_month_weak[0]
  uo_month_weak = com_month_weak[1]
  taux_month_weak = com_month_weak[2]

  var_ssh_weak = com_month_weak[3]
  var_uo_weak = com_month_weak[4]
  var_taux_weak = com_month_weak[5]
  

  com_month_nor  = composite_mean(ssh_month_12,uo_month_12, taux_month_12, ind_tnor1)

  ssh_month_nor = com_month_nor[0]
  uo_month_nor = com_month_nor[1]
  taux_month_nor = com_month_nor[2]

  var_ssh_nor = com_month_nor[3]
  var_uo_nor = com_month_nor[4]
  var_taux_nor = com_month_weak[5]


  ssh_month_str = (/ssh_month_str - ssh_month_nor/)
  ssh_month_weak = (/ssh_month_weak - ssh_month_nor/)

  uo_month_str = (/uo_month_str - uo_month_nor/)
  uo_month_weak = (/uo_month_weak - uo_month_nor/)

  taux_month_str = (/taux_month_str - taux_month_nor/)
  taux_month_weak = (/taux_month_weak - taux_month_nor/)
;;;;;

printMinMax(taux_month_str,1)
printMinMax(taux_month_weak,1)
; exit()
;;;; 
  
  prob_ssh = get_ttest(ssh_month_str,ssh_month_nor,ssh_month_weak,var_ssh_str,var_ssh_nor,var_ssh_weak)
  prob_ssh_str = prob_ssh[0]
  prob_ssh_weak = prob_ssh[1]

  prob_uo = get_ttest(uo_month_str,uo_month_nor,uo_month_weak, var_uo_str,var_uo_nor,var_uo_weak)
  prob_uo_str = prob_uo[0]
  prob_uo_weak = prob_uo[1]


  prob_taux = get_ttest(taux_month_str,taux_month_nor,taux_month_weak, var_taux_str,var_taux_nor,var_taux_weak)
  prob_taux_str = prob_taux[0]
  prob_taux_weak = prob_taux[1]



    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)



;;**************************************************************
;画图
;;***************************************************************
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False
  
  res@lbLabelBarOn=False

  res@cnFillPalette = "BlueWhiteOrangeRed"
  
  res@cnFillOn=True
  res@cnFillDrawOrder="PreDraw"
  res@cnLinesOn=True
  res@gsnContourNegLineDashPattern = 1
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn=False

  res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
  res@cnMinLevelValF       = -0.01               ; min level
  res@cnMaxLevelValF       =  0.01           ; max leve
  res@cnLevelSpacingF      = 0.002


  res@gsnLeftString = "(b1) TAUX' "
  res@gsnRightString = "El Ni"+ntilde+"o (Weak - Nor monsoon)"

  ; res@vpWidthF  = 0.6
  ; res@vpHeightF = 0.4

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues= fspan(135.25, 270.25, 4)
  res@tmXBLabels = (/"135E","180","135W","90W"/)
  res@tmXBMinorValues = fspan(135.25, 270.25, 10)

  res@tmXBLabelFontHeightF = 0.023
  res@tmYLLabelFontHeightF = 0.023

  res@gsnLeftStringFontHeightF   = 0.023
  res@gsnRightStringFontHeightF  = 0.023
  
  res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  ; res@tmYLLabelsOn=True
  res@tmYUseLeft=True

  res@tmYLValues    = ispan(0,11,1)
  res@tmYLLabels = (/"M","A","M","J","J","A","S","O","N","D","J(+1)","F(+1)"/)
  
  ; res@tiYAxisString = "Month" 
  ; res@tiYAxisFontHeightF=0.02
  res@tmXBLabelsOn = False
 
  plot(0) =  gsn_csm_hov(wks,taux_month_weak, res)  
  res@gsnLeftString = "(b2) UO' "
  res@cnMinLevelValF       = -0.1               ; min level
  res@cnMaxLevelValF       =  0.1           ; max leve
  res@cnLevelSpacingF      = 0.02
  
  res@tmYLLabelsOn = False

  plot(1) =  gsn_csm_hov(wks,uo_month_weak , res)  
  
  res@gsnLeftString = "(b3) SSH' "
  res@cnMinLevelValF       = -2.0               ; min level
  res@cnMaxLevelValF       =  2.0           ; max leve
  res@cnLevelSpacingF      = .5

  plot(2) =  gsn_csm_hov(wks,ssh_month_weak , res)  
  

  res@gsnLeftString = "(c1) TAUX' "
  res@gsnRightString = "El Ni"+ntilde+"o (Str - Nor monsoon)"
  res@lbLabelBarOn=True
  res@pmLabelBarWidthF = 0.5
  res@pmLabelBarHeightF = 0.15
  ; res@lbAutoManage = False
  res@lbLabelAutoStride = False
  res@lbLabelStride = 2

  res@cnMinLevelValF       = -0.01               ; min level
  res@cnMaxLevelValF       =  0.01           ; max leve
  res@cnLevelSpacingF      = 0.002
  
   res@tmYLLabelsOn = True
   res@tmXBLabelsOn = True

  plot(3) =  gsn_csm_hov(wks,taux_month_str, res)  
  res@gsnLeftString = "(c2) UO' "
  res@cnMinLevelValF       = -0.1               ; min level
  res@cnMaxLevelValF       =  0.1           ; max leve
  res@cnLevelSpacingF      = 0.02

    res@tmYLLabelsOn = False

  plot(4) =  gsn_csm_hov(wks,uo_month_str, res)  
  
  res@gsnLeftString = "(c3) SSH' "
  res@cnMinLevelValF       = -2.0               ; min level
  res@cnMaxLevelValF       =  2.0           ; max leve
  res@cnLevelSpacingF      = .5
  plot(5) =  gsn_csm_hov(wks,ssh_month_str , res)  
 

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

    
    iplot(0) = gsn_csm_contour(wks,prob_taux_weak,res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_uo_weak,res_pro)
    iplot(2) = gsn_csm_contour(wks,prob_ssh_weak,res_pro)

    iplot(3) = gsn_csm_contour(wks,prob_taux_str,res_pro)
    iplot(4) = gsn_csm_contour(wks,prob_uo_str,res_pro)
    iplot(5) = gsn_csm_contour(wks,prob_ssh_str,res_pro)
    

    do i=0,5
    overlay(plot(i), iplot(i))
    end do
   
    ; overlay(plot2(i), iplot2(i))
    delete(prob_taux_str)
    delete(prob_uo_str)
    delete(prob_ssh_str)

    delete(prob_taux_weak)
    delete(prob_uo_weak)
    delete(prob_ssh_weak)
  

  pres = True
  pres@gsnPanelLabelBar  = False
  pres@lbLabelStride=1
  gsn_panel(wks, plot, (/3,3/), pres)

end do

end 