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

;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 187101
    endmon = 201012
    startyear = 1871
    endyear = 2010

    ; diri = "/home/ys17-23/Extension/SODA_region/"
    ; fils1=systemfunc("ls "+diri+"soda_TEMP_1871-2010.nc")

    ; f1=addfile(fils1, "r")
    ; date:=cd_calendar(f1->TIME, 1)
    ; timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    ; sst = f1->TEMP(timeind,0,:,:)

    ; ; sst=lonFlip(f1->sst(timeind,:,:))
    ; printVarSummary(sst)

    ; ssta=(/rmMonAnnCycTLL(sst)/)
    ; copy_VarCoords(sst, ssta)


    ; ;;detrend 
    ; ; ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ; ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&TIME),1),ssta,False,True,0)/)
    ; ;;;
    ; delete(date)
    ; delete(timeind)
    
    ; ssta!0 = "time"
   
    ; nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    ; nino34_djf = month_to_season(nino34, "DJF")

    ; nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    ; nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    ; nino34_12_r =  split_12(nino34)

    ; nino34_jja = nino34_djf
    ; nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)

    ; nino34_djf = nino34_jja

    ; std_nino = dim_stddev(nino34_djf)

;===============================================================================     


;;;读入HadiSST海温数据
    diri="/home/ys17-23/Extension/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst_had = lonFlip(f1->sst(timeind,:,:))
    

    ssta_had =(/rmMonAnnCycTLL(sst_had)/)
    copy_VarCoords(sst_had, ssta_had)
    ssta = ssta_had

    ;;detrend 
    ssta_had = (/dtrend_msg_n(ispan(1,dimsizes(ssta_had&time),1),ssta_had,False,True,0)/)
    ;;;    
    nino34 =  dim_avg_n_Wrap(ssta_had(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_jja = nino34_djf
    nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)

    nino34_djf = nino34_jja

    std_nino = dim_stddev(nino34_djf)
;;



;;;
    sst_total_year = ispan(startyear, endyear, 1)

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))
;;;


;;;read data 
    fils1 = "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/SODA/uo_his_mixed_anom.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    startmon = startmon 
    endmon = endmon                 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)    

    uo = f1->uo_his_mixed_anom(timeind1,:,:)

    fils2 = "/home/ys17-23/Extension/SODA_region/soda_ssh_1871-2010.nc"
    f2 = addfile(fils2, "r")
    ssh = f2->SSH(timeind1,:,:)

    fils3 = "/home/ys17-23/Extension/SODA_region/soda_taux_1871-2010.nc"
    f3 = addfile(fils3, "r")
    taux = f3->TAUX(timeind1,:,:)
;;;;;


;;;;anomaly
    ssh_CLM = clmMonTLL(ssh)
    ssha = calcMonAnomTLL(ssh, ssh_CLM)
    ssha = (/dtrend_msg_n(ispan(1,dimsizes(ssha&TIME),1),ssha,False,True,0)/)

    ssha = (/ssha*100.0/) ;;convert to cm
    
    uo_CLM = clmMonTLL(uo)
    uoa = calcMonAnomTLL(uo, uo_CLM)
    uoa = (/dtrend_msg_n(ispan(1,dimsizes(uoa&time),1),uoa,False,True,0)/)
    uoa!0 = "time"

    taux_CLM = clmMonTLL(taux)
    tauxa = calcMonAnomTLL(taux, taux_CLM)
    tauxa = (/dtrend_msg_n(ispan(1,dimsizes(tauxa&TIME),1),tauxa,False,True,0)/)
    tauxa!0 = "time"
;;;;


;;;AISMR DATA 
    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("../All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:2010}) 

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2010-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;

;;;strong and weak monsoon year 
    total_year = ispan(startyear,2009,1)   ;;;1872-2011
    DMI_ctrl_JJA_2 = DMI_ctrl_JJA(0:dimsizes(DMI_ctrl_JJA)-2)  ;;;1871-2011
    std = dim_stddev_Wrap(DMI_ctrl_JJA)

    inds = ind(DMI_ctrl_JJA_2.gt.0.43*std)
    indw = ind(DMI_ctrl_JJA_2.lt.-0.43*std)
    indnor = ind(DMI_ctrl_JJA_2.gt.-0.43*std.and.DMI_ctrl_JJA_2.lt.0.43*std)

    strmonrainfall = dim_avg_n_Wrap(DMI_ctrl_JJA_2(inds),0)/120.0
    print(strmonrainfall)
    

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
;;;;


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal


 nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-1-(endyear-2009))

 wet_warm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 wet_cold_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 dry_warm_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 dry_cold_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 dry_norm_year =  total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 
 nor_warm_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 nor_cold_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 nor_norm_year =  total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))
;;;;

    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)


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
       ssh_month_12(i,:,:,:) = ssha(i+2:dimsizes(ssha&TIME)-1-12:12,:,:)
       uo_month_12(i,:,:,:) = uoa(i+2:dimsizes(uoa&time)-1-12:12,:,:)
       taux_month_12(i,:,:,:) = tauxa(i+2:dimsizes(tauxa&time)-1-12:12,:,:)
     else
       ssh_month_12(i,:,:,:) = ssha(i+2:dimsizes(ssha&TIME)-1:12,:,:)
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

 wks = gsn_open_wks("pdf", "HADISST-"+stringlesft(k)+"-TAUXSSHUO-minusnormal-diffmonsoon")
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
  res@cnMinLevelValF       = -0.02               ; min level
  res@cnMaxLevelValF       =  0.02           ; max leve
  res@cnLevelSpacingF      = 0.004


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
  res@cnMinLevelValF       = -0.2               ; min level
  res@cnMaxLevelValF       =  0.2           ; max leve
  res@cnLevelSpacingF      = 0.04
  
  res@tmYLLabelsOn = False

  plot(1) =  gsn_csm_hov(wks,uo_month_weak , res)  
  
  res@gsnLeftString = "(b3) SSH' "
  res@cnMinLevelValF       = -4.0               ; min level
  res@cnMaxLevelValF       =  4.0           ; max leve
  res@cnLevelSpacingF      = 1.0

  plot(2) =  gsn_csm_hov(wks,ssh_month_weak , res)  
  

  res@gsnLeftString = "(c1) TAUX' "
  res@gsnRightString = "El Ni"+ntilde+"o (Str - Nor monsoon)"
  res@lbLabelBarOn=True
  res@pmLabelBarWidthF = 0.5
  res@pmLabelBarHeightF = 0.15
  ; res@lbAutoManage = False
  res@lbLabelAutoStride = False
  res@lbLabelStride = 2

  res@cnMinLevelValF       = -0.02               ; min level
  res@cnMaxLevelValF       =  0.02           ; max leve
  res@cnLevelSpacingF      = 0.004
  
   res@tmYLLabelsOn = True
   res@tmXBLabelsOn = True

  plot(3) =  gsn_csm_hov(wks,taux_month_str, res)  
  res@gsnLeftString = "(c2) UO' "
  res@cnMinLevelValF       = -0.2               ; min level
  res@cnMaxLevelValF       =  0.2           ; max leve
  res@cnLevelSpacingF      = 0.04

    res@tmYLLabelsOn = False

  plot(4) =  gsn_csm_hov(wks,uo_month_str, res)  
  
  res@gsnLeftString = "(c3) SSH' "
  res@cnMinLevelValF       = -4.0               ; min level
  res@cnMaxLevelValF       =  4.0           ; max leve
  res@cnLevelSpacingF      = 1.0
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