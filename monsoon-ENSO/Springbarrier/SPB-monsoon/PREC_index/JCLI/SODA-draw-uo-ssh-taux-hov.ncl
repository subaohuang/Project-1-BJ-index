;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
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
function composite_mean(ssta_month_12,uo_month_12,taux_month_12 ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
   ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,{lat1:lat2},:),(/1,2/))

   uo_month_str = dim_avg_n_Wrap(uo_month_12(:,ind1,{lat1:lat2},:),(/1,2/))

   taux_month_str = dim_avg_n_Wrap(taux_month_12(:,ind1,{lat1:lat2},:),(/1,2/))

return [/ssta_month_str/]
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

;;===============================================================================     


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
    sst_total_year = ispan(startyear, endyear, 1)
;;;




;;;
    sst_total_year = ispan(startyear, endyear, 1)

    ; nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    ; lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))
;;;


;;;read data 
    fils1 = "/home/ys17-23/Extension/SODA_region/UVEL_mixed_soda.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->TIME
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
    uoa = (/dtrend_msg_n(ispan(1,dimsizes(uoa&TIME),1),uoa,False,True,0)/)
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

;;;split to 12 month 
  dim1  = dimsizes(uo)
  yearnum = dim1(0)/12 -1
  uo_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)


  uo_str_evo = new((/18,dim1(2)/),float)
  uo_weak_evo = new((/18,dim1(2)/),float)
  uo_nor_evo = new((/18,dim1(2)/),float)
  
  copy_VarAtts(uo(0:17,1,:), uo_str_evo)
  copy_VarAtts(uo(0:17,1,:), uo_weak_evo)
  copy_VarAtts(uo(0:17,1,:), uo_nor_evo)

 do i=0,11   
     if i.le.9 then
       uo_month_12(i,:,:,:) = uo(i+2:dimsizes(uo&time)-1-12:12,:,:)   ;;;3月开始，到次年2月
     else
        uo_month_12(i,:,:,:) = uo(i+2:dimsizes(uo&time)-1:12,:,:)
     end if
 end do 
;;======================================================
  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  stringlesft = (/"El Nino"," La Nina"," Normal"/)
  plot = new(3,graphic)
  midplot = plot

  stss = ispan(1,3,1)+""
do k=0,2 
  print("***********start*********************"+k)

   case_need:= case[k]

   ind_ts1 :=  get1Dindex(sst_year,case_need[0])
   ind_tw1_or :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])

   ind_tnor1_r := ind_tnor1(0:dimsizes(ind_tnor1)-2)

   ind_tw1 :=  ind_tw1_or(0:dimsizes(ind_tw1_or)-2)  ;;去掉最后一年


  com_month_str  = composite_mean(uo_month_12,ind_ts1)
  uo_str_evo(0:11,:) = com_month_str[0]
 
  com_month_weak  = composite_mean(uo_month_12,ind_tw1)
  uo_weak_evo(0:11,:) = com_month_weak[0]

  com_month_nor  = composite_mean(uo_month_12,ind_tnor1_r)
  uo_nor_evo(0:11,:) = com_month_nor[0]

  com_month_str  = composite_mean(uo_month_12(0:5,:,:,:),ind_ts1+1)
  uo_str_evo(12:,:) = com_month_str[0]

  com_month_weak  = composite_mean(uo_month_12(0:5,:,:,:),ind_tw1+1)
  uo_weak_evo(12:,:) = com_month_weak[0]

  com_month_nor  = composite_mean(uo_month_12(0:5,:,:,:),ind_tnor1_r+1)
  uo_nor_evo(12:,:) = com_month_nor[0]

    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
 printMinMax(uo_str_evo, 1)


;;**************************************************************
;画图
;;***************************************************************
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False
  
  res@lbLabelBarOn=False
  

  res@cnFillOn=True
  res@cnFillDrawOrder="PreDraw"
  res@cnLinesOn=False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn=False

  res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
  res@cnMinLevelValF       = -0.02               ; min level
  res@cnMaxLevelValF       =  0.02           ; max leve
  res@cnLevelSpacingF      = 0.004


  res@gsnLeftString = "(g) TAUX' for El Nino (Str monsoon)"
  res@gsnRightString = ""

  ; res@vpWidthF  = 0.6
  ; res@vpHeightF = 0.4

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues= ispan(140, 280, 20)
  res@tmXBLabels = ispan(140, 280, 20) +"E"

  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025

  res@gsnLeftStringFontHeightF   = 0.025
  ; res@gsnRightStringFontHeightF  = 0.04
  
  res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  ; res@tmYLLabelsOn=True
  res@tmYUseLeft=True

  res@tmYLValues    = ispan(0,11,1)
  res@tmYLLabels = (/"M","A","M","J","J","A","S","O","N","D","J(+1)","F(+1)"/)
  
  ; res@tiYAxisString = "Month" 
  ; res@tiYAxisFontHeightF=0.02

  wks = gsn_open_wks("pdf", "HADISST-"+stringlesft(k)+"-TAUX-anomalies-evolution-during-diffmonsoon")
  ; wks = gsn_open_wks("pdf", stringlesft(k)+"-TAUX-anomalies-evolution-during-diffmonsoon")
  plot(0) =  gsn_csm_hov(wks,uo_str_evo(0:11,:) , res)  
    res@gsnLeftString = "(h) TAUX' for El Nino (Weak monsoon)"
  plot(1) =  gsn_csm_hov(wks,uo_weak_evo(0:11,:) , res)  
    res@gsnLeftString = "(i) TAUX' for El Nino (Normal monsoon)"
  plot(2) =  gsn_csm_hov(wks,uo_nor_evo(0:11,:) , res)  
  pres = True
  pres@gsnPanelLabelBar  = True
  pres@lbLabelStride=1
  gsn_panel(wks, plot, (/1,3/), pres)

end do

end 