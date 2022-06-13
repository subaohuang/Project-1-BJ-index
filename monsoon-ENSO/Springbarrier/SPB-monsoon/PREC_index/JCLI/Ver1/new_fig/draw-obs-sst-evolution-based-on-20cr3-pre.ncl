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
function composite_mean(ssta_month_12,ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,{lat1:lat2},{lon1:lon2}),(/1,2,3/))
 

return [/ssta_month_str/]
end 

function Lanczos_highpass(var)
begin
    month1 = 9*12*1.0
    ; month2 = 12.0
    nwt = 49    ;;; loose 24 months each end
    fca = 1./month1 ;;;通过的频率
    ; fcb = 1./month2
    ihp = 1   ;;;1 high pass  0:low pass         
    nsigma = 1.  ; Lanczos sigma

    wgts = filwgts_lanczos (nwt, ihp, fca, -999, nsigma)  

    var  = (/wgt_runave_n_Wrap(var, wgts,0,0)/)
 return var
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
    startmon = 189101
    endmon = 201912
    startyear = 1891
    endyear = 2019

;;;读入HadiSST海温数据
    diri="/home/ys17-23/Extension/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst=lonFlip(f1->sst(timeind,:,:))
   

; ;;;;;;

; ;;;读入ERSST海温数据
;     diri="/home/ys17-23/Extension/SST/NOAA_ERSST/"
;     fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")
;     f1=addfile(fils1, "r")
;     date:=cd_calendar(f1->time, 1)
;     lat = f1->lat
;     lon = f1->lon
;     timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019
;     sst=f1->sst(timeind,:,:)
;     printVarSummary(sst)
; ;;;;


    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)


    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;

    ;;; high pass 9years, lose each end 2years
    ; ssta  = (/Lanczos_highpass(ssta)/)
    ; ssta  = (/bw_bandpass(ssta)/)
    ;;;

    delete(date)
    delete(timeind)
     
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_jja = nino34_djf
    nino34_jja = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)  ;;;;JJA > 0.43*STD

    nino34_djf = nino34_jja

    std_nino = dim_stddev(nino34_djf)
    sst_total_year = ispan(startyear, endyear, 1)

    ; nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    ; lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))
;;;

;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    endmonyear = 2015

    pre_jjas = jjas({startyear:endmonyear}) 

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,endmonyear-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
   
   ; DMI_ctrl_JJA = (/bw_bandpass_filter(DMI_ctrl_JJA, 1/9.0, 1/2.1, False, 0)/)
     ;;; high pass 9years, lose each end 2years
    ; DMI_ctrl_JJA  = (/Lanczos_highpass_year(DMI_ctrl_JJA)/)
    ; DMI_ctrl_JJA_11avg = runave_n_Wrap(DMI_ctrl_JJA, 7, 0, 0)
    ; DMI_ctrl_JJA = (/ DMI_ctrl_JJA - DMI_ctrl_JJA_11avg/)
;;;


;;;
; ;;;读入ERA wind data
;     fils1 = "/home/ys17-23/Extension/GPCC/precip.mon.total.v2018.nc"  ;;;; 1891-2016
;     f1 = addfile(fils1, "r")
;     time1 = f1->time
;     date1=cd_calendar(time1, 1)

;     timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  
;     p  = f1->precip(timeind1,:,:)
;     p = (/p/30.0/) ;; mm/day
; ;;;;;

;读入ERA wind data
    ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/ys17-23/Extension/20thCRV2/prate.mon.mean.nc"

    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)
   
    timeind1=ind(date1.ge.startmon.and.date1.le.194712)  

    fils2 = "/home/ys17-23/Extension/NCEP-NCAR/monthly/prate.sfc.mon.mean.nc"
    f2 = addfile(fils2, "r")
    time2 = f2->time
    date2=cd_calendar(time2, 1)
    timeind2=ind(date2.ge.194801.and.date2.le.201512)  
    p_1 =  f1->prate(timeind1,:,:)
    p_2 =  f2->prate(timeind2,:,:)
    p = array_append_record(p_1, p_2, 0)  ;;;18 - 2010
    

;     fils2 = "/home/ys17-23/Extension/NCEP-NCAR/monthly/prate.sfc.mon.mean.nc"
;     f2 = addfile(fils2, "r")
;     time2 = f2->time
;     date2=cd_calendar(time2, 1)
;     timeind2=ind(date2.ge.197001.and.date2.le.endmon)  
    
;     p  =  f2->prate(timeind2,:,:)
;      p = (/p*86400/) ;; mm/day
; ;;;
; ;;;读入ERA wind data
;     ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
;     fils1 = "/home/ys17-23/Extension/20thCRV3/prate.mon.mean.nc"
;     f1 = addfile(fils1, "r")
;     time1 = f1->time
;     date1=cd_calendar(time1, 1)

;     endmon = endmon                 ;;; 1871-2015
;     startmon = startmon 
;     timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  
; ;;;;



;;;;prate
    ; p  = f1->prate(timeind1,:,:)
    p = (/p*86400/) ;; mm/day
;;;;prate


    
    p = (/rmMonAnnCycTLL(p)/)

    ;;detrend 
    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
     ;;;;filter 
     ; p = (/Lanczos_highpass(p)/)
     ; p = (/bw_bandpass(p)/)
    ;;;
    ;;;

    p_indian = dim_avg_n_Wrap(p(:,{10:30},{70:90}), (/1,2/))

    p_indian_jjas = month_to_season(p_indian, "JJA")
    p_indian_12 = split_12(p_indian)
    p_indian_jjas = (/dim_avg_n_Wrap(p_indian_12(5:8,:),0)/)

    DMI_ctrl_JJA = p_indian_jjas
;;;



;;;strong and weak monsoon year 
    total_year = ispan(startyear,endmonyear,1)
    std = dim_stddev_Wrap(DMI_ctrl_JJA)

    DMI_ctrl_JJA_2 = DMI_ctrl_JJA(0:dimsizes(DMI_ctrl_JJA)-1)

    inds = ind(DMI_ctrl_JJA.gt.0.43*std)
    indw = ind(DMI_ctrl_JJA.lt.-0.43*std)
    indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal


 nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-1-(endyear-endmonyear))



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

print(wet_warm_year)
print(wet_cold_year)
print(wet_norm_year)
print(dry_warm_year)
print(dry_cold_year)
print(dry_norm_year)
print(nor_warm_year)
print(nor_cold_year)
print(nor_norm_year)


;;;

    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)


;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]


  
  ntilde  = "n~H-13V2F35~D~FV-2H3~"

  stringlesft = (/"(a) El Ni"+ntilde+"o","(b) La Ni"+ntilde+"a","(c) Non-ENSO"/)

  plot = new(3,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "necp1-20crbase-BWhighpass-ssta-evolution-during-diffmonsoon")
   
  stss = ispan(1,3,1)+""


do k=0,2 
  print("***********start*********************"+k)

   case_need:= case[k]

   ind_ts1 :=  get1Dindex(sst_year,case_need[0])
   ind_tw1 :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1 :=  get1Dindex(sst_year,case_need[2])



;;;split to 12 month 
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 -1
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)


  ssta_str_evo = new((/18/),float)
  ssta_weak_evo = new((/18/),float)
  ssta_nor_evo = new((/18/),float)



 do i=0,11   
     if i.le.9 then
       ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1-12:12,:,:)   ;;;3月开始，到次年2月
     else
        ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1:12,:,:)
     end if
 end do 



    com_month_str  = composite_mean(ssta_month_12,ind_ts1)
    
    ssta_str_evo(0:11) = com_month_str[0]


    com_month_weak  = composite_mean(ssta_month_12,ind_tw1)
    
    ssta_weak_evo(0:11) = com_month_weak[0]


      ; ind_tnor1_r := ind_tnor1(0:dimsizes(ind_tnor1)-2)
      ; ind_unor1_r := ind_unor1(0:dimsizes(ind_unor1)-2)
    ; end if
    
    com_month_nor  = composite_mean(ssta_month_12,ind_tnor1)
   
    ssta_nor_evo(0:11) = com_month_nor[0]


    com_month_str  = composite_mean(ssta_month_12(0:5,:,:,:),ind_ts1+1)
    
    ssta_str_evo(12:) = com_month_str[0]


    com_month_weak  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tw1+1)
    
    ssta_weak_evo(12:) = com_month_weak[0]


    com_month_nor  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tnor1+1)
    
    ssta_nor_evo(12:) = com_month_nor[0]

   

    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
 

;;**************************************************************
;画图
;;***************************************************************
  num_x = 16

  res = True
  
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  ; res@tmXTOn         = False
  ; res@tmYROn = False


  res@trXMinF        = 0
  res@trXMaxF        = num_x-1



  res@gsnYRefLine = 0.0


  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025

  res@gsnLeftStringFontHeightF   = 0.025
  res@gsnRightStringFontHeightF  = 0.025

  res@tmBorderThicknessF = 2.0

  res@tmXBMode = "Explicit"
  ; res@tmXBValues = ispan(0, num_x-1, 1)
  res@tmXBValues = ispan(0, num_x-1, 3)
  res@tmXBMinorValues = ispan(0, num_x-1, 1)

  

  res@tmXBLabels = (/"Mar","Jun","Sep","Dec","Mar(1)","Jun(1)"/)
  ; res@tmXBLabels = (/"M","M","J","A","S","N","D","F(+1)","M(+1)","M(+1)",/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J(+1)","F(+1)","M(+1)","A(+1)","M(+1)"/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J1","F1","M1","A1","M1","J1","J1","A1"/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
  res@tmYLMode = "Explicit"



  if k.eq.0 then 
   ; res@trYMinF        = -0.6
   ; res@trYMaxF        = 1.6

   res@trYMinF        = -1.0
   res@trYMaxF        = 1.6

   ; res@tmYLValues := ispan(-6, 16, 2)*0.1
   ; res@tmYLLabels := ispan(-6, 16, 2)*0.1
  end if 

  if k.eq.1 then 
   ; res@trYMinF        = -1.2
   ; res@trYMaxF        = 0.4
   res@trYMinF        = -1.6
   res@trYMaxF        = 1.0

   ; res@tmYLValues := ispan(-12, 4, 2)*0.1
   ; res@tmYLLabels := ispan(-12, 4, 2)*0.1
  end if 

  if k.eq.2 then 

   res@trYMinF        = -0.6
   res@trYMaxF        = 0.6

   ; res@tmYLValues := ispan(-6, 6, 2)*0.1
   ; res@tmYLLabels := ispan(-6, 6, 2)*0.1
  end if 

  if k.eq.0 then 
  res@tiYAxisString = "SSTA [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.025
  else
   res@tiYAxisString = "" 
  end if
  ; res@tiXAxisString = "Month"
  ; res@tiXAxisFontHeightF = 0.025


  res@xyMarkLineModes = (/"MarkLines","MarkLines","MarkLines"/)
  res@xyMarkers = (/1,1,1/)
  res@xyMarkerSizes = (/0.03,0.03,0.03/)
  res@xyMarkerColors= (/"blue","red","black"/)


  res@xyMarkerThicknesses = (/4.0,4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0,2.0/)
  res@xyLineColors       = (/"blue","red","black"/)
  res@xyDashPatterns     = (/1,1,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = stringlesft(k)
  res@gsnRightString = "OBS"

  if k.ne.2 then 
   plot(k)  = gsn_csm_xy (wks, ispan(0,15,1), (/ssta_str_evo(0:15),ssta_weak_evo(0:15),ssta_nor_evo(0:15)/), res)
  else
   plot(k)  = gsn_csm_xy (wks, ispan(0,15,1), (/ssta_str_evo(0:15),ssta_weak_evo(0:15)/), res)
  end if 


  respl = True
  respl@gsLineColor = "gray"
  respl@gsLineDashPattern = 1.0

  x =  1.0*res@tmXBValues
  y = x

  if k.eq.0 then 
  y = 0.43*std_nino
  end if 

  if k.eq.1 then 
  y = -0.43*std_nino
  end if 
  
  if k.ne.2 then 
  plot@$stss(k)$ =  gsn_add_polyline(wks, plot(k),x, y, respl)
  end if



  ;;添加legend 
    gres = True
    ; if k.eq.0 then 
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ; else
    ; gres@YPosPercent = 25.    ; expressed as %, 0->100, sets position of top border of legend 
    ; end if
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors

    if k.ne.2 then 
    lineres@lgLineColors       = (/"red","black","blue"/)
    else
    lineres@lgLineColors := (/"red","blue"/)
    end if
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True
    
    if k.ne.2 then 
    textres@lgLabels := (/"Dry","Nor","Wet"/)
    else
    textres@lgLabels := (/"Dry","Wet"/)
    end if
    textres@lgLabelFontHeightF = 0.03
    ; textres@lgLabelOpacity = (/0/) 
    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   
   

end do
pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/1,3/), pres)
end 