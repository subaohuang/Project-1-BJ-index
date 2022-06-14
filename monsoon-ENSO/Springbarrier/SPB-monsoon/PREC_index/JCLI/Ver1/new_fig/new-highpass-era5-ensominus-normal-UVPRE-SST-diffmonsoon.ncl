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

function one_test(pv_ano,pv_var,dim)
 begin
    pv_east=pv_ano
    n=dim
    var_pv_ano_east    =  pv_var
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
function composite_mean(ssta_month_12,p_month_12,u_850_month_12,v_850_month_12,ind1,ind2,ind3)
begin 
    
    if dimsizes(ind1).eq.1 then 

    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0/))
    p_month_str = dim_avg_n_Wrap(p_month_12(:,ind2,:,:),(/0/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0/))
  else
    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0,1/))
    p_month_str = dim_avg_n_Wrap(p_month_12(:,ind2,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/)) 
  end if 
    var_ssta_str = var_get(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)  ;;; get variance
    var_p_str = var_get(dim_avg_n_Wrap(p_month_12(:,:,:,:),0), ind2)
    var_u_str = var_get(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    var_v_str = var_get(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,p_month_str,u_month_str,v_month_str,var_ssta_str,var_p_str,var_u_str,var_v_str/]
end 
function Lanczos_highpass(var)
begin
    month1 = 10*12*1.0
    ; month2 = 12.0
    nwt = 49     ;;; loose 24 months each end
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
    startmon = 197901
    endmon = 202012
    startyear = 1979
    endyear = 2020

;;;读入HadiSST海温数据
    diri="/home/ys17-23/Extension/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude
    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019
    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
;;;;

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
;;;;;
    ; startmon = 187101
    ; endmon = endyear212
    ; startyear = 1871
    ; endyear = endyear2

    ; fils1 = "/home/ys17-23/Extension/20thCRV3/air.2m.mon.mean.nc"
    ; f1 = addfile(fils1, "r")
    ; time = f1->time
    ; date=cd_calendar(time, 1)

    ; endmon = endmon                 ;;; 1871-endyear2
    ; startmon = startmon 
    ; timeind=ind(date.ge.startmon.and.date.le.endmon) 
    ; sst :=f1->air(timeind,:,:)
;;;;


    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;

    ;;; interannual variance
    ; ssta = (/Lanczos_highpass(ssta)/)
    ; ssta = (/bw_bandpass(ssta)/)
    ;;;

    delete(date)
    delete(timeind)
     
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

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


;;;读入ERA wind data
    ; ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    ; fils1 = "/home/ys17-23/Extension/20thCRV2/prate.mon.mean.nc"

    ; f1 = addfile(fils1, "r")
    ; time1 = f1->time
    ; date1=cd_calendar(time1, 1)
   
    ; timeind1=ind(date1.ge.startmon.and.date1.le.194712)  

    ; fils2 = "/home/ys17-23/Extension/NCEP-NCAR/monthly/prate.sfc.mon.mean.nc"
    ; f2 = addfile(fils2, "r")
    ; time2 = f2->time
    ; date2=cd_calendar(time2, 1)
    ; timeind2=ind(date2.ge.194801.and.date2.le.endmon)  
    ; p_1 =  f1->prate(timeind1,:,:)
    ; p_2 =  f2->prate(timeind2,:,:)
     ; p = array_append_record(p_1, p_2, 0)
    

    fils2 = "/home/ys17-23/Extension/ERA5/ERA5-monthly/surface/tp_mon_1x1_1979_2020.nc"
    f2 = addfile(fils2, "r")
    time2 = f2->time
    date2=cd_calendar(time2, 1)
    timeind2=ind(date2.ge.197901.and.date2.le.endmon)  
    
    p  = short2flt(f2->tp(timeind2,:,:))



;;;;prate
    p = (/p*100.0*10/) ;; mm/day
    p = (/rmMonAnnCycTLL(p)/)

    ;;detrend 
     p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
     ;;;;filter 
     ; p = (/Lanczos_highpass(p)/)
     ; p = (/bw_bandpass(p)/)
    ;;;
    ;;;

    p_indian = dim_avg_n_Wrap(p(:,{10:30},{60:100}), (/1,2/))

    p_indian_jjas = month_to_season(p_indian, "JJA")
    p_indian_12 = split_12(p_indian)
    p_indian_jjas = (/dim_avg_n_Wrap(p_indian_12(5:8,:),0)/)
;;;



;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("../../../All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)
    endyear2 = 2016

    pre_jjas = jjas({startyear:endyear2})  ;;;1871-endyear2

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,endyear2-startyear+1,1),pre_jjasa,False,True,0)/)

    ; pre_jjasa = (/Lanczos_highpass(pre_jjasa)/)
    ; pre_jjasa = (/bw_bandpass(pre_jjasa)/)

    ;;;!!!! alert for change to pre_jjasa
    ; p_indian_jjas = (/pre_jjasa/)

;;;
 

;;;读入ERA wind data
    diri_in = "/home/ys17-23/Extension/ERA5/ERA5-monthly/pressure/"
    fils1 = diri_in + "u_mon_1x1_1979_2020.nc"

    f1 = addfile(fils1, "r")
    time1 := f1->time
    date1:=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1871-endyear2
    startmon = startmon 
    timeind1:=ind(date1.ge.startmon.and.date1.le.endmon)  
;;;;
 
;;;;U,V

    u_200  = short2flt(f1->u(timeind1,{200},:,:))
    u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)

    printVarSummary(u_850)
   
;;;


;;;
    fils2 = diri_in + "v_mon_1x1_1979_2020.nc"
    f2 = addfile(fils2, "r")
    time2 := f2->time
    date2 := cd_calendar(time2, 1)

    timeind2 := ind(date2.ge.startmon.and.date2.le.endmon)  ;;;风场数据从1979年开始

;;;;U,V

    v_200  = short2flt(f2->v(timeind2,{200},:,:))
    v_850  = short2flt(f2->v(timeind2,{850},:,:))

    v_850 = (/rmMonAnnCycTLL(v_850)/)
    v_200 = (/rmMonAnnCycTLL(v_200)/)
   
    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)

  ;;;filter 
   ; u_850 = (/Lanczos_highpass(u_850)/)
   ; v_850 = (/Lanczos_highpass(v_850)/)

   ; u_850 = (/bw_bandpass(u_850)/)
   ; v_850 = (/bw_bandpass(v_850)/)
  

;;;wind index :difference in the 850-hPa zonal winds 
;;;between a southern region (5–15N, 40–80E) and a northern region (20–30N, 70–90E)
   
  u_850_diff = dim_avg_n_Wrap(u_850(:,{5:15},{40:80}), (/1,2/)) - dim_avg_n_Wrap(u_850(:,{20:30},{70:90}), (/1,2/))

  copy_VarCoords(u_850(:,1,1), u_850_diff)

  u_indian_jjas = month_to_season(u_850_diff, "JJA")
  
;   ;;;!! alert 
;   p_indian_jjas = u_indian_jjas
; ;;;;


;;;

    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)
    dim_u = dimsizes(u_850)

;;;;;

;;;split to 12 month    3月开始，第一年是1871年
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 - 1
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

  dim2  = dimsizes(u_850)
  yearnum := dim2(0)/12 - 1 
  u_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)
  v_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)
  
  dim3  = dimsizes(p)
  yearnum := dim3(0)/12 - 1 

  p_month_12 = new((/12,yearnum,dim3(1),dim3(2)/),float)
printVarSummary(p_month_12)
printVarSummary(u_850_month_12)
printVarSummary(p)
printVarSummary(ssta)



 do i=0,11   
     if i.le.9 then
       p_month_12(i,:,:,:) = p(i+2:dimsizes(p&time)-1-12:12,:,:)
       ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1-12:12,:,:) 

       u_850_month_12(i,:,:,:) = u_850(i+2:dimsizes(u_850&time)-1-12:12,:,:)
       v_850_month_12(i,:,:,:) = v_850(i+2:dimsizes(v_850&time)-1-12:12,:,:)
     else
        p_month_12(i,:,:,:) = p(i+2:dimsizes(p&time)-1:12,:,:)
        ssta_month_12(i,:,:,:) = ssta(i+2:dimsizes(ssta&time)-1:12,:,:) 
        u_850_month_12(i,:,:,:) = u_850(i+2:dimsizes(u_850&time)-1:12,:,:)
        v_850_month_12(i,:,:,:) = v_850(i+2:dimsizes(v_850&time)-1:12,:,:)
     end if
 end do 


;;;;;
    sst_year = ispan(startyear,endyear,1) 
    u_year = ispan(startyear,endyear2,1) 
;;;;;;;

;;;mon_elnino_years 

    DMI_ctrl_JJA_r = pre_jjasa 
    total_year = ispan(startyear,endyear2,1)   ;;;1872-endyear2-1

    DMI_ctrl_JJA = DMI_ctrl_JJA_r(0:dimsizes(DMI_ctrl_JJA_r)-1)  ;;;1871-endyear2-1
    std = dim_stddev_Wrap(DMI_ctrl_JJA_r)
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-1-(endyear-endyear2))

    wet_warm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

    wet_cold_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

    wet_norm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


    dry_warm_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

    dry_cold_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

    dry_norm_year =  total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


    nor_warm_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

    nor_cold_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

    nor_norm_year =  total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))
;;;;;

;;****************************************************************************************************************************************************
  ;;;FOR 20CR based Monsoon Index
    DMI_ctrl_JJA_r = p_indian_jjas(0:dimsizes(p_indian_jjas)-1 - (endyear - endyear2) )
    DMI_ctrl_JJA = DMI_ctrl_JJA_r(0:dimsizes(DMI_ctrl_JJA_r)-1)  ;;;1871-endyear2-1
    std = dim_stddev_Wrap(DMI_ctrl_JJA_r)
;;;

    wet_warm_year_20cr = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

    wet_cold_year_20cr = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

    wet_norm_year_20cr = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


    dry_warm_year_20cr = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

    dry_cold_year_20cr = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

    dry_norm_year_20cr =  total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


    nor_warm_year_20cr = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

    nor_cold_year_20cr = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

    nor_norm_year_20cr =  total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))
;;;;;


;;;******************************************************************
;;;;;开始画图
;;;******************************************************************
  ntilde  = "n~H-13V2F35~D~FV-2H3~"
  wksstring = (/"El Ni"+ntilde+"o","La Ni"+ntilde+"a","Normal"/)
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  case1_20cr = [/wet_warm_year_20cr,dry_warm_year_20cr,nor_warm_year_20cr/]   ;;;;;

  case2_20cr = [/wet_cold_year_20cr ,dry_cold_year_20cr,nor_cold_year_20cr/]

  case3_20cr = [/wet_norm_year_20cr,dry_norm_year_20cr,nor_norm_year_20cr/]

  case_20cr = [/case1_20cr,case2_20cr,case3_20cr/]

  imsr_inten = new(6,float)

print(wet_warm_year)
print(wet_cold_year)
print(wet_norm_year)
print(dry_warm_year)
print(dry_cold_year)
print(dry_norm_year)
print(nor_warm_year)
print(nor_cold_year)
print(nor_norm_year)

; exit()
 ;;;;;
   do j=0,2

    case_test = case[j]
    ind1 :=  get1Dindex(u_year,case_test[0])
    ind2 :=  get1Dindex(u_year,case_test[1])
    ind3 :=  get1Dindex(u_year,case_test[2])

    p_elnino_str =  avg(pre_jjasa(ind1)) - avg(pre_jjasa (ind3))
    p_elnino_weak =  avg(pre_jjasa(ind2)) - avg(pre_jjasa (ind3))

    imsr_inten(j*2) = p_elnino_str
    imsr_inten(j*2+1) = p_elnino_weak
    if j.eq.1 then 
      print(ind2)
      print(ind3)
    end if
   end do 


do k=0,2

  ; case_need:= case[k]
  case_need:= case_20cr[k]

  ind_ts1 :=  get1Dindex(u_year,case_need[0])
  ind_tw1 :=  get1Dindex(u_year,case_need[1])
  ind_tnor1 :=  get1Dindex(u_year,case_need[2])

  case_need2:= case_20cr[k]
  ; case_need2:= case[k]
  ind_us1 :=  get1Dindex(u_year,case_need2[0])
  ind_uw1 :=  get1Dindex(u_year,case_need2[1])
  ind_unor1 :=  get1Dindex(u_year,case_need2[2])

  wks = gsn_open_wks("pdf", "new-era5-diff-JJASOND-UVPRE-"+" "+wksstring(k))



  month = (/"JJAS","OND"/) 

  order = (/"(a)","(b)","(c)","(d)",\
           "(e)","(f)","(g)","(h)"/)
 
  
  if k.eq. 0 then 
  stringr = (/"(EN-Dry - EN-Nor)","(EN-Wet - EN-Nor)"/)
  end if 

  if k.eq. 1 then 
  stringr = (/"(LN-Dry - LN-Nor)","(LN-Wet - LN-Nor)"/)
  end if 

  if k.eq. 2 then 
  stringr = (/"(NE-Dry)","(NE-Wet)"/)
  end if 


  month_str = ispan(1,4,1) + ""

  plot = new(dimsizes(month)*2,graphic)
  plot1 = new(dimsizes(month)*2,graphic)
  plot2 = new(dimsizes(month)*2,graphic)
  iplot = plot1
  iplot1 = plot
  iplot2 = plot


  
do i=0,dimsizes(month)-1
;;; strong and weak monsoon year p , u ,v 
   if i.eq.0 then 
    index1 = i*3+3
    index2 = i*3+3+3
    h = 0
   end if

   if i.eq.1 then 
    index1 = i*3 + 3 +1
    index2 = i*3 + 3 + 3
    h = 0
   end if


    com_month_str  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_ts1+h,ind_us1+h,ind_us1+h)
    
    ssta_month_str = com_month_str[0]
    p_month_str = com_month_str[1]
    u_month_str = com_month_str[2]
    v_month_str = com_month_str[3]
    
    ; var_ssta_str = com_month_str[4]
    ; var_p_str = com_month_str[5]
    ; var_u_str = com_month_str[6]
    ; var_v_str = com_month_str[7]

    com_month_weak  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                     v_850_month_12(index1:index2,:,:,:),ind_tw1+h,ind_uw1+h,ind_uw1+h)
    

    ssta_month_weak = com_month_weak[0]
    p_month_weak = com_month_weak[1]
    u_month_weak = com_month_weak[2]
    v_month_weak = com_month_weak[3]
    
    ; var_ssta_weak = com_month_weak[4]
    ; var_p_weak = com_month_weak[5]
    ; var_u_weak = com_month_weak[6]
    ; var_v_weak = com_month_weak[7]

    com_month_nor  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_tnor1+h,ind_unor1+h,ind_unor1+h)
    
    ssta_month_nor = com_month_nor[0]
    p_month_nor = com_month_nor[1]
    u_month_nor = com_month_nor[2]
    v_month_nor = com_month_nor[3]
    
    ; var_ssta_nor = com_month_nor[4]
    ; var_p_nor = com_month_nor[5]
    ; var_u_nor = com_month_nor[6]
    ; var_v_nor = com_month_nor[7]
    ; p_month_str = where(var_p_str.le.0.1, p_month_str, p_month_str@_FillValue)
    ; p_month_weak = where(var_p_weak.le.0.1, p_month_weak, p_month_weak@_FillValue)
;;;;;;

;;;; significant test
  ;;;; p
  ; var1 := (/var_p_str/)
  ; var2 := (/var_p_nor/)
  ; var3 := (/var_p_weak/)

  ; s1 = var_p_str@dim
  ; s2 = var_p_nor@dim
  ; s3 = var_p_weak@dim

  ; printVarSummary(var1)
  printVarSummary(p_month_str)


  ; if k.ne.2 then 

  ; prob_p_str = (/ttest(p_month_str, var1, s1, p_month_nor, var2, s2, False, False)/)
  ; prob_p_weak = (/ttest(p_month_weak, var3, s3, p_month_nor, var2, s2, False, False)/)
  ; copy_VarCoords(p_month_str,prob_p_str)
  ; copy_VarCoords(p_month_str,prob_p_weak)
  ; else 
  ; prob_p_str = one_test(p_month_str,var1,s1)
  ; prob_p_weak = one_test(p_month_weak,var3,s3)
  ; printVarSummary(prob_p_str)
 
  ; end if 

  
  ; var1_u := (/var_u_str/)
  ; var2_u := (/var_u_nor/)
  ; var3_u := (/var_u_weak/)

  ; s1_u = var_u_str@dim
  ; s2_u = var_u_nor@dim
  ; s3_u = var_u_weak@dim


  ; if k.ne.2 then 
  ; prob_u_str = (/ttest(u_month_str, var1_u,s1_u , u_month_nor, var2_u, s2_u, False, False)/)
  ; prob_u_weak = (/ttest(u_month_weak, var3_u, s3_u, u_month_nor, var2_u, s2_u, False, False)/)
  ; copy_VarCoords(u_month_str,prob_u_str)
  ; copy_VarCoords(u_month_str,prob_u_weak)
  ; else 
  ; prob_u_str = one_test(u_month_str,var1_u,s1_u)
  ; prob_u_weak = one_test(u_month_weak,var3_u,s3_u)
  ; end if 



  ; var1_v := (/var_v_str/)
  ; var2_v := (/var_v_nor/)
  ; var3_v := (/var_v_weak/)

  ; s1_v = var_v_str@dim
  ; s2_v = var_v_nor@dim
  ; s3_v = var_v_weak@dim
  
  ; if k.ne.2 then
  ; prob_v_str = (/ttest(v_month_str, var1_v,s1_v , v_month_nor, var2_v, s2_v, False, False)/)
  ; prob_v_weak = (/ttest(v_month_weak, var3_v, s3_v, v_month_nor, var2_v, s2_v, False, False)/)

  ; copy_VarCoords(v_month_str,prob_v_str)
  ; copy_VarCoords(v_month_str,prob_v_weak)
  ; else 
  ; prob_v_str = one_test(v_month_str,var1_v,s1_v)
  ; prob_v_weak = one_test(v_month_weak,var3_v,s3_v)
  ; end if 

  ; ;;;

  ; ;;; ssta 
  ; var1 := (/var_ssta_str/)
  ; var2 := (/var_ssta_nor/)
  ; var3 := (/var_ssta_weak/)

  ; s1 := var_ssta_str@dim
  ; s2 := var_ssta_nor@dim
  ; s3 := var_ssta_weak@dim
  
  ; if k.ne.2 then
  ; prob_ssta_str = (/ttest(ssta_month_str, var1,s1 , ssta_month_nor, var2, s2, False, False)/)
  ; prob_ssta_weak = (/ttest(ssta_month_weak, var3, s3, ssta_month_nor, var2, s2, False, False)/)
  ; copy_VarCoords(ssta_month_str,prob_ssta_str)
  ; copy_VarCoords(ssta_month_str,prob_ssta_weak)
  ; else
  ; prob_ssta_str = one_test(ssta_month_str,var1,s1)
  ; prob_ssta_weak = one_test(ssta_month_weak,var3,s3)
  ; end if 
 
;;;;get diff data 
;;;;;;
  if k.ne.2 then 
 
 
  u_month_str = (/u_month_str - u_month_nor/)
  v_month_str = (/v_month_str - v_month_nor/)
  p_month_str = (/p_month_str - p_month_nor/)
  ssta_month_str = (/ssta_month_str - ssta_month_nor/)

  u_month_weak = (/u_month_weak - u_month_nor/)
  v_month_weak = (/v_month_weak - v_month_nor/)
  p_month_weak = (/p_month_weak - p_month_nor/)
  ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)
  end if 
;;;;;


  ; if i.eq.1 then 
  ; prob_u_str({0:30},{90:150}) = 0.05
  ; end if
  ; u_month_str = where(prob_u_str.le.0.1.or.prob_v_str.le.0.1, u_month_str, u_month_str@_FillValue)
  ; u_month_weak = where(prob_u_weak.le.0.1.or.prob_v_weak.le.0.1, u_month_weak, u_month_weak@_FillValue)

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
  res@mpMinLatF             = -20
  res@mpMaxLatF             = 40
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
  res@cnLevels := (/-5,-4,-3,-2,-1,1,2,3,4,5/)*0.1

  res@gsnScalarContour        = False

  ;;;;;;;LineArrow
  res@vcMinDistanceF           = 0.015           ; thin vectors
  ; res@vcMinDistanceF           = 0.016          ; thin vectors
  ; res@vcMinDistanceF           = 0.016          ; thin vectors
  res@vcMinMagnitudeF          =0.1
  ; res@vcRefMagnitudeF          = 3.0
  res@vcRefMagnitudeF          = 2.0
  res@vcRefLengthF             = 0.05     ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF = 1.0
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006
  

  res@vcRefAnnoOn = False

  res@gsnLeftString = order(i*2+1)+" "+month(i)
  res@gsnLeftStringFontHeightF = 0.02
  res@gsnRightStringFontHeightF = 0.02

  res@gsnRightString = "SST" +" "+ stringr(0)

  res@tmYLLabelsOn  = False
  res@tmYRLabelsOn  = False

  res2 = res
  res2@gsnLeftString = order(i*2+4+1)+" "+month(i)
  res2@gsnRightString = "SST" +" "+stringr(1)
 

  if i.eq.1 then 
  res2@tmXBLabelsOn = True
  res2@lbLabelBarOn=True
  res2@lbOrientation="Horizontal"
  res2@pmLabelBarWidthF = 0.6
  res2@pmLabelBarHeightF =0.05
  res2@pmLabelBarOrthogonalPosF = 0.35
  res2@lbLabelStride = 1
  res2@lbLabelAutoStride = False
  res2@lbAutoManage = False
  end if 

  ;;;; for SST
  plot(i*2+1) = gsn_csm_contour_map(wks,ssta_month_weak,res)
  plot1(i*2+1) = gsn_csm_contour_map(wks,ssta_month_str,res2)


  ;;;; for PRE
  res@tmYLLabelsOn  = True
  res@gsnLeftString = order(i*2)+" "+month(i)
  res@gsnRightString = "PRE&UV850"  + " "+stringr(0)
  

  res@gsnScalarContour        = True
  res@cnFillPalette = "precip_diff_12lev"
  ; res@cnMinLevelValF       =   -2.
  ; res@cnMaxLevelValF       =   2.
  ; res@cnLevelSpacingF      =  0.4
  ; res@cnLevels := (/-2.0,-1.6,-1.2,-0.8,-0.1,0.1,0.8,1.2,1.6,2.0/)
  ; res@cnLevels := (/-2.0,-1.6,-1.2,-0.8,-0.4,0.4,0.8,1.2,1.6,2.0/)
  res@cnLevels := (/-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
  res@cnFillColors = (/2,3,4,5,6,8,9,10,11,13,14/) - 2 

  res2 = res
  res2@gsnLeftString = order(i*2+4)+" "+month(i)
  res2@gsnRightString = "PRE&UV850"+ " "+ stringr(1)

  if i.eq.1 then 
  res2@tmXBLabelsOn = True
  res2@vcRefAnnoOn = True
  res2@vcRefAnnoString1 = "2.5 m/s"
  res2@vcRefAnnoOrthogonalPosF = 0.3
  res2@vcRefAnnoFontHeightF = 0.01
  res2@vcRefAnnoPerimOn = False

  res2@lbLabelBarOn=True
  res2@lbOrientation="Horizontal"
  res2@pmLabelBarWidthF = 0.6
  res2@pmLabelBarHeightF =0.05
  res2@pmLabelBarOrthogonalPosF = 0.35
  res2@lbLabelStride = 1
  res2@lbLabelAutoStride = False
  res2@lbAutoManage = False
  end if 
  
  plot(i*2) = gsn_csm_vector_scalar_map(wks,u_month_weak, v_month_weak, p_month_weak, res)
  plot1(i*2) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, p_month_str, res2)

  delete(res@cnFillColors)
  delete(res2@cnFillColors)
   

  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.5
  pgres@gsLineDashPattern=1
  ;;;画box
  lat1 = -5
  lat2 = 5
  lon1 = 190
  lon2 = 240

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;
  ; if i.eq.1 then 
   plot@$month_str(i*2+1)$ = gsn_add_polyline(wks, plot(i*2+1), pgx_pool, pgy_pool, pgres)
   plot1@$month_str(i*2+1)$ = gsn_add_polyline(wks, plot1(i*2+1), pgx_pool, pgy_pool, pgres)
  ; end if 
  
  lat1 = 5
  lat2 = 30
  lon1 = 60
  lon2 = 100

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  
  ; if i.eq.0 then 
   plot@$month_str(i*2)$ = gsn_add_polyline(wks, plot(i*2), pgx_pool, pgy_pool, pgres)
   plot1@$month_str(i*2)$ = gsn_add_polyline(wks, plot1(i*2), pgx_pool, pgy_pool, pgres)
  ; end if 

   ; plot2@$month(i)$ = gsn_add_polyline(wks2, plot2(i), pgx_pool, pgy_pool, pgres)

   ; str1 = unique_string("dum1")  ; "unique_string" will return a unique    

   ; ;;;打点
   ;  res_pro=True
   ;  res_pro@gsnAddCyclic = False
   ;  res_pro@gsnDraw=False
   ;  res_pro@gsnFrame=False
   ;  res_pro@cnFillOn = True

   ;  res_pro@cnLevelSelectionMode = "ExplicitLevels"
   ;  res_pro@cnLevels = (/0,0.1/)

   ;  if k.eq.2.and.i.eq.1 then 
   ;  res_pro@cnLevels = (/0,0.2/)
   ;  else
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  end if


   ;  res_pro@lbLabelBarOn = False
   ;  res_pro@cnLinesOn=False
   ;  ; res_pro@cnMonoLineColor =False
   ;  ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

   ;  res_pro@cnMonoFillPattern = True
   ;  res_pro@cnFillPattern = 17
   ;  res_pro@cnFillDotSizeF = 0.001
   ;  res_pro@cnMonoFillColor = False
   ;  res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
   ;  ; res_pro@cnLevels = ispan(-6,6,1)
   ;  res_pro@cnLineLabelsOn=False
   ;  res_pro@tmYRMode = "Automatic"
   ;  res_pro@cnInfoLabelOn=False

    
   ;  iplot(i*2+1) = gsn_csm_contour(wks,prob_ssta_weak,res_pro)
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  iplot1(i*2+1) = gsn_csm_contour(wks,prob_ssta_str,res_pro)

   ;  iplot(i*2) = gsn_csm_contour(wks,prob_p_weak,res_pro)
   ;  iplot1(i*2) = gsn_csm_contour(wks,prob_p_str,res_pro)
   ;  ; iplot2(i) = gsn_csm_contour(wks2,prob_p_nor,res_pro)

   ;  overlay(plot(i*2), iplot(i*2))
   ;  overlay(plot1(i*2), iplot1(i*2))
   ;  overlay(plot(i*2+1), iplot(i*2+1))
   ;  overlay(plot1(i*2+1), iplot1(i*2+1))
    ; overlay(plot2(i), iplot2(i))
    ; delete(prob_p_str)
    ; delete(prob_p_weak)
    ; delete(prob_u_str)
    ; delete(prob_u_weak)
    ; delete(prob_ssta_str)
    ; delete(prob_ssta_weak)

    ;;;;plot for strong monsoon 
    ;;;;plot1 for weak monsoon 
       
end do


pres = True
pres@gsnPanelLabelBar  = False
pres@pmLabelBarOrthogonalPosF = -0.02
pres@gsnPanelMainString =  wksstring(k) +" events "
plot_total = array_append_record(plot, plot1, 0)
gsn_panel(wks, plot_total, (/4,2/), pres)

end do

end 