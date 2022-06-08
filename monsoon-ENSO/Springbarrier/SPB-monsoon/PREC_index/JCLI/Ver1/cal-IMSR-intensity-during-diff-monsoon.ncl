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


;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 187101
    endmon = 201912
    startyear = 1871
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
   
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
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



endyear2 = 2012

;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:endyear2}) 

    pre_jjas = pre_jjas*0.1/122.0 ;;convert to mm/day 
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,endyear2-startyear+1,1),pre_jjasa,False,True,0)/)
;;;
    DMI_ctrl_JJA = pre_jjasa
;;;

;;;strong and weak monsoon year 
    total_year = ispan(startyear,endyear2,1)
    std = dim_stddev_Wrap(DMI_ctrl_JJA)

    inds = ind(DMI_ctrl_JJA.gt.0.43*std)
    indw = ind(DMI_ctrl_JJA.lt.-0.43*std)
    indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
;;;;


corr = escorc_n(nino34_djf(0:dimsizes(nino34_djf)-1-(endyear- endyear2)), pre_jjasa,0, 0)


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal


 nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-1-(endyear- endyear2))

 wet_warm_year = ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino)

 wet_cold_year = ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino)

 wet_norm_year = ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino)


 dry_warm_year = ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino)

 dry_cold_year =ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino)
 
 dry_norm_year =  ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino)

 
 nor_warm_year = ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino)

 nor_cold_year = ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino)
 
 nor_norm_year =  ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino)
;;;;

  yearnum = dimsizes(nino34)/12 -1
  dim_sst = dimsizes(ssta)


ismr_wet_warm = dim_avg_n_Wrap(DMI_ctrl_JJA(wet_warm_year),0)
ismr_wet_cold = dim_avg_n_Wrap(DMI_ctrl_JJA(wet_cold_year),0)
ismr_wet_norm = dim_avg_n_Wrap(DMI_ctrl_JJA(wet_norm_year),0)

ismr_dry_warm = dim_avg_n_Wrap(DMI_ctrl_JJA(dry_warm_year),0)
ismr_dry_cold = dim_avg_n_Wrap(DMI_ctrl_JJA(dry_cold_year),0)
ismr_dry_norm = dim_avg_n_Wrap(DMI_ctrl_JJA(dry_norm_year),0)


ismr_nor_warm = dim_avg_n_Wrap(DMI_ctrl_JJA(nor_warm_year),0)
ismr_nor_cold = dim_avg_n_Wrap(DMI_ctrl_JJA(nor_cold_year),0)
ismr_nor_norm = dim_avg_n_Wrap(DMI_ctrl_JJA(nor_norm_year),0)


imsr_inten = (/ismr_wet_warm,ismr_dry_warm ,ismr_wet_cold,ismr_dry_cold,ismr_wet_norm, ismr_dry_norm /)


opt="OBS-ISMR_inten_diff_mons.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",imsr_inten))  

end 
