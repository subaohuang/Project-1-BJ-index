;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
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

function Lanczos_highpass(var)
begin
    month1 = 10*12*1.0
    nwt = 49     ;;; loose 24 months each end
    fca = 1./month1 ;;;通过的频率
    ihp = 1     ;;;high pass           
    nsigma = 1.  ; Lanczos sigma

    wgts = filwgts_lanczos (nwt, ihp, fca, -999., nsigma)  

    var  = (/wgt_runave_Wrap(var, wgts, 0)/)
 return var
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

; ;;;读入海温数据
    diri="/home/ys17-23/CMIP6/1.0x1.0/CESM2/piControl/Omon/"
    fils1=systemfunc("ls "+diri+"tos_Omon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")

    f1=addfiles(fils1, "r")
    date:=cd_calendar(f1[:]->time, 1)

    lat := f1[0]->lat
    lon := f1[0]->lon

    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    sst=lonFlip(f1[:]->tos(timeind,{lat1:lat2},{lon1:lon2}))
    ; sst = (/sst*86400/)
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
   
    print("************** finish SST data*****************")

 ;;;; high pass 10years, lose each end 2years
    ; ssta  = (/Lanczos_highpass(ssta)/)
 ;;;
    ca    = 10.0*12.0        ; band start (longer period)
    cb    = 12.0        ; band end

    fca   = 1.0/ca      ; 'left'  frequency
    fcb   = 1.0/cb      ; 'right' frequency

    opt   = False       ; use default options (time series of filtered
                        ; values will be returned)

    ssta = (/bw_bandpass_filter (ssta,fca,fcb,opt,0)/)   ; 

    system("rm ssta_picontrol_highpass_bw.nc")
    ncout = addfile("ssta_picontrol_highpass_bw.nc", "c")
    ncout->ssta = ssta
    exit()
;;;;;;


;;;;;
    nino34 = dim_avg_n_Wrap(ssta(:,{-5:5},{190:240}), (/1,2/))
    nino34_12_r =  split_12(nino34)

    nino34_djf = month_to_season(nino34, "DJF")
    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)
    std_nino = dim_stddev(nino34_djf)

    sst_total_year = ispan(startyear, endyear+2, 1) ;;;0-199

; ;;;

;;读入ua
    diri2 = "/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
    fils2 = systemfunc("ls "+diri2+"ua_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")
  
    f2 := addfiles(fils2, "r")
    date:=cd_calendar(f2[:]->time, 1)

    lat := f2[0]->lat
    lon := f2[0]->lon

    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    u_850 = lonFlip(f2[:]->ua(timeind,{85000},{lat1:lat2},{lon1:lon2}))
    printVarSummary(u_850)
;;;;

;;读入va
    diri3 = "/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
    fils3 = systemfunc("ls "+diri2+"va_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")
  
    f3 := addfiles(fils3, "r")
    date:=cd_calendar(f3[:]->time, 1)

    lat := f3[0]->lat
    lon := f3[0]->lon

    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    v_850 = lonFlip(f3[:]->va(timeind,{85000},{lat1:lat2},{lon1:lon2}))
    printVarSummary(v_850)
;;;;
    u_850 = (/rmMonAnnCycTLL(u_850)/)
    v_850 = (/rmMonAnnCycTLL(v_850)/)

    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)

    u_850 = (/Lanczos_highpass(u_850)/)
    v_850 = (/Lanczos_highpass(v_850)/)
    

    system("rm uv_picontrol_highpass.nc")
    ncout = addfile("uv_picontrol_highpass.nc", "c")
    ncout->u_850 = u_850
    ncout->v_850 = v_850
;;;;


;;;PRE DATA 
    diri="/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
    fils2=systemfunc("ls "+diri+"pr_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")

    f2:=addfiles(fils2, "r")
    date:=cd_calendar(f2[:]->time, 1)
    lat := f2[0]->lat
    lon := f2[0]->lon
    
    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    pre = lonFlip(f2[:]->pr(timeind,:,:))
  
    prea = pre
    
    prea = (/rmMonAnnCycTLL(pre)/)
    delete(pre)
    prea = (/dtrend_msg_n(ispan(1,dimsizes(prea&time),1),prea,False,True,0)/)
;;;


;;;; write precipitation data 
    prea = (/Lanczos_highpass(prea)/)

    system("rm pre_picontrol_highpass.nc")
    ncout = addfile("pre_picontrol_highpass.nc", "c")
    ncout->pre = prea 
    delete(prea)
;;;

;;;PRE DATA 
    diri="/home/yangsong3/CMIP6/1.0x1.0/CESM2/piControl/"
    fils2:=systemfunc("ls "+diri+"pr_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")

    f2:=addfiles(fils2, "r")
    date:=cd_calendar(f2[:]->time, 1)
    lat := f2[0]->lat
    lon := f2[0]->lon
    
    timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    pre:= lonFlip(f2[:]->pr(timeind,:,:))
  
    prea := pre
    
    prea = (/rmMonAnnCycTLL(pre)/)
    delete(pre)
    prea = (/dtrend_msg_n(ispan(1,dimsizes(prea&time),1),prea,False,True,0)/)
;;;

;;;
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
   total_year = ispan(startyear+2,endyear,1)   ;;1-399

   std = dim_stddev_Wrap(DMI_ctrl_JJA)

  
   ; print(monstryear)
   ; print(monwekyear)
   ; print(noryear)

   
;;;;;  


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
 nino34_djf_2 = nino34_djf(2:dimsizes(nino34_djf)-3)  ;;;少取一年

 DMI_ctrl_JJA_2 = DMI_ctrl_JJA(2:dimsizes(DMI_ctrl_JJA)-3)

 nino_year = total_year(ind(nino34_djf_2.gt.0.43*std_nino))
 lanina_year = total_year(ind(nino34_djf_2.le.-0.43*std_nino))
 normal_year = total_year(ind(nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 monstryear = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std))
 monwekyear = total_year(ind(DMI_ctrl_JJA_2.le.-0.43*std))
 noryear = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std_nino.and.DMI_ctrl_JJA_2.le.0.43*std_nino))



 wet_warm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 wet_cold_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 dry_warm_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 dry_cold_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 dry_norm_year =  total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 
 nor_warm_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 nor_cold_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 nor_norm_year =  total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 print(wet_warm_year)
 print(wet_cold_year)
 print(wet_norm_year)
 
 print(dry_warm_year)
 print(dry_cold_year)
 print(dry_norm_year)

 print(nor_warm_year)
 print(nor_cold_year)
 print(nor_norm_year)


    system("rm piControl_monsoon_and_ensoyear_highpass.nc")
    year_out = addfile("piControl_monsoon_and_ensoyear_highpass.nc", "c")

    year_out->nino_year = nino_year
    year_out->lanina_year = lanina_year
    year_out->normal_year = normal_year

    year_out->monstryear = monstryear
    year_out->monwekyear = monwekyear
    year_out->noryear = noryear

    year_out->wet_warm_year = wet_warm_year
    year_out->wet_cold_year = wet_cold_year
    year_out->wet_norm_year = wet_norm_year

    year_out->dry_warm_year = dry_warm_year
    year_out->dry_cold_year = dry_cold_year
    year_out->dry_norm_year = dry_norm_year

    year_out->nor_warm_year = nor_warm_year
    year_out->nor_cold_year = nor_cold_year
    year_out->nor_norm_year = nor_norm_year

print("end manuscript")
end 