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
    ; u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,{lat1:lat2},{lon1:lon2}),(/1,2,3/))
    ; v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,{lat1:lat2},{lon1:lon2}),(/1,2,3/))

return [/ssta_month_str/]
end 


;;*************************************************
;; main program 
;;*************************************************************
begin
    ; startmon = 101
    ; endmon = 89712
    ; startyear = 1
    ; endyear = 897

    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999
;;;读入海温数据
    diri="/home/ys17-23/CMIP6/1.0x1.0/CESM2/piControl/Omon/"
    fils1=systemfunc("ls "+diri+"tos_Omon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")

    f1=addfiles(fils1, "r")
    date:=cd_calendar(f1[:]->time, 1)

    lat = f1[0]->lat
    lon = f1[0]->lon

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1[:]->tos(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
    nino34_djf = month_to_season(nino34, "DJF")

    printVarSummary(nino34_djf)

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_12_r =  split_12(nino34)

    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)


    std_nino = dim_stddev(nino34_djf)

    sst_total_year = ispan(startyear, endyear, 1) ;;;0-199

    ; nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    ; lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

    nino_year = sst_total_year(ind(nino34_djf.gt.0.43*std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.lt.-0.43*std_nino))
    normal_year = sst_total_year(ind(nino34_djf.gt.-0.43*std_nino.and.nino34_djf.lt.0.43*std_nino))

    print("************** finish SST data*****************")
;;;



;;;PRE DATA 
    diri="/home/ys17-23/CMIP6/1.0x1.0/CESM2/piControl/"
    fils2=systemfunc("ls "+diri+"pr_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")

    f2=addfiles(fils2, "r")
    date:=cd_calendar(f2[:]->time, 1)
    lat = f2[0]->lat
    lon = f2[0]->lon
    
    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    pre = lonFlip(f2[:]->pr(timeind,:,:))
  
    prea = pre
    
    prea = (/rmMonAnnCycTLL(pre)/)
    delete(pre)

    prea = (/dtrend_msg_n(ispan(1,dimsizes(prea&time),1),prea,False,True,0)/)

    prea_indian = dim_avg_n_Wrap(prea(:,{10:30},{70:90}), (/1,2/))
    delete(prea)

    prea_indian_12 = split_12(prea_indian)

    pre_jjasa = month_to_season(prea_indian, "MAM")
    delete(prea_indian)

    pre_jjasa = (/dim_avg_n_Wrap(prea_indian_12(5:8,:), 0)/)
    print("************** finish PRE data*****************")
;;;

   
;;;
    DMI_ctrl_JJA_r = pre_jjasa
    DMI_ctrl_JJA = DMI_ctrl_JJA_r(0:dimsizes(DMI_ctrl_JJA_r)-2)  
    std = dim_stddev_Wrap(DMI_ctrl_JJA_r)
;;;

;;;strong and weak monsoon year 
   total_year = ispan(startyear,endyear-1,1)   ;;1-399

   inds = ind(DMI_ctrl_JJA.ge.0.43*std)
   indw = ind(DMI_ctrl_JJA.le.-0.43*std)

   ; print(DMI_ctrl_clm)

   indnor = ind(DMI_ctrl_JJA.gt.-0.43*std.and.DMI_ctrl_JJA.lt.0.43*std)
 
   monstryear = total_year(inds)
   monwekyear = total_year(indw)
   noryear = total_year(indnor)

   ; print(monstryear)
   ; print(monwekyear)
   ; print(noryear)

   
;;;;;  


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal


 nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-2)  ;;;少取一年

 wet_warm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 wet_cold_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 dry_warm_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 dry_cold_year = total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 dry_norm_year =  total_year(ind(DMI_ctrl_JJA.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 
 nor_warm_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino))

 nor_cold_year = total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino))
 
 nor_norm_year =  total_year(ind(DMI_ctrl_JJA.ge.-0.43*std.and.DMI_ctrl_JJA.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 print(wet_warm_year)
 print(wet_cold_year)
 print(wet_norm_year)
 
 print(dry_warm_year)
 print(dry_cold_year)
 print(dry_norm_year)

 print(nor_warm_year)
 print(nor_cold_year)
 print(nor_norm_year)

   system("rm piControl_monsoon_and_ensoyear.nc")
    year_out = addfile("piControl_monsoon_and_ensoyear.nc", "c")

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

   ; exit()
end 