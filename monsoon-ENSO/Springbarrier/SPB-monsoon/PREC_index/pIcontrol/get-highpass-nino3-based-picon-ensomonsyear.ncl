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


;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999

; ;;;nino3
;     lat1 = -5
;     lat2 = 5
;     lon1 = 210
;     lon2 = 270
; ;;;;

;;;nino3
    lat1 = -5
    lat2 = 5
    lon1 = 150
    lon2 = 200
;;;;



; ;;;读入海温数据
    ncout = addfile("ssta_picontrol_highpass_bw.nc", "r")
    ssta = ncout->ssta 
;;


;;;;;
    nino34 = dim_avg_n_Wrap(ssta(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    nino34_12_r =  split_12(nino34)
    nino34_djf = month_to_season(nino34, "DJF")
    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)

    std_nino = dim_stddev(nino34_djf)   
    sst_total_year = ispan(startyear, endyear, 1) ;;;0-199
;;;;


;;;;
    ncout = addfile("uv_picontrol_highpass_bw.nc", "r")
    u_850 = ncout->u_850
    v_850 = ncout->v_850
;;;;

;;;
    ncout = addfile("pre_picontrol_highpass_bw.nc", "r")
    prea = ncout->pre 
    ; delete(prea)
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
   total_year = ispan(startyear,endyear-1,1)   ;;1-399
   std = dim_stddev_Wrap(DMI_ctrl_JJA)
;;;;;  


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
 nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-3)  ;;;少取一年

 DMI_ctrl_JJA_2 = DMI_ctrl_JJA(0:dimsizes(DMI_ctrl_JJA)-3)

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


    system("rm piControl_monsoon_and_ensoyear_highpass_bw_nino4.nc")
    year_out = addfile("piControl_monsoon_and_ensoyear_highpass_bw_nino4.nc", "c")

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