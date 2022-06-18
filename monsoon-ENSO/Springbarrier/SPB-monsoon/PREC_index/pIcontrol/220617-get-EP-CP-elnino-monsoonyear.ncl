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

;;;nino3
    lat1 = -5
    lat2 = 5
    lon1 = 190
    lon2 = 240
;;;;

; ;;;nino4
;     lat1 = -5
;     lat2 = 5
;     lon1 = 150
;     lon2 = 200
; ;;;;

; ;;;读入海温数据
    ncout = addfile("ssta_picontrol_highpass_bw.nc", "r")
    ssta = ncout->ssta 
;;


;;;;;
    nino34 = dim_avg_n_Wrap(ssta(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
    nino34_12_r =  split_12(nino34)
    nino34_djf = month_to_season(nino34, "DJF")
    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)  ;;;; JJAS 

    std_nino = dim_stddev(nino34_djf)   
    sst_total_year = ispan(startyear, endyear, 1) ;;;0-199
;;;;

;;;;;
    nino3 = dim_avg_n_Wrap(ssta(:,{-5:5},{210:270}), (/1,2/))
    nino3_12_r =  split_12(nino3)
    nino3_djf = month_to_season(nino3, "DJF")
    nino3_djf = (/dim_avg_n_Wrap(nino3_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

;;;;;
    nino4 = dim_avg_n_Wrap(ssta(:,{-5:5},{160:210}), (/1,2/))
    nino4_12_r =  split_12(nino4)
    nino4_djf = month_to_season(nino4, "DJF")
    nino4_djf = (/dim_avg_n_Wrap(nino4_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;



; ;;;;
;     ncout = addfile("uv_picontrol_highpass_bw.nc", "r")
;     u_850 = ncout->u_850
;     v_850 = ncout->v_850
; ;;;;

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

 nino3_djf_2 = nino3_djf(0:dimsizes(nino3_djf)-3)  ;;;少取一年
 nino4_djf_2 = nino4_djf(0:dimsizes(nino4_djf)-3)  ;;;少取一年
 DMI_ctrl_JJA_2 = DMI_ctrl_JJA(0:dimsizes(DMI_ctrl_JJA)-3)
;;;;;


;;;;; IF EP (NINO3>NINO4) Considered as EP Elnino, versus CP El Nino events 
    en_year = nino34_djf_2
    la_year = en_year
    en_year = where(nino3_djf_2.ge.nino4_djf_2, 1,0)
    la_year = where(nino4_djf_2.ge.nino3_djf_2, 1,0)


    nino_ep_year = total_year(ind(nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.1))
    nino_cp_year = total_year(ind(nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.0))
    lanina_ep_year = total_year(ind(nino34_djf_2.le.-0.43*std_nino.and.la_year.eq.1))
    lanina_cp_year = total_year(ind(nino34_djf_2.le.-0.43*std_nino.and.la_year.eq.0))

    normal_year = total_year(ind(nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))
;;;;


;;;;
 wet_warm_ep_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.1))
 wet_cold_ep_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino.and.la_year.eq.1))
 
 wet_warm_cp_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.0))
 wet_cold_cp_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino.and.la_year.eq.0))

 wet_norm_year = total_year(ind(DMI_ctrl_JJA_2.gt.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))



 dry_warm_ep_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.1))
 dry_cold_ep_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino.and.la_year.eq.1))
 
 dry_warm_cp_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.0))
 dry_cold_cp_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.lt.-0.43*std_nino.and.la_year.eq.0))

 dry_norm_year = total_year(ind(DMI_ctrl_JJA_2.lt.-0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))

 nor_warm_ep_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.1))
 nor_warm_cp_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.gt.0.43*std_nino.and.en_year.eq.0))


 nor_cold_ep_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino.and.la_year.eq.1))
 nor_cold_cp_year = total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.lt.-0.43*std_nino.and.la_year.eq.0))


 nor_norm_year =  total_year(ind(DMI_ctrl_JJA_2.ge.-0.43*std.and.DMI_ctrl_JJA_2.le.0.43*std.and.nino34_djf_2.ge.-0.43*std_nino.and.nino34_djf_2.le.0.43*std_nino))


 print(wet_warm_ep_year)
 print(wet_cold_ep_year)
 
 print(wet_warm_cp_year)
 print(wet_cold_cp_year)

 print(dry_warm_ep_year)
 print(dry_cold_ep_year)
 
 print(dry_warm_cp_year)
 print(dry_cold_cp_year)

 
 print(nor_warm_ep_year)
 print(nor_cold_ep_year)
 print(nor_warm_cp_year)
 print(nor_cold_cp_year)


    system("rm piControl_ep-cp_monsoon_and_ensoyear_highpass_bw.nc")
    year_out = addfile("piControl_ep-cp_monsoon_and_ensoyear_highpass_bw.nc", "c")

    year_out->nino_ep_year = nino_ep_year
    year_out->nino_cp_year = nino_cp_year
    
    year_out->lanina_ep_year = lanina_ep_year
    year_out->lanina_cp_year = lanina_cp_year

    year_out->normal_year = normal_year


    year_out->wet_warm_ep_year = wet_warm_ep_year
    year_out->wet_cold_ep_year = wet_cold_ep_year
    year_out->wet_warm_cp_year = wet_warm_cp_year
    year_out->wet_cold_cp_year = wet_cold_cp_year

    year_out->wet_norm_year = wet_norm_year

    year_out->dry_warm_ep_year = dry_warm_ep_year
    year_out->dry_cold_ep_year = dry_cold_ep_year
    year_out->dry_warm_cp_year = dry_warm_cp_year
    year_out->dry_cold_cp_year = dry_cold_cp_year

    year_out->dry_norm_year = dry_norm_year

    year_out->nor_warm_ep_year = nor_warm_ep_year
    year_out->nor_warm_cp_year = nor_warm_cp_year

    year_out->nor_cold_ep_year = nor_cold_ep_year
    year_out->nor_cold_cp_year = nor_cold_cp_year

    year_out->nor_norm_year = nor_norm_year

print("end manuscript")
end 