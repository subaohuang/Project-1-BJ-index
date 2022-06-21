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

undef("remove_enso_signal")
function remove_enso_signal(nino3_4,  data)
  begin  
   ;;;;;
      ; data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_de  = data
      data_rc = regCoef_n(nino3_4,data_de, 0, 0)
      copy_VarCoords(data(1,:,:), data_rc)
      data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
      m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
      nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))
  
      data_relate_enso  = nino3_4_3d*data_rc_3d + m 
      data_new = data 
      data_new = (/data - data_relate_enso/)
  return data_new
end 
undef("enso_eof")
function enso_eof(data)
  begin  
     ;;;;EOF 
    optEOF      = True
    optEOF@jopt = 0
    neval  = 3                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
    ev_sst = eofunc_n_Wrap(data,neval,optEOF,0)    ;;;得到的标准化后的场，要使得不标准化
    ev_sst&lon@units = "degrees_east"
    ev_sst&lat@units = "degrees_north"
    ev_timef = eofunc_ts_n_Wrap(data, ev_sst, False,0)

    do ne=0,neval-1
       ev_sst(ne,:,:) = ev_sst(ne,:,:)*sqrt(ev_sst@eval(ne)) ;;;使得量级能够相同
    end do

    eval=ev_sst@eval
    ts_model = ev_timef(0,:)/eval(0)
    ; ts_model = ev_timef(0,:)
    ; ts_model=ev_timef(1,:)/eval(1)
    ; tos_enino = ev_sst(0,:,:)
    ts_model_std = dim_standardize(ts_model, 0)
    tos_enino  =  regCoef_n(ts_model_std, data, 0, 0)
    tos_enino@var = ev_sst@pcvar
    copy_VarCoords(data(1,:,:), tos_enino)
    
  return [/tos_enino,ts_model_std/]
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

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)  ;;;; JJAS 

    std_nino = dim_stddev(nino34_djf)   
    sst_total_year = ispan(startyear, endyear, 1) ;;;0-199
;;;;

;;;;;
    nino3 = dim_avg_n_Wrap(ssta(:,{-5:5},{210:270}), (/1,2/))
    nino3_12_r =  split_12(nino3)
    nino3_djf = month_to_season(nino3, "DJF")

    nino3_djf(0:dimsizes(nino3_djf)-2) = nino3_djf(1:) 
    nino3_djf(dimsizes(nino3_djf)-1) = nino3(dimsizes(nino3)-1) 

    nino3_djf = (/dim_avg_n_Wrap(nino3_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

;;;;;
    nino4 = dim_avg_n_Wrap(ssta(:,{-5:5},{160:210}), (/1,2/))
    nino4_12_r =  split_12(nino4)
    nino4_djf = month_to_season(nino4, "DJF")

    nino4_djf(0:dimsizes(nino4_djf)-2) = nino4_djf(1:) 
    nino4_djf(dimsizes(nino4_djf)-1) = nino4(dimsizes(nino4)-1) 

    nino4_djf = (/dim_avg_n_Wrap(nino4_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

;;;;;
    nino1_2 = dim_avg_n_Wrap(ssta(:,{-10:0},{270:280}), (/1,2/))
    nino1_2_12_r =  split_12(nino1_2)
    nino1_2_djf = month_to_season(nino1_2, "DJF")

    nino1_2_djf(0:dimsizes(nino1_2_djf)-2) = nino1_2_djf(1:) 
    nino1_2_djf(dimsizes(nino1_2_djf)-1) = nino1_2(dimsizes(nino1_2)-1) 

    nino1_2_djf = (/dim_avg_n_Wrap(nino1_2_12_r(5:8,:), 0)/)  ;;;; JJAS 
;;;;

; ;;; Kao and Yu 2009

;    ;;;去掉Nino4去信号
    
;     ssta_rm_nino1_2 = ssta(:,{-20:20},{120:280})
;     ssta_rm_nino1_2 = (/remove_enso_signal(nino1_2,ssta(:,{-20:20},{120:280}))/)

;     ssta_rm_nino4 = ssta(:,{-20:20},{120:280})
;     ssta_rm_nino4 = (/remove_enso_signal(nino4,ssta(:,{-20:20},{120:280}))/)
;   ;;;
;    print("******start EOF*************")
;      ;;;;EOF 
;     eof_cp = enso_eof(ssta_rm_nino1_2)
;     eof_cp_pattern = eof_cp[0]
;     ; eof_cp_pattern = (/-1*eof_cp_pattern/)
;     pc_cp = eof_cp[1]
;     var1 = eof_cp_pattern@var

;     eof_ep = enso_eof(ssta_rm_nino4)
;     eof_ep_pattern = eof_ep[0]
;     pc_ep = eof_ep[1]
;     var2 = eof_ep_pattern@var

;     copy_VarCoords(nino1_2, pc_ep)
;     copy_VarCoords(nino1_2, pc_cp)
;     ;;;;
;   ;;输出PC的数据
;     opt="CP_PC1_picontrol_bw.ascii"    
;     asciiwrite(opt, sprintf ("%6.2f,",pc_cp))  

;     opt="EP_PC1_picontrol_bw.ascii"    
;     asciiwrite(opt, sprintf ("%6.2f,",pc_ep))  

;  print("******end EOF*************")
; ;;;;


  ;;输入PC的数据
    opt="CP_PC1_picontrol_bw.ascii"    
    pc_cp =  asciiread(opt, -1,"float")  

    opt="EP_PC1_picontrol_bw.ascii"    
    pc_ep =  asciiread(opt, -1, "float")  

    pc_cp = -1*pc_cp
    pc_ep = -1*pc_ep

    tos_enino_ep  =  regCoef_n(pc_ep, ssta, 0, 0)
    copy_VarCoords(ssta(1,:,:), tos_enino_ep)

    tos_enino_cp  =  regCoef_n(pc_cp, ssta, 0, 0)
    copy_VarCoords(ssta(1,:,:), tos_enino_cp)
    
    pc_ep_jjas  = nino34_djf
    pc_ep_12 =  split_12(pc_ep)
    pc_ep_jjas  =  (/dim_avg_n_Wrap(pc_ep_12(5:8,:), 0)/)  ;;;; JJAS 
    
    pc_cp_jjas  =  nino34_djf
    pc_cp_12 =  split_12(pc_cp)
    pc_cp_jjas  =  (/dim_avg_n_Wrap(pc_cp_12(5:8,:), 0)/)  ;;;; JJAS 
;;;    
   

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
   
    ; ;;;
    ; en_year = where(nino3_djf_2.ge.nino4_djf_2, 1,0)
    ; la_year = where(nino4_djf_2.ge.nino3_djf_2, 1,0)
    
    ;;for Kao and Yu 
    pc_ep_jjas_2 = pc_ep_jjas(0:dimsizes(pc_ep_jjas)-3)  ;;;少取一年
    pc_cp_jjas_2 = pc_cp_jjas(0:dimsizes(pc_cp_jjas)-3)  ;;;少取一年


    en_year = where(pc_ep_jjas_2.ge.pc_cp_jjas_2, 1,0)
    la_year = where(pc_cp_jjas_2.ge.pc_ep_jjas_2, 1,0)

    ;;;


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


    system("rm piControl_ep-cp_monsoon_and_ensoyear_highpass_bw_eof.nc")
    year_out = addfile("piControl_ep-cp_monsoon_and_ensoyear_highpass_bw_eof.nc", "c")

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