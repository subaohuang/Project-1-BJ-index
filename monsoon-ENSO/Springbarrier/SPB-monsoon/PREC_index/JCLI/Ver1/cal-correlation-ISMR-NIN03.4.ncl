;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var = dim_variance_n(pv_ano(east_time,:,:), 0)
    var  =where(var.eq.0, var@_FillValue, var)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
    var@dim = n
 return var
 end
function composite_mean(ssta_month_12,u_850_month_12,v_850_month_12,ind1,ind2,ind3)
begin 

    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/))

    prob_ssta_str = one_test(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)
    prob_u_str = one_test(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    prob_v_str = one_test(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,u_month_str,v_month_str,prob_ssta_str,prob_u_str,prob_v_str/]
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
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99712
    startyear = 1
    endyear = 997

    lat1 = -50
    lat2 = 50
    lon1 = -180
    lon2 = 180

; ;;;读入海温数据
    ; diri="/home/ys17-23/CMIP6/1.0x1.0/CESM2/piControl/Omon/"
    ; fils1=systemfunc("ls "+diri+"tos_Omon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")

    ; f1=addfiles(fils1, "r")
    ; date:=cd_calendar(f1[:]->time, 1)

    ; lat = f1[0]->lat
    ; lon = f1[0]->lon

    ; timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

    ; sst=lonFlip(f1[:]->tos(timeind,{lat1:lat2},{lon1:lon2}))
    ; ; sst = (/sst*86400/)
    ; printVarSummary(sst)
    
    ; ssta=(/rmMonAnnCycTLL(sst)/)
    ; copy_VarCoords(sst, ssta)

    ; ;;detrend 
    ; ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ; ;;;
    ; delete(date)
    ; delete(timeind)
    
    ; ; nino34_12_r =  split_12(nino34)
    ; ; nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)
    ; ; std_nino = dim_stddev(nino34_djf)
    ; sst_total_year = ispan(startyear, endyear+2, 1) ;;;0-199
    ; print("************** finish SST data*****************")

    ; ncout = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/JCLI/ssta_picontrol.nc", "c")
    ; ncout->ssta = ssta
    ; exit()
; ; ;;;
  f1 = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/ssta_picontrol.nc", "r")
  nino34 = dim_avg_n_Wrap(f1->ssta(:,{-5:5},{190:240}),(/1,2/))

  nino34_12_r =  split_12(nino34)
  nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)

;;;;;;

; ;;读入ua
;     diri2 = "/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
;     fils2 = systemfunc("ls "+diri2+"ua_Amon_CESM2_piControl_r1i1p1f1_gn_0{0..9}"+"*")
  
;     f2 = addfiles(fils2, "r")
;     date:=cd_calendar(f2[:]->time, 1)

;     lat = f2[0]->lat
;     lon = f2[0]->lon

;     timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

;     u_850 = lonFlip(f2[:]->ua(timeind,{85000},{lat1:lat2},{lon1:lon2}))
;     printVarSummary(u_850)
 
; ;;;;

; ;;读入va
;     diri3 = "/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
;     fils3 = systemfunc("ls "+diri2+"va_Amon_CESM2_piControl_r1i1p1f1_gn_0{0..9}"+"*")
  
;     f3 = addfiles(fils3, "r")
;     date:=cd_calendar(f3[:]->time, 1)

;     lat = f3[0]->lat
;     lon = f3[0]->lon

;     timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

;     v_850 = lonFlip(f3[:]->va(timeind,{85000},{lat1:lat2},{lon1:lon2}))
;     printVarSummary(v_850)
; ;;;;
;     u_850 = (/rmMonAnnCycTLL(u_850)/)
;     v_850 = (/rmMonAnnCycTLL(v_850)/)

;     u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
;     v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)
;     ncout = addfile("uv_picontrol.nc", "c")
;     ncout->u_850 = u_850
;     ncout->v_850 = v_850
; ; ;;;;
    f2 = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/pre_picontrol.nc", "r")
    
  pindian = dim_avg_n_Wrap(f2->pre(:,{10:30},{60:100}),(/1,2/))

  pindian_12_r =  split_12(pindian)
  pindian_djf = (/dim_avg_n_Wrap(pindian_12_r(5:8,:), 0)/)

;;;;

corr = escorc_n(nino34_djf, pindian_djf, 0, 0)
print(corr)
exit()
  

end 

