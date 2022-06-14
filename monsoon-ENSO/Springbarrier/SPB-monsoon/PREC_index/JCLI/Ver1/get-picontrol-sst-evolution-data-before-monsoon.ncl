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
  

    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999
   
    ; year_out = addfile("../60-100_piControl_monsoon_and_ensoyear.nc", "r")
    year_out = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/piControl_monsoon_and_ensoyear_highpass_bw.nc", "r")

    nino_year = year_out->nino_year 
    lanina_year = year_out->lanina_year 
    normal_year = year_out->normal_year 

    monstryear = year_out->monstryear 
    monwekyear = year_out->monwekyear
    noryear = year_out->noryear
 
    wet_warm_year =  year_out->wet_warm_year
    wet_cold_year = year_out->wet_cold_year
    wet_norm_year = year_out->wet_norm_year

    dry_warm_year = year_out->dry_warm_year
    dry_cold_year = year_out->dry_cold_year
    dry_norm_year = year_out->dry_norm_year

    nor_warm_year = year_out->nor_warm_year
    nor_cold_year = year_out->nor_cold_year
    nor_norm_year = year_out->nor_norm_year

    ; f3 = addfile("../ssta_picontrol.nc", "r")
    f3 = addfile("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/ssta_picontrol_highpass_bw.nc", "r")
    ssta = f3->ssta
    

;;;;;
   sst_year = ispan(startyear,endyear,1) 
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]

  stringlesft = (/"(a) El Nino","(b) La Nina","(c) Normal"/)
   str1= (/"El Nino","La Nina","Normal"/)



do k=0,2 
  print("***********start*********************"+k)

   case_need:= case[k]

   ind_ts1_r :=  get1Dindex(sst_year,case_need[0])
   ind_tw1_r :=  get1Dindex(sst_year,case_need[1])
   ind_tnor1_r :=  get1Dindex(sst_year,case_need[2])

   ind_tnor1 := ind_tnor1_r(1:dimsizes(ind_tnor1_r)-1)
   ind_tw1 := ind_tw1_r(1:dimsizes(ind_tw1_r)-1)
   ind_ts1 := ind_ts1_r(1:dimsizes(ind_ts1_r)-1)

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

    com_month_str  = composite_mean(ssta_month_12,ind_ts1-1)
    ssta_str_evo(0:11) = com_month_str[0]
  
    com_month_weak  = composite_mean(ssta_month_12,ind_tw1-1)
    ssta_weak_evo(0:11) = com_month_weak[0]
    
   
   

    com_month_nor  = composite_mean(ssta_month_12,ind_tnor1_r-1)
   
    ssta_nor_evo(0:11) = com_month_nor[0]


    com_month_str  = composite_mean(ssta_month_12(0:5,:,:,:),ind_ts1)
    
    ssta_str_evo(12:) = com_month_str[0]
  

    com_month_weak  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tw1)
    
    ssta_weak_evo(12:) = com_month_weak[0]
    
    
    com_month_nor  = composite_mean(ssta_month_12(0:5,:,:,:),ind_tnor1_r)
    
    ssta_nor_evo(12:) = com_month_nor[0]

    ;;;;写入数据
    ; ssta_month_str = where(prob_ssta_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)
    ; ssta_month_weak = where(prob_ssta_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)
    data := (/ssta_str_evo,ssta_weak_evo,ssta_nor_evo/)
    
    ; opt1="60-100_highpass_"+str1(k)+""+"diffmonsoon-SSTevolution.txt"
    opt1="highpass_bw_beforemonsoon_"+str1(k)+""+"diffmonsoon-SSTevolution.txt"
    ; write_matrix(data, "12f7.2",opt1)
    asciiwrite (opt1, data)

    
end do 

end 