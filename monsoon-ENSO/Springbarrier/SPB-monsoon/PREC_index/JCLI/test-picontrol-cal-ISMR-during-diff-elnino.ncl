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

    ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,:,:),(/0,1/))
    p_month_str = dim_avg_n_Wrap(p_month_12(:,ind1,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/))
    
    var_ssta_str = var_get(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)  ;;; get variance
    var_p_str = var_get(dim_avg_n_Wrap(p_month_12(:,:,:,:),0), ind1)
    var_u_str = var_get(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    var_v_str = var_get(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,p_month_str,u_month_str,v_month_str,var_ssta_str,var_p_str,var_u_str,var_v_str/]
end 
;;*************************************************
;; main program 
;;*************************************************************
begin
    startmon = 101
    endmon = 99912
    startyear = 1
    endyear = 999

    lat1 = -50
    lat2 = 50
    lon1 = -180
    lon2 = 180

; ;;;读入海温数据
;     diri="/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
;     fils1=systemfunc("ls "+diri+"pr_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")
;     print(fils1)

;     f1=addfiles(fils1, "r")
;     date:=cd_calendar(f1[:]->time, 1)

;     lat = f1[0]->lat
;     lon = f1[0]->lon

;     timeind:=ind(date.le.endmon+200.and.date.ge.startmon)

;     sst=lonFlip(f1[:]->pr(timeind,{lat1:lat2},{lon1:lon2}))
;     sst = (/sst*86400/)
;     printVarSummary(sst)
    
;     ssta=(/rmMonAnnCycTLL(sst)/)
;     copy_VarCoords(sst, ssta)

;     ;;detrend 
;     ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
;     ;;;
;     delete(date)
;     delete(timeind)
    
;     ; nino34_12_r =  split_12(nino34)
;     ; nino34_djf = (/dim_avg_n_Wrap(nino34_12_r(5:8,:), 0)/)
;     ; std_nino = dim_stddev(nino34_djf)
;     sst_total_year = ispan(startyear, endyear+2, 1) ;;;0-199
;     print("************** finish SST data*****************")

;     ncout = addfile("pre_picontrol.nc", "c")
;     ncout->pre = ssta
; ;;;
      f1 = addfile("pre_picontrol.nc", "r")
      p = f1->pre 
;;;;;

; ;;读入ua
;     diri2 = "/home/yangsong3/CMIP6/2.5x2.5/CESM2/piControl/"
;     fils2 = systemfunc("ls "+diri2+"ua_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")
  
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
;     fils3 = systemfunc("ls "+diri2+"va_Amon_CESM2_piControl_r1i1p1f1_gn_0[0-9]"+"*")
  
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

    f2 = addfile("uv_picontrol.nc", "r")
    u_850 = f2->u_850
    v_850 = f2->v_850

;;;


    f3 = addfile("ssta_picontrol.nc", "r")
    ssta = f3->ssta
    


;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    year_out = addfile("piControl_monsoon_and_ensoyear.nc", "r")

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


    yearnum = endyear - startyear +1 -1
    dim_sst = dimsizes(ssta)
    dim_u = dimsizes(u_850)


;;;;;
   sst_year = ispan(startyear,endyear,1) 
   u_year = sst_year
;;;;

;;;split to 12 month    3月开始，第一年是1900年
  dim1  = dimsizes(ssta)
  yearnum = dim1(0)/12 - 1
  ssta_month_12 = new((/12,yearnum,dim1(1),dim1(2)/),float)

  dim2  = dimsizes(u_850)
  yearnum := dim2(0)/12 - 1 
  u_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)
  v_850_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)
  
  p_month_12 = new((/12,yearnum,dim2(1),dim2(2)/),float)

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



;;;******************************************************************
;;;;;开始画图
;;;******************************************************************

  wksstring = (/"El Nino","La Nina","Normal"/)
;;;;


  case1 = [/wet_warm_year,dry_warm_year,nor_warm_year/]   ;;;;;

  case2 = [/wet_cold_year,dry_cold_year,nor_cold_year/]

  case3 = [/wet_norm_year,dry_norm_year,nor_norm_year/]

  case = [/case1,case2,case3/]
  
 imsr_inten = new(6,float)


do k=0,2

  case_need:= case[k]
  ind_us1 :=  get1Dindex(u_year,case_need[0])
  ind_uw1 :=  get1Dindex(u_year,case_need[1])
  ind_unor1 :=  get1Dindex(u_year,case_need[2])

 ; month = (/"JJA","SON","DJF(+1)","MAM(+1)"/)  
  ; month = (/"MAM","JJA","SON","DJF(+1)","MAM(+1)","JJA(+1)"/) 

  month = (/"MAM","JJAS","OND","JF(+1)"/) 

  order = (/"(a)","(b)","(c)","(d)",\
           "(e)","(f)","(g)","(h)"/)

  month_str = ispan(1,8,1) + ""

  plot = new(dimsizes(month)*2,graphic)
  plot1 = new(dimsizes(month)*2,graphic)
  plot2 = new(dimsizes(month)*2,graphic)
  iplot = plot1
  iplot1 = plot
  iplot2 = plot
   
 

do i=0,dimsizes(month)-1
;;; strong and weak monsoon year p , u ,v 
   if i.eq.0 then 
    index1 = i*3
    index2 = i*3+2
    h = 0
   end if

   if i.eq.1 then 
    index1 = i*3
    index2 = i*3+3
    h = 0
   end if

   if i.eq.2 then 
    index1 = i*3+1
    index2 = i*3+3
    h = 0
   end if

   if i.eq.3 then 
    index1 = i*3+1
    index2 = i*3+2
    h = 0
   end if


    com_month_str  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_us1+h,ind_us1+h,ind_us1+h)
    
    ssta_month_str = com_month_str[0]
    p_month_str = com_month_str[1]
    u_month_str = com_month_str[2]
    v_month_str = com_month_str[3]
    
    var_ssta_str = com_month_str[4]
    var_p_str = com_month_str[5]
    var_u_str = com_month_str[6]
    var_v_str = com_month_str[7]

    com_month_weak  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                     v_850_month_12(index1:index2,:,:,:),ind_uw1+h,ind_uw1+h,ind_uw1+h)
    

    ssta_month_weak = com_month_weak[0]
    p_month_weak = com_month_weak[1]
    u_month_weak = com_month_weak[2]
    v_month_weak = com_month_weak[3]
    
    var_ssta_weak = com_month_weak[4]
    var_p_weak = com_month_weak[5]
    var_u_weak = com_month_weak[6]
    var_v_weak = com_month_weak[7]

    com_month_nor  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_unor1+h,ind_unor1+h,ind_unor1+h)
    
    ssta_month_nor = com_month_nor[0]
    p_month_nor = com_month_nor[1]
    u_month_nor = com_month_nor[2]
    v_month_nor = com_month_nor[3]
    
    var_ssta_nor = com_month_nor[4]
    var_p_nor = com_month_nor[5]
    var_u_nor = com_month_nor[6]
    var_v_nor = com_month_nor[7]
    ; p_month_str = where(var_p_str.le.0.1, p_month_str, p_month_str@_FillValue)
    ; p_month_weak = where(var_p_weak.le.0.1, p_month_weak, p_month_weak@_FillValue)
;;;;;;


;;;;get data 
;;;;;;

  u_month_str = (/u_month_str - u_month_nor/)
  v_month_str = (/v_month_str - v_month_nor/)
  p_month_str = (/p_month_str - p_month_nor/)
  ssta_month_str = (/ssta_month_str - ssta_month_nor/)

  u_month_weak = (/u_month_weak - u_month_nor/)
  v_month_weak = (/v_month_weak - v_month_nor/)
  p_month_weak = (/p_month_weak - p_month_nor/)
  ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)
  
;;;;;
if i.eq.1 then 
  p_month_str_indian =  avg(p_month_str({10:30},{60:100}))
  p_month_weak_indian =  avg(p_month_weak({10:30},{60:100}))
 print("*******************"+wksstring(k)+"******************************")
 print("Str_monsoon="+p_month_str_indian)
 print("Weak_monsoon="+p_month_weak_indian)
 
 imsr_inten(k*2) = p_month_str_indian
 imsr_inten(k*2+1) = p_month_weak_indian
end if 
end do
end do 

opt="ISMR_inten_diff_mons.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",imsr_inten))  


end 