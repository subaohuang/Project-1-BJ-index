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
    p_month_str = dim_avg_n_Wrap(p_month_12(:,ind2,:,:),(/0,1/))
    u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,:,:),(/0,1/))
    v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,:,:),(/0,1/))
    
    var_ssta_str = var_get(dim_avg_n_Wrap(ssta_month_12(:,:,:,:),0), ind1)  ;;; get variance
    var_p_str = var_get(dim_avg_n_Wrap(p_month_12(:,:,:,:),0), ind2)
    var_u_str = var_get(dim_avg_n_Wrap(u_850_month_12(:,:,:,:),0), ind2)
    var_v_str = var_get(dim_avg_n_Wrap(v_850_month_12(:,:,:,:),0), ind3)

return [/ssta_month_str,p_month_str,u_month_str,v_month_str,var_ssta_str,var_p_str,var_u_str,var_v_str/]
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
    printVarSummary(sst)
;;;;
    ; startmon = 187101
    ; endmon = 201512
    ; startyear = 1871
    ; endyear = 2015

    ; fils1 = "/home/ys17-23/Extension/20thCRV3/air.2m.mon.mean.nc"
    ; f1 = addfile(fils1, "r")
    ; time = f1->time
    ; date=cd_calendar(time, 1)

    ; endmon = endmon                 ;;; 1871-2015
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
     ssta = (/bw_bandpass(ssta)/)
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
    ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/ys17-23/Extension/20thCRV3/prate.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1871-2015
    startmon = startmon 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  

;;;;prate
    p  = f1->prate(timeind1,:,:)
    p = (/p*86400/) ;; mm/day
    p = (/rmMonAnnCycTLL(p)/)

    ;;detrend 
    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
     ;;;;filter 
     ; p = (/Lanczos_highpass(p)/)
      p = (/bw_bandpass(p)/)
    ;;;
    ;;;

    p_indian = dim_avg_n_Wrap(p(:,{10:30},{70:90}), (/1,2/))
    p_indian_jjas = month_to_season(p_indian, "JJA")
    p_indian_12 = split_12(p_indian)
    p_indian_jjas = (/dim_avg_n_Wrap(p_indian_12(5:8,:),0)/)
;;;



;;;AISMR DATA 

    ncols = 18 ;;;":YEAR JAN FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC    JF   MAM  JJAS   OND   ANN"
    nrows=146
    ; data = readAsciiTable("All-indian-monsoon-rainfall-index.txt", ncol, "float", 2)
    data = asciiread("/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/All-indian-monsoon-rainfall-index.txt", (/nrows,ncols/), "float")
    jjas = data(:,15) ;;sum of jjas precipitation  1871-2016
    
    jjas!0 = "time"
    jjas&time = ispan(1871,2016,1)

    pre_jjas = jjas({startyear:2015})  ;;;1871-2015

    pre_jjas = pre_jjas*0.1 ;;convert to mm
    clm = dim_avg_n_Wrap(pre_jjas, 0)
    pre_jjasa = dim_rmvmean(pre_jjas)
    pre_jjasa = (/dtrend_msg_n(ispan(1,2015-startyear+1,1),pre_jjasa,False,True,0)/)

    ; pre_jjasa = (/Lanczos_highpass(pre_jjasa)/)

    ; ;;;!!!! alert for change to pre_jjasa
    ; p_indian_jjas = (/pre_jjasa/)

;;;
 

;;;读入ERA wind data
    ; fils1 = "/home/ys17-23/Extension/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    fils1 = "/home/ys17-23/Extension/20thCRV3/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1871-2015
    startmon = startmon 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  
;;;;
 
;;;;U,V

    ; u_200  = short2flt(f1->u(timeind1,{200},:,:))
    ; u_850  = short2flt(f1->u(timeind1,{850},:,:))

    u_200  = f1->uwnd(timeind1,{200},:,:)
    u_850  = f1->uwnd(timeind1,{850},:,:)

    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;


;;;
    fils2 = "/home/ys17-23/Extension/20thCRV3/vwnd.mon.mean.nc"
    f2 = addfile(fils2, "r")
    time2 = f2->time
    date2 = cd_calendar(time2, 1)

    timeind2 = ind(date2.ge.startmon.and.date2.le.endmon)  ;;;风场数据从1979年开始

;;;;U,V
    v_200  = f2->vwnd(timeind1,{200},:,:)
    v_850  = f2->vwnd(timeind1,{850},:,:)

    v_850 = (/rmMonAnnCycTLL(v_850)/)
    v_200 = (/rmMonAnnCycTLL(v_200)/)
   
    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)

  ;;;filter 
  
   u_850 = (/bw_bandpass(u_850)/)
   v_850 = (/bw_bandpass(v_850)/)

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

;;;;;
    sst_year = ispan(startyear,endyear,1) 
    u_year = ispan(startyear,2015,1) 
;;;;;;;

;;;mon_elnino_years 

    DMI_ctrl_JJA_r = pre_jjasa 
    total_year = ispan(startyear,2014,1)   ;;;1872-2014

    DMI_ctrl_JJA = DMI_ctrl_JJA_r(0:dimsizes(DMI_ctrl_JJA_r)-2)  ;;;1871-2014
    std = dim_stddev_Wrap(DMI_ctrl_JJA_r)
;;;;

;;;nine composite  wet-warm, wet-cold, wet-normal, dry-warm, dry-cold, dry-normal, normal-warm, normal-cold, normal-normal
    nino34_djf_2 = nino34_djf(0:dimsizes(nino34_djf)-1-(endyear-2014))

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
    DMI_ctrl_JJA_r = p_indian_jjas
    DMI_ctrl_JJA = DMI_ctrl_JJA_r(0:dimsizes(DMI_ctrl_JJA_r)-2)  ;;;1871-2014
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
  
 imsr_inten = new(6,float)
 sst_inten = new(6,float)
 u_inten = new(6,float)


do k=0,2

  case_need:= case[k]
  ind_ts1 :=  get1Dindex(u_year,case_need[0])
  ind_tw1 :=  get1Dindex(u_year,case_need[1])
  ind_tnor1 :=  get1Dindex(u_year,case_need[2])

  case_need2:= case_20cr[k]
  ind_us1 :=  get1Dindex(u_year,case_need2[0])
  ind_uw1 :=  get1Dindex(u_year,case_need2[1])
  ind_unor1 :=  get1Dindex(u_year,case_need2[2])

  ; wks = gsn_open_wks("pdf", "highpass-Combinedindex-20CRV3-diff-9month-seasonalmean-UVPRE-Str-normon-"+" "+wksstring(k))
  ; wks1 = gsn_open_wks("pdf", "highpass-Combinedindex-20CRV3-diff-9month-seasonalmean-UVPRE-Weak-normon-"+" "+wksstring(k))  

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
                                    v_850_month_12(index1:index2,:,:,:),ind_ts1+h,ind_us1+h,ind_us1+h)
    
    ssta_month_str = com_month_str[0]
    p_month_str = com_month_str[1]
    u_month_str = com_month_str[2]
    v_month_str = com_month_str[3]
    
    var_ssta_str = com_month_str[4]
    var_p_str = com_month_str[5]
    var_u_str = com_month_str[6]
    var_v_str = com_month_str[7]

    com_month_weak  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                     v_850_month_12(index1:index2,:,:,:),ind_tw1+h,ind_uw1+h,ind_uw1+h)
    

    ssta_month_weak = com_month_weak[0]
    p_month_weak = com_month_weak[1]
    u_month_weak = com_month_weak[2]
    v_month_weak = com_month_weak[3]
    
    var_ssta_weak = com_month_weak[4]
    var_p_weak = com_month_weak[5]
    var_u_weak = com_month_weak[6]
    var_v_weak = com_month_weak[7]

    com_month_nor  = composite_mean(ssta_month_12(index1:index2,:,:,:),p_month_12(index1:index2,:,:,:),u_850_month_12(index1:index2,:,:,:),\
                                    v_850_month_12(index1:index2,:,:,:),ind_tnor1+h,ind_unor1+h,ind_unor1+h)
    
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

  ; u_month_str = (/u_month_str - u_month_nor/)
  ; v_month_str = (/v_month_str - v_month_nor/)
  ; p_month_str = (/p_month_str - p_month_nor/)
  ; ssta_month_str = (/ssta_month_str - ssta_month_nor/)

  ; u_month_weak = (/u_month_weak - u_month_nor/)
  ; v_month_weak = (/v_month_weak - v_month_nor/)
  ; p_month_weak = (/p_month_weak - p_month_nor/)
  ; ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)
  

  if k.ne.2 then 

    ; p_month_str = (/p_month_str - p_month_nor/)
    ; p_month_weak = (/p_month_weak - p_month_nor/)

    ssta_month_str = (/ssta_month_str - ssta_month_nor/)
    ssta_month_weak = (/ssta_month_weak - ssta_month_nor/)
    
    u_month_str = (/u_month_str - u_month_nor/)
    u_month_weak = (/u_month_weak - u_month_nor/)

  end if 

;;;;;
if i.eq.1 then 
  p_month_str_indian =  avg(p_month_str({10:30},{70:90}))
  p_month_weak_indian =  avg(p_month_weak({10:30},{70:90}))

  u_month_str_elnino =  avg(u_month_str({-5:5},{130:160}))
  u_month_weak_elnino   =  avg(u_month_weak({-5:5},{130:160}))

    ; u_month_str_elnino =  avg(u_month_str({0:10},{130:160}))
    ; u_month_weak_elnino   =  avg(u_month_weak({0:10},{130:160}))

    ; p_month_str_indian =  avg(p_month_str({10:30},{60:100}))
    ; p_month_weak_indian =  avg(p_month_weak({10:30},{60:100}))

    print("*******************"+wksstring(k)+"******************************")
    print("Str_monsoon="+p_month_str_indian)
    print("Weak_monsoon="+p_month_weak_indian)

    imsr_inten(k*2) = p_month_str_indian
    imsr_inten(k*2+1) = p_month_weak_indian

    u_inten(k*2) = u_month_str_elnino
    u_inten(k*2+1) = u_month_weak_elnino
end if 

if i.eq.2 then 
    ssta_month_str_elnino =  avg(ssta_month_str({-5:5},{190:240}))
    ssta_month_weak_elnino   =  avg(ssta_month_weak({-5:5},{190:240}))
    ; p_month_str_indian =  avg(p_month_str({10:30},{60:100}))
    ; p_month_weak_indian =  avg(p_month_weak({10:30},{60:100}))

    print("*******************"+wksstring(k)+"******************************")
    print("Str_monsoon="+ssta_month_str_elnino)
    print("Weak_monsoon="+ssta_month_weak_elnino)

    sst_inten(k*2) = ssta_month_str_elnino
    sst_inten(k*2+1) = ssta_month_weak_elnino
end if 

end do
end do 

; opt="OBS-ISMR_inten_diff_mons.ascii"    
; asciiwrite(opt, sprintf ("%6.2f,",imsr_inten))

; opt="OBS-SST_inten_diff_mons.ascii"    
; asciiwrite(opt, sprintf ("%6.2f,",sst_inten))

opt="OBS-u_inten_diff_mons.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",u_inten))

end 