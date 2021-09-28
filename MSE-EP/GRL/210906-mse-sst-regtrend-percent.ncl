function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data,0,0) 
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

function part_reg(var1,var2,var3)
begin 
  ;;;remove var3 influence 
  ;;;var1 for y var2 for x1 
  ;;;variance
  variance1 = dim_stddev_n_Wrap(var1, 0)
  ; reg = regCoef_n(var3, var2,0, 0)
  ; var2_left = var2
  ; var2_left = var2 - (reg*var3 +reg@yintercept)
  ; variance2 = dim_stddev_n_Wrap(var2_left, 0)
   variance2 = dim_stddev_n_Wrap(var2, 0)

  r12 = escorc_n(var1, var2, 0, 0)
  ; copy_VarCoords(var2(0,:,:), rc_12)  ;;;这是二维
  r13 = escorc_n(var1, var3, 0, 0)  ;;;这个是1dim
  r23 = escorc_n(var2, var3, 0, 0)   ;;;这个是2dim
  ; r12_3 = (/(r12 - r13*r23)/sqrt((1- r13^2)*(1- r23^2))/)
  r12_3 = (/(r12 - r13*r23)/(1- r23^2)/)
  rc_mam = r12
  rc_mam = (/r12_3*variance1/variance2/)
return rc_mam
end 

function normal_reg(var1,var2)
begin 
  ;;;var1 for x var2 for y
  ;;;remove var3 influence 
  ;;;variance
  rc_mam = regCoef_n(var1, var2, 0, 0)
return rc_mam
end 
function part_corr(var1,var2,var3)
begin 
  ;;;remove var3 influence 
  r12 = escorc_n(var1, var2, 0, 0)
  ; copy_VarCoords(var2(0,:,:), r12)  ;;;这是二维
  r13 = escorc_n(var1, var3, 0, 0)  ;;;这个是1dim
  r23 = escorc_n(var2, var3, 0, 0)   ;;;这个是2dim
  r12_3 = (/(r12 - r13*r23)/sqrt((1- r13^2)*(1- r23^2))/)
  rc_mam = r12
  rc_mam = (/r12_3/)
return rc_mam
end 

function part_reg_4(var1,var2,var3,var4)
begin 
;;;;去掉局地海温的影响,偏相关性
  variance1 = dim_stddev_n_Wrap(var1, 0)
  variance2 = dim_stddev_n_Wrap(var2, 0)
;;; 1 for MSE, 2 sst region , 3 for nino3_4 ,4 for atlantic SST
  r12no3 = part_corr(var1,var2,var3) ;;;12.3  2dim
  ; rc_mam_rmaltantic = part_corr(MSE_nep_mam,ssta_mam,nino3_4) ;;;12.4 
  r12no4 = part_corr(var1,var2, var4) ;;; 2dim

  r13no4 = part_corr(var1,var3, var4)  ;;; 1dim

  r23no4 = part_corr(var2,var3,var4)    ;;; 1dim
 
  r13no4_23no4 = r13no4*r23no4
  r13no4_23no4_2dim  = r12no4
  r13no4_23no4_2dim  = r13no4_23no4
  rc_mam = (/r12no4 - r13no4_23no4_2dim/)/(1-r23no4^2)
  rc_mam = rc_mam*(/variance1/variance2/)
return rc_mam
end 
begin
    ;;;time span
    startmon =198001
    endmon = 202012

    running_time = 15

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;; 
;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ssta = sst
    ssta=(/rmMonAnnCycTLL(sst)/)

    ssta_extr = ssta
  
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

    copy_VarCoords(sst, ssta)
   ;;
;;;
    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_ncep1_50-20.nc", "r")

    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    startmon = 198001
    endmon = 202012

    startyear = 1980
    endyear = 2020

    timeind :=ind(date.le.endmon.and.date.ge.startmon)

    MSE = f1->MSE(timeind,:,:,:)
   

    ;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)

    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    MSE_col_extr = MSE_col
    MSE_nep_extr = MSE_col_extr
    MSE_nep_extr = (/MSE_col_extr/10^6/)

    MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)
   

    ; MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_nep = MSE_col
    MSE_nep = (/MSE_col/10^6/)


;;;;;;split to 12 month  is a array size (12,yearnum)
    MSE_nep_12 = split_12(MSE_nep)
    ssta_12 = split_12(ssta) 

    ssta_extr_12 = split_12(ssta_extr) 
    MSE_nep_extr_12  = split_12(MSE_nep_extr) 

;;;;;

   
;;;年平均的数据
    MSE_nep_annual = month_to_annual(MSE_nep, 1)
    MSE_nep_annual!0 = "time"
    MSE_nep_annual&time = ispan(startyear,endyear,1)

    ssta_annual = month_to_annual(ssta, 1)
    ssta_annual!0 = "time"
    ssta_annual&time = ispan(startyear,endyear,1)
    
    ssta_annual = (/dim_avg_n(ssta_12(0:4,:,:,:),0)/)
    MSE_nep_annual = (/dim_avg_n(MSE_nep_12(0:4,:,:,:),0)/)

    ssta_mam = ssta_annual
    

    ssta_annual_expr = ssta_annual
    ssta_annual_expr = (/dim_avg_n_Wrap(ssta_extr_12(0:4,:,:,:),0)/)

    MSE_nep_annual_expr = MSE_nep_annual
    MSE_nep_annual_expr = (/dim_avg_n_Wrap(MSE_nep_extr_12(0:4,:,:,:),0)/)

;;;;;

;;;
 
    MSE_nep_mam = dim_avg_n_Wrap(MSE_nep_annual(:,{0:20},{250:280}), (/1,2/))
    MSE_nep_annual_expr_mean = dim_avg_n_Wrap(MSE_nep_annual_expr(:,{0:20},{250:280}), (/1,2/))
    ;;;局地海温
    ; ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:20},{240:280}), (/1,2/))

    ;  ;;nino1+2:
    ;   ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{270:280}), (/1,2/))
    ; ;;;
    
    ;  ;;nino3:
    ;   ssta_mam_nep = dim_avg_n_Wrap(ssta_mam(:,{0:-10},{150:270}), (/1,2/))
    ; ;;;
     
     ;;;nino4 
      nino3_4 = dim_avg_n_Wrap(ssta_mam(:,{-5:5},{190:240}), (/1,2/))
      ssta_mam_atlantic = dim_avg_n_Wrap(ssta_mam(:,{10:30},{280:330}), (/1,2/))

      ssta_mam_indian = dim_avg_n_Wrap(ssta_mam(:,{-10:30},{40:110}), (/1,2/))

      ssta_mam_local = dim_avg_n_Wrap(ssta_mam(:,{0:20},{250:280}), (/1,2/))
    ;;;;

    ; ;;;;去掉ENSO海温的影响,偏相关性
      rc_mam_rmenso_indian = part_reg(MSE_nep_mam,ssta_mam_indian, nino3_4) ;;;12.3
      ; rc_mam_rmenso_atlantic = part_reg(MSE_nep_mam, ssta_mam_atlantic,nino3_4) ;;;12.3
    ;   rc_mam_rmenso_local = part_reg(MSE_nep_mam,ssta_mam_local, nino3_4) ;;;12.3
    ; ;;;

       ;;;;去掉ENSO海温的影响同时去掉Atlantic偏相关性
      rc_mam_rmenso_atlantic = part_reg_4(MSE_nep_mam, ssta_mam_atlantic,nino3_4,ssta_mam_indian) ;;;12.3
      ; rc_mam_rmenso_indian = part_reg_4(MSE_nep_mam,ssta_mam_indian, nino3_4,ssta_mam_atlantic) ;;;12.3
      ; rc_mam_rmenso_atlantic = part_reg_4(MSE_nep_mam, ssta_mam_atlantic,nino3_4,ssta_mam_atlantic) ;;;12.3
      rc_mam_rmenso_local = part_reg_4(MSE_nep_mam,ssta_mam_local, nino3_4,ssta_mam_indian) ;;;12.3
    ;;;

    ;;;;

    ; ;;没有偏相关的情况
    ;   rc_mam_rmenso_indian = normal_reg(ssta_mam_indian,MSE_nep_mam) ;;;12.3
    ;   rc_mam_rmenso_atlantic = normal_reg(ssta_mam_atlantic,MSE_nep_mam) ;;;12.3
    ;   rc_mam_rmenso_local = normal_reg(ssta_mam_local,MSE_nep_mam) ;;;12.3
    ; ;;;
    
      ssta_mam_atlantic_expr = dim_avg_n_Wrap(ssta_annual_expr(:,{10:30},{280:330}), (/1,2/))

      ssta_mam_indian_expr = dim_avg_n_Wrap(ssta_annual_expr(:,{-10:30},{40:110}), (/1,2/))

      ssta_mam_local_expr = dim_avg_n_Wrap(ssta_annual_expr(:,{0:20},{250:280}), (/1,2/))

    ;;;;linear trend 
      
      trend_mse  = mse_trend(MSE_nep_annual_expr_mean,startyear,endyear)

      trend_indian  = mse_trend(ssta_mam_indian_expr,startyear,endyear)

      trend_atlantic  = mse_trend(ssta_mam_atlantic_expr,startyear,endyear)

      trend_local  = mse_trend(ssta_mam_local_expr,startyear,endyear)
    ;;;;


     mse_trend_indian = rc_mam_rmenso_indian*trend_indian

     mse_trend_atlantic = rc_mam_rmenso_atlantic*trend_atlantic

     mse_trend_local = rc_mam_rmenso_local*trend_local


     print("mse_trend="+trend_mse)

     print("mse_trend_atlantic="+mse_trend_atlantic)

     print("mse_trend_local="+mse_trend_local)

     print("mse_trend_indian="+mse_trend_indian)


     print("percent_atlantic="+mse_trend_atlantic/trend_mse)
     print("percent_indian="+mse_trend_indian/trend_mse)
     print("percent_local="+mse_trend_local/trend_mse)

 end 