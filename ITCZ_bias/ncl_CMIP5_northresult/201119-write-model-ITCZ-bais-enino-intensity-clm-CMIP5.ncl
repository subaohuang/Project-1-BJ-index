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


begin
;读入数据

    read_name = asciiread("/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/info.txt", -1,"string")
    print(read_name)
    model_name = read_name 

    dim_model = dimsizes(model_name)

    plot = new(dim_model,graphic)
    iplot= plot
    plot1 = plot

;;;读入GPCP的数据,计算其春季气候态
    f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.200412)
    p=f1->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    
    p_annual = month_to_season(p, "MAM")
    p_annual = (/dtrend_msg_n(ispan(1,dimsizes(p_annual&time),1),p_annual,False,True,0)/)
    p_clm = dim_avg_n_Wrap(p_annual, 0) 
    printVarSummary(p_clm)
    printMinMax(p_clm, 1)
;;;


;;;
    yearstart = 195001
    yearend = 200412
    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)
    lanina_intensity = new(dim_model,float)
;;;



; ;;;读入HadiSST海温数据
      diri="/home/yangsong3/data-observation/linshh/data/sst/"
      fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
      f1=addfile(fils1, "r")
      date:=cd_calendar(f1->time, 1)

      timeind:=ind(date.le.201912.and.date.ge.195001)

      sst =  lonFlip(f1->sst(timeind,:,:))

      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)

      
      ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

      ; ssta_index = dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
      ssta_index = dim_avg_n_Wrap(ssta(:,{-5:5},{210:270}), (/1,2/)) ;;nino3
      ; ssta_index = dim_avg_n_Wrap(ssta(:,{-5:5},{160:210}), (/1,2/)) ;;nino3
      ; ssta_index = dim_avg_n_Wrap(ssta(:,{-5:5},{250:280}), (/1,2/)) ;;nino3


      printVarSummary(ssta_index)

      obs_nino34_r = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
      
      obs_nino34 = obs_nino34_r(1:)
      printVarSummary(obs_nino34)
      ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1
      
      obs_nino34_std = dim_stddev_n_Wrap(obs_nino34_r, 0)
      ; ind_obs_nino = ind(obs_nino34.ge.obs_nino34_std)
      ; ind_obs_nino = ely - 1950
      ; obs_nino_intensity = avg(obs_nino34(ind_obs_nino))
       
      ind_obs_nino = ind(obs_nino34_r.ge.0.5*obs_nino34_std)
      ind_obs_lanina = ind(obs_nino34_r.le.-0.5*obs_nino34_std)

      obs_nino_intensity = avg(obs_nino34_r(ind_obs_nino))
      obs_lanina_intensity = avg(obs_nino34_r(ind_obs_lanina))
       
    
   
do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
   
    fils1 =  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/pr_chazhi/pr_Amon_"+model+"_historical_r1i1p1_"+"*")

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    pr = lonFlip(f1->pr(timeind,:,:))

    pr = (/dtrend_msg_n(ispan(1,dimsizes(pr&time),1),pr,False,True,0)/)
    
    pr = (/pr*86400/) ;;conver to mm/day
    printVarSummary(pr)
    printMinMax(pr, 1)

    pr_annual = month_to_season(pr, "MAM")
    PR_clm = dim_avg_n_Wrap(pr_annual, 0)

    diff = PR_clm
    diff = (/PR_clm - p_clm/) 

    ; ;;;0 -20 100-150
    ; siindex = avg(diff({0:-20},{100:150}))
    ; ;;

    ;;;0 -20 150-
    siindex =  avg(diff({0:10},{210:280}))
    ;;;0 -20 150-
    ; siindex =  avg(diff({0:-10},{210:280})) - avg(diff({0:10},{210:280}))
    ; siindex =  avg(diff({0:-10},{210:280})) 
    ; siindex = avg(diff({10:-10},{180:270}))
    ; siindex = avg(diff({0:-30},{210:260}))
    ;;
    dbitcz(k) = siindex 

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    ; timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    timeind:=ind(date.le.200412.and.date.ge.195001)

    delete(date)     
    tos = lonFlip(f2[:]->tos(timeind,:,:))
    printVarSummary(tos)

   
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos = (/rmMonAnnCycTLL(tos)/)
    tos_djf = month_to_season(tos, "DJF")

    tos_djf(0:dimsizes(tos&time)/12 -2,:,:) = tos_djf(1:,:,:)
    tos_djf(dimsizes(tos&time)/12 -1,:,:) = tos(dimsizes(tos&time)-1,:,:)

    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{190:240}), (/1,2/))  ;;;nino3_4
    nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{210:270}), (/1,2/)) ;;nino3
    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{160:210}), (/1,2/)) ;;nino3
     ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-10:0},{270:280}), (/1,2/)) ;;nino1-2
    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{250:280}), (/1,2/)) ;;nino3
    nino3_4_std = dim_stddev_n_Wrap(nino3_4, 0)

    print(nino3_4_std)

    yearind := ind(nino3_4.ge.0.5*nino3_4_std)
    yearlaind := ind(nino3_4.le.-0.5*nino3_4_std)

    tos_nino := nino3_4(yearind)
    nino_intensity = avg(tos_nino)

    tos_lanina := nino3_4(yearlaind)
    la_intensity = avg(tos_lanina)


    enino_intensity(k) = (/nino_intensity - obs_nino_intensity/)
    lanina_intensity(k) = (/la_intensity - obs_lanina_intensity/)


end do 

    ; asciiwrite("Nino4_bais_CMIP5", enino_intensity)
    ; asciiwrite("Nino4_bais_CMIP5", enino_intensity)
    asciiwrite("Nino3_bais_0.5std_CMIP5", enino_intensity)
    asciiwrite("Lanina3_bais_0.5std_CMIP5", lanina_intensity)


    ; asciiwrite("enino_intensity_30model", enino_intensity)
    ; asciiwrite("DJF-itcz_bais(0-10_210-280)_CMIP5", dbitcz)
    ; exit()
end 