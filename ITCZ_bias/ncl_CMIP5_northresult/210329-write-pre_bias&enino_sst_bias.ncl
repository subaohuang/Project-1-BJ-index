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
    p=f1->precip(timeind,{-30:30},:)
   
    delete(p&time)
    p&time=date(timeind)
    
    p_annual = month_to_season(p, "DJF")
    ; p_annual = month_to_annual(p, 1)

    p_annual = (/dtrend_msg_n(ispan(1,dimsizes(p_annual&time),1),p_annual,False,True,0)/)
    p_clm = dim_avg_n_Wrap(p_annual, 0) 
    printVarSummary(p_clm)
    printMinMax(p_clm, 1)

   dim = dimsizes(p_clm)
;;;



; ;;;读入HadiSST海温数据
      diri="/home/yangsong3/data-observation/linshh/data/sst/"
      fils1=systemfunc("ls "+diri+"HadISST_sst_2.5.nc")
      f1=addfile(fils1, "r")
      date:=cd_calendar(f1->time, 1)

      timeind:=ind(date.le.201912.and.date.ge.195001)

      sst =  lonFlip(f1->sst(timeind,{-30:30},:))

      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)
      printVarSummary(ssta)
     

      
      ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)

      ssta_djf_r = month_to_season(ssta, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
      
      ssta_djf = ssta_djf_r(1:,:,:)  ;;;;50年的DJF
      printVarSummary(ssta_djf)


      ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1
      
      ; obs_nino34_std = dim_stddev_n_Wrap(obs_nino34, 0)
      ; ind_obs_nino = ind(obs_nino34.ge.obs_nino34_std)

      ind_obs_nino = ely - 1950

      obs_nino_intensity = dim_avg_n_Wrap(ssta_djf(ind_obs_nino,:,:), 0)

    
      dim_sst = dimsizes(sst(1,:,:))
;;;
    yearstart = 195001
    yearend = 200412
    pre_bias = new((/dim_model,dim(0),dim(1)/),float)
    enino_bias = new((/dim_model,dim_sst(0),dim_sst(1)/),float)
;;;

   
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
    pr = lonFlip(f1->pr(timeind,{-30:30},:))

    pr = (/dtrend_msg_n(ispan(1,dimsizes(pr&time),1),pr,False,True,0)/)
    
    pr = (/pr*86400/) ;;conver to mm/day
    printVarSummary(pr)
    printMinMax(pr, 1)

    pr_annual = month_to_season(pr, "DJF")
    ; pr_annual = month_to_annual(pr, 1)

    PR_clm = dim_avg_n_Wrap(pr_annual, 0)

    pre_diff = PR_clm
    pre_diff = (/PR_clm - p_clm/) 

    pre_bias(k,:,:) = pre_diff
    

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    ; timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    timeind:=ind(date.le.200412.and.date.ge.195001)

    delete(date)     
    tos = lonFlip(f2[:]->tos(timeind,{-30:30},:))
    printVarSummary(tos)

   
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos = (/rmMonAnnCycTLL(tos)/)
    tos_djf = month_to_season(tos, "DJF")

    nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{190:240}), (/1,2/))  ;;;nino3_4
    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{210:270}), (/1,2/)) ;;nino3
     ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-10:0},{270:280}), (/1,2/)) ;;nino1-2
    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{250:280}), (/1,2/)) ;;nino3

    nino3_4_std = dim_stddev_n_Wrap(nino3_4, 0)



    yearind := ind(nino3_4.ge.0.8*nino3_4_std)


    tos_nino := dim_avg_n_Wrap(tos_djf(yearind,:,:),0)
    enino_bias(k,:,:) = tos_nino
    enino_bias(k,:,:) = (/tos_nino - obs_nino_intensity/)
  
     ;;;
end do 
enino_bias!0 = "model"
pre_bias!0 = "model"
printVarSummary(enino_bias)
printVarSummary(pre_bias)

    system("/bin/rm -f pre_and_enino_bias.nc")
    out = addfile("DJF_pre_and_enino_bias.nc", "c")
    out->pre_bias = pre_bias
    out->enino_bias = enino_bias



end 