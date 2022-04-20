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

;;;
    yearstart = 195001
    yearend = 200412
    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)
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
      ; ssta_index = dim_avg_n_Wrap(ssta(:,{-5:5},{250:280}), (/1,2/)) ;;nino3


      printVarSummary(ssta_index)

      obs_nino34_r = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
      
      obs_nino34 = obs_nino34_r(1:)
      printVarSummary(obs_nino34)


      ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1
      
      ; obs_nino34_std = dim_stddev_n_Wrap(obs_nino34, 0)
      ; ind_obs_nino = ind(obs_nino34.ge.obs_nino34_std)

      ind_obs_nino = ely - 1950

      obs_nino_intensity = avg(obs_nino34(ind_obs_nino))

      print("obs_nino_intensity="+obs_nino_intensity)



do k = 0,dim_model-1
;;;;

    model = model_name(k)


    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    ; timeind:=ind(date.le.200412.and.date.ge.195001)

    delete(date)     
    tos = lonFlip(f2[:]->tos(timeind,:,:))
    printVarSummary(tos)

   
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)

    tos_clm  = dim_avg_n_Wrap(tos, 0)
    
    ;;;大西洋和太平洋的海温差

    tos_diff  = avg(tos_clm({0:-20},{270:280})) - avg(tos_clm({0:20},{290:330}))

    enino_intensity(k) = tos_diff
    ;;;
end do 
    asciiwrite("SST-DIFF-PAC-ATL_CMIP5", enino_intensity)
    
  
end 