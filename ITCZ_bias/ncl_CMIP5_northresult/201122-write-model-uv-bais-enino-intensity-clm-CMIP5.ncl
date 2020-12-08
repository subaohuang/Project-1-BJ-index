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



;;;读入v的数据,计算其春季气候态
    f1=addfile("/home/yangsong3/data-observation/ERA-interim/monthly/vwnd_2.5.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.200412)
    vwnd = lonFlip(short2flt(f1->v(timeind,{1000},:,:)))
    printVarSummary(vwnd)

    
    vwnd_annual = month_to_season(vwnd, "MAM")
    vwnd_annual = (/dtrend_msg_n(ispan(1,dimsizes(vwnd_annual&time),1),vwnd_annual,False,True,0)/)
    vwnd_clm = dim_avg_n_Wrap(vwnd_annual, 0) 
    printVarSummary(vwnd_clm)
    printMinMax(vwnd_clm, 1)
;;;


;;;
    yearstart = 197901
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
      ; ssta_index = dim_avg_n_Wrap(ssta(:,{-5:5},{210:270}), (/1,2/)) ;;nino3
      ssta_index = dim_avg_n_Wrap(ssta(:,{-5:5},{250:280}), (/1,2/)) ;;nino3


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

    fils1 := systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/v_chazhi/va_Amon_"+model+"_historical_"+"*")
   print(fils1)

    if dimsizes(fils1).ne.1 then 
      f1 := addfiles(fils1, "r")
      time:=f1[:]->time    
    else
     f1 := addfile(fils1, "r")
      time:=f1->time    
    end if

    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    
   if dimsizes(fils1).ne.1 then 
      va = lonFlip(f1[:]->va(timeind,{1000},:,:))
   else
      va = lonFlip(f1->va(timeind,{1000},:,:))
   end if
    va = (/dtrend_msg_n(ispan(1,dimsizes(va&time),1),va,False,True,0)/)
  
    va_annual = month_to_season(va, "MAM")

    va_clm = dim_avg_n_Wrap(va_annual, 0)

    diff = va_clm

    printVarSummary(va_clm)
    printVarSummary(vwnd_clm)
    
    diff = (/va_clm - vwnd_clm/) 

    ; ;;;0 -20 100-150
    ; siindex = avg(diff({0:-20},{100:150}))
    ; ;;

    ;;;0 -20 150-
    printVarSummary(diff)
  
    siindex =  avg(diff({-10:10},{210:270}))
    ;;;0 -20 150-
    ; siindex =  avg(diff({0:-10},{210:280})) - avg(diff({0:10},{210:280}))
    ; siindex =  avg(diff({0:-10},{210:280})) 
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
    ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{210:270}), (/1,2/)) ;;nino3
     ; nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-10:0},{270:280}), (/1,2/)) ;;nino1-2
    nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{250:280}), (/1,2/)) ;;nino3

    nino3_4_std = dim_stddev_n_Wrap(nino3_4, 0)

    print(nino3_4_std)


    yearind := ind(nino3_4.ge.0.8*nino3_4_std)


    tos_nino := nino3_4(yearind)


    nino_intensity = avg(tos_nino)

    enino_intensity(k) = (/nino_intensity - obs_nino_intensity/)
    ; enino_intensity(k) = nino_intensity
    ;;;
end do 
    ; asciiwrite("diff-Nino(250-280)_enino_intensity_30model", enino_intensity)
    ; asciiwrite("enino_intensity_30model", enino_intensity)
    asciiwrite("vep_north_eastern_pacific_bais_30model", dbitcz)
    ; exit()
end 