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

    yearstart = 197901
    yearend = 200412
    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)

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
    ; siindex =  avg(diff({0:10},{210:280}))
    ;;;0 -20 150-
    siindex =  avg(diff({0:-10},{210:280})) - avg(diff({0:10},{210:280}))
    ;;
    dbitcz(k) = siindex 

    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    tos = lonFlip(f2[:]->tos(timeind,:,:))

   
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos = (/rmMonAnnCycTLL(tos)/)
    tos_djf = month_to_season(tos, "DJF")


    nino3_4 = dim_avg_n_Wrap(tos_djf(:,{-5:5},{190:240}), (/1,2/))

    nino3_4_std = dim_stddev_n_Wrap(nino3_4, 0)

    yearind := ind(nino3_4.ge.nino3_4_std)

    tos_nino := nino3_4(yearind)

    nino_intensity = avg(tos_nino)

    enino_intensity(k) = nino_intensity

    ;;;
end do 

    asciiwrite("enino_intensity_30model", enino_intensity)
     asciiwrite("north_eastern_pacific_baise_gradient_30model", dbitcz)
  exit()
end 