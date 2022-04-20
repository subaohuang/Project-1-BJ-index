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

    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
    p_clm = clmMonTLL(p)

    printVarSummary(p_clm)
    printMinMax(p_clm, 1)
    dim = dimsizes(p_clm)
;;;

;;;读入ERA的数据
    f1=addfile("/home/yangsong3/data-observation/ERA-interim/monthly/uwnd_2.5.mon.mean.nc", "r")
    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.200412)
    uwnd = lonFlip(short2flt(f1->u(timeind,{1000},{-30:30},:)))
    printVarSummary(uwnd)

    uwnd= (/dtrend_msg_n(ispan(1,dimsizes(uwnd&time),1),uwnd,False,True,0)/)
    uwnd_clm = clmMonTLL(uwnd)
    printVarSummary(uwnd_clm)
    printMinMax(uwnd_clm, 1)

    f2=addfile("/home/yangsong3/data-observation/ERA-interim/monthly/vwnd_2.5.mon.mean.nc", "r")
    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.200412)
    vwnd = lonFlip(short2flt(f2->v(timeind,{1000},{-30:30},:)))
    printVarSummary(vwnd)

    vwnd= (/dtrend_msg_n(ispan(1,dimsizes(vwnd&time),1),vwnd,False,True,0)/)
    vwnd_clm = clmMonTLL(vwnd)
    printVarSummary(vwnd_clm)
    printMinMax(vwnd_clm, 1)
;;;
;;;



;;;
    yearstart = 195001
    yearend = 200412
    pre_bias = new((/dim_model,12,dim(1),dim(2)/),float)
    u_bias = new((/dim_model,12,dim(1),dim(2)/),float)
    v_bias = new((/dim_model,12,dim(1),dim(2)/),float)
  
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

    PR_clm = clmMonTLL(pr)

    pre_diff = PR_clm
    pre_diff = (/PR_clm - p_clm/) 
    printVarSummary(pre_bias)
    printVarSummary(pre_diff)
   
    pre_bias(k,:,:,:) = pre_diff


  ;;;;read  ua
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/u_chazhi/ua_Amon_"+model+"_historical_r1i1p1_"+"*")
    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    ; timeind:=ind(date.le.yearend.and.date.ge.yearstart)s
    timeind:=ind(date.le.200412.and.date.ge.195001)

    delete(date)     
    ua = lonFlip(f2[:]->ua(timeind,0,{-30:30},:))
    printVarSummary(ua)

   
    ua = (/dtrend_msg_n(ispan(1,dimsizes(ua&time),1),ua,False,True,0)/)
    ua = (/ua*86400/) ;;conver to mm/day

    ua_clm = clmMonTLL(ua)

    ua_diff = ua_clm
 
    ua_diff = (/ua_clm - uwnd_clm/) 
    u_bias(k,:,:,:) = ua_diff

   ;;;;read va

    fils3 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/v_chazhi/va_Amon_"+model+"_historical_r1i1p1_"+"*")
    f3 := addfiles(fils3, "r")

    time:=f3[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    ; timeind:=ind(date.le.yearend.and.date.ge.yearstart)s
    timeind:=ind(date.le.200412.and.date.ge.195001)

    delete(date)     
    va = lonFlip(f3[:]->va(timeind,0,{-30:30},:))
    printVarSummary(va)

    va = (/dtrend_msg_n(ispan(1,dimsizes(va&time),1),va,False,True,0)/)
    va = (/va*86400/) ;;conver to mm/day

    va_clm = clmMonTLL(va)

    va_diff = va_clm
    va_diff = (/va_clm - uwnd_clm/) 
    v_bias(k,:,:,:) = va_diff
end do 
    pre_bias!0 = "model"
    u_bias! 0 = "model"
    v_bias! 0 = "model"
    
    printVarSummary(pre_bias)
    printVarSummary(u_bias)


    system("/bin/rm -f all_month_pre_and_1000uv_bias.nc")
    out = addfile("all_month_pre_and_1000uv_bias.nc", "c")
    out->pre_bias = pre_bias
    out->ua_bias = u_bias
    out->va_bias = v_bias



end 