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


; ;;;读入HadiSST海温数据
      diri="/home/yangsong3/data-observation/linshh/data/sst/"
      fils1=systemfunc("ls "+diri+"HadISST_sst_2.5.nc")
      f1=addfile(fils1, "r")
      date:=cd_calendar(f1->time, 1)

      timeind:=ind(date.le.201912.and.date.ge.195001)

      sst =  lonFlip(f1->sst(timeind,{-30:30},:))

      sst = (/dtrend_msg_n(ispan(1,dimsizes(sst&time),1),sst,False,True,0)/)

      sst_clm = clmMonTLL(sst)
      
     dim_sst  = dimsizes(sst_clm)
;;;
    yearstart = 195001
    yearend = 200412
  
    SST_bias = new((/dim_model,dim_sst(0),dim_sst(1),dim_sst(2)/),float)
;;;

   
do k = 0,dim_model-1
;;;;

    model = model_name(k)

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
    
    tos_clm = clmMonTLL(tos)

    tos_diff = tos_clm
    tos_diff = (/tos_clm - sst_clm/)
    SST_bias(k,:,:,:) = tos_diff
 
end do 
    SST_bias!0 = "model"

    system("/bin/rm -f all_month_sst_bias.nc")
    out = addfile("all_month_sst_bias.nc", "c")
    out->SST_bias = SST_bias
   


end 