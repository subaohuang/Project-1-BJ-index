
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
    yearstart = 190001
    yearend = 200412
    hcindex_all = new(dim_model,float)
    std_diff_all = new(dim_model,float)
    dsst_diff_all = new(dim_model,float)
;;;

   
do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart1 = 190001
    yearend1 = 194912

    yearstart2 = 196001
    yearend2 = 200412



    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP5-derived/msf/1850-1949/msf_hc_"+model+"_history_r1i1p1_185001-194912.nc",\
               "/home/yangsong3/data-observation/linshh/CMIP5-derived/msf/1950-2004/msf_hc_"+model+"_history_r1i1p1_195001-200412.nc"/)
    f1 := addfiles(fils1, "r")
    time:=f1[:]->time

    date:=cd_calendar(time, 1)
    delete(time)
    timeind1:=ind(date.le.yearend1.and.date.ge.yearstart1)
    timeind2:=ind(date.le.yearend2.and.date.ge.yearstart2)
    delete(date)     
    

    msf1= f1[:]->msf(timeind1,:,:,:)
    msf2= f1[:]->msf(timeind2,:,:,:)

    msf1 = (/msf1/10^9/)
    msf2 = (/msf2/10^9/)


    msf_diff = dim_avg_n_Wrap(msf1, 0)
    msf_diff = dim_avg_n_Wrap(msf2, 0) - dim_avg_n_Wrap(msf1, 0)
 
 

    printVarSummary(msf_diff)
  
    ; hcindex =  avg(dim_avg_n(msf_diff({1000:150},{-10:10},{210:270}),(/1,2/))) ;;;HC 强度.nino3 区
    hcindex =  avg(dim_avg_n(msf_diff({1000:150},{-10:10},{240:270}),(/1,2/))) ;;;HC 强度.nino3 区


    hcindex_all(k) = hcindex

    delete(hcindex)



    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind1:=ind(date.le.yearend1.and.date.ge.yearstart1)
    timeind2:=ind(date.le.yearend2.and.date.ge.yearstart2)
 

    delete(date)     
    tos1 = lonFlip(f2[:]->tos(timeind1,:,:))
    tos2 = lonFlip(f2[:]->tos(timeind2,:,:))

    tos_clm1 = dim_avg_n_Wrap(tos1, 0)
    tos_clm2 = dim_avg_n_Wrap(tos2, 0)

    dsst_clm1  = avg(tos_clm1({-5:5},{200:280})) - avg(tos_clm1({-5:5},{80:160}))
    dsst_clm2  = avg(tos_clm2({-5:5},{200:280})) - avg(tos_clm2({-5:5},{80:160}))

    dsst_diff = dsst_clm2 
    dsst_diff = dsst_clm2 - dsst_clm1

    dsst_diff_all(k) = dsst_diff
    delete(dsst_diff)

    tos1 = (/rmMonAnnCycTLL(tos1)/)
    tos2 = (/rmMonAnnCycTLL(tos2)/)


    ;;;计算方差
     std1  = dim_stddev_n_Wrap(tos1, 0)
     std2  = dim_stddev_n_Wrap(tos2, 0)
     
     ; std1_avg = avg(std1({-5:5},{210:270}))  ;;;;;nino3
     ; std2_avg = avg(std2({-5:5},{210:270}))  ;;;;;nino3

     std1_avg = avg(std1({-5:5},{240:270}))  ;;;;;nino3
     std2_avg = avg(std2({-5:5},{240:270}))  ;;;;;nino3
 
    std_diff = dim_avg_n_Wrap(std1_avg, 0)
    std_diff = dim_avg_n_Wrap(std2_avg, 0) - dim_avg_n_Wrap(std1_avg, 0)

    std_diff_all(k) = std_diff
    delete(std_diff)
end do 
    ; asciiwrite("diff-Nino(250-280)_enino_intensity_30model", enino_intensity)
    asciiwrite("HC_change_30model_240-270", hcindex_all)
    asciiwrite("sta_std_change_30model_240-270",std_diff_all)
    asciiwrite("dsst_change_30model",dsst_diff_all)
    ; exit()
end 