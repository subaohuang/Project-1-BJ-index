
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



    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP5-derived/msf_wc/1850-1949/msf_wc_"+model+"_history_r1i1p1_185001-194912.nc",\
               "/home/yangsong3/data-observation/linshh/CMIP5-derived/msf_wc/1950-2004/msf_wc_"+model+"_history_r1i1p1_195001-200412.nc"/)
    f1 := addfiles(fils1, "r")
    time:=f1[:]->time

    date:=cd_calendar(time, 1)
    delete(time)
    timeind1:=ind(date.le.yearend1.and.date.ge.yearstart1)
    timeind2:=ind(date.le.yearend2.and.date.ge.yearstart2)
    delete(date)     
    

    msf1= dble2flt(f1[:]->msf(timeind1,:,:,:))
    msf2= dble2flt(f1[:]->msf(timeind2,:,:,:))

    msf1 = (/tofloat(msf1/10^9)/)
    msf2 = (/tofloat(msf2/10^9)/)



    msf_diff = dim_avg_n_Wrap(msf1, 0)
    msf_diff = dim_avg_n_Wrap(msf2, 0) - dim_avg_n_Wrap(msf1, 0)
 
 

    printVarSummary(msf_diff)
  
    ; hcindex =  avg(dim_avg_n(msf_diff({1000:150},{-10:10},{210:270}),(/1,2/))) ;;;HC 强度.nino3 区
    hcindex =  avg(dim_avg_n(msf_diff({1000:150},{-5:5},{150:210}),(/1,2/))) ;;;HC 强度.nino3 区


    hcindex_all(k) = hcindex

    delete(hcindex)

end do 

    asciiwrite("WC_change_30model_150-210", hcindex_all)
 
end 