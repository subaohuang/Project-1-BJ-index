
;;**************************************************************
;;;main program
;;***************************************************************
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

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;     "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;     "GFDL-CM4", "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         "HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM"/)
model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-CSM2-MR",\
              "CAMS-CSM1-0","CanESM5", "CAS-ESM2-0" ,"CESM2","CESM2-WACCM","CIESM","CMCC-CM2-SR5","CMCC-ESM2","CNRM-CM6-1"\
   , "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
     "GISS-E2-1-G","HadGEM3-GC31-LL","IPSL-CM6A-LR","MIROC-ES2L","MIROC6", \
    "MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3",  "NorESM2-LM"/)

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5"/)

; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(dim_model,graphic)

std_all = new((/dim_model,24,144/),float)

do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start read"+" "+model+"**********************************************")
    yearstart = 195001
    yearend = 199912

    yearstart2 = 205001
    yearend2 = 209912

    fils1 := systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/historical/tos_Omon_"+"*")
    fils2 := systemfunc("ls /home/yangsong3/CMIP6/hank_2.5x2.5/"+model+"/ssp585/tos_Omon_"+"*")
    ; print(fils1)

   if dimsizes(fils1).eq.1 then 
        f1 := addfile(fils1, "r")

        time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)

        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        delete(date)     
        sst_his := lonFlip(f1->tos(timeind,{-30:30},:))
        delete(timeind)
    else 
        f1 := addfiles(fils1, "r")
        time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)
        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        delete(date)     
        sst_his := lonFlip(f1[:]->tos(timeind,{-30:30},:))
        delete(timeind)
    end if


    if dimsizes(fils2).eq.1 then 

        f2 := addfile(fils2, "r")

        time:=f2->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)

        timeind:=ind(date.le.yearend2.and.date.ge.yearstart2)
        delete(date)     
        sst_ssp := lonFlip(f2->tos(timeind,{-30:30},:))
        delete(timeind)
    else
        f2 := addfiles(fils2, "r")
        time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)
        timeind:=ind(date.le.yearend2.and.date.ge.yearstart2)
        delete(date)     
        sst_ssp := lonFlip(f2[:]->tos(timeind,{-30:30},:))
        delete(timeind)
    end if
;;;;    
        ssta_his = sst_his


        ssta_his :=(/rmMonAnnCycTLL(sst_his)/)

        ssta_his = (/dtrend_msg_n(ispan(1,dimsizes(sst_his&time),1),ssta_his,False,True,0)/)
        copy_VarCoords(sst_his, ssta_his)
        ; printVarSummary(ssta_his)

        ssta_ssp = sst_ssp
        ssta_ssp :=(/rmMonAnnCycTLL(sst_ssp)/)

        ssta_ssp = (/dtrend_msg_n(ispan(1,dimsizes(sst_ssp&time),1),ssta_ssp,False,True,0)/)
        copy_VarCoords(sst_ssp, ssta_ssp)
        ; printVarSummary(ssta_ssp)
;;;;

;;;去掉7年以上信号
        ssta_his = (/ ssta_his - runave_n(ssta_his, 7*12,0,0)/)
        ssta_ssp = (/ ssta_ssp - runave_n(ssta_ssp, 7*12,0,0)/)

;;;


;;;;计算标准偏差(5S-5N)

        std1 = dim_stddev_n_Wrap(ssta_his, 0)

        std2 = dim_stddev_n_Wrap(ssta_ssp, 0)
;;;;

        std_diff  = std1

        std_diff = (/std2 - std1/)

        std_all(k,:,:) = std_diff
        delete(std1)
        delete(std2)
end do 

print("************"+"写入nc文件"+"*************")
std_all!0  = "model_name"


systemfunc("rm -f std_diff_SSP-his-2050-2099minus1950-1999.nc")
out = addfile("std_diff_SSP-his-2050-2100minus1950-2000.nc", "c")
out->std_diff = std_all

end
