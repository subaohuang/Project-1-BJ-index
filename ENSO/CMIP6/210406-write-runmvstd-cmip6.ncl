
;;**************************************************************
;;;main program
;;***************************************************************
; undef("rmstd") ;;;滑动std
; function rmstd(x,nave)
; begin
; dim = dimsizes(x)
; length=dim(0)
; std = new((/length-nave+1, dimsizes(x&lat), dimsizes(x&lon)/),typeof(x))
; do t=nave,length
;    x(t-nave:t-1,:,:) = (/dtrend_msg_n(ispan(1,nave,1),x(t-nave:t-1,:,:),False,True,0)/)
;    std_mv = dim_stddev_n_Wrap(x(t-nave:t-1,:,:),0)
;    std(t-nave,:,:) = std_mv
; end do 
; return std
; end 

;;;;for each month
undef("rmstd") ;;;滑动std
function rmstd(x,nave)
begin
dim = dimsizes(x)
length=dim(0)/12
std = new((/length-nave+1, dimsizes(x&lat), dimsizes(x&lon)/),typeof(x))

;;;
do t=nave,length

   m = x((t-nave)*12:t*12-1,:,:)
    ;;;rm annual cycle in each running period followed as kim and cai,2013 ncc
   m = (/rmMonAnnCycTLL(m)/)
    ;;; detrend 
   m = (/dtrend_msg_n(ispan(1,nave*12,1),m,False,True,0)/)
    ;;; 7-year running mean 
   m = (/ m - runave_n(m, 7*12,0,0)/)
   std_mv = dim_stddev_n_Wrap(m,0)
   std(t-nave,:,:) = std_mv
   delete(m)
end do 
return std
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
; model_name = (/"CESM2-WACCM","CIESM"/)

; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(dim_model,graphic)

yeartotal = (2099-1850+1)

; std_all = new((/dim_model,yeartotal-50+1,24,144/),float)
std_all = new((/dim_model,yeartotal-50+1,8,68/),float)

do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start read"+" "+model+"**********************************************")
    yearstart = 185001
    yearend = 201412

    yearstart2 = 201501
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
        print(date(ind(date.eq.yearstart)))
        delete(date)     
        sst_his := lonFlip(f1->tos(timeind,{-10:10},:))
        delete(timeind)
    else 
        f1 := addfiles(fils1, "r")
        time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)
        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        print(date(ind(date.eq.yearstart)))
        delete(date)     
        sst_his := lonFlip(f1[:]->tos(timeind,{-10:10},:))
        delete(timeind)
    end if


    if dimsizes(fils2).eq.1 then 

        f2 := addfile(fils2, "r")

        time:=f2->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)

        timeind:=ind(date.le.yearend2.and.date.ge.yearstart2)
        print(date(ind(date.eq.yearend2)))
        delete(date)     
        sst_ssp := lonFlip(f2->tos(timeind,{-10:10},:))
        delete(timeind)
    else
        f2 := addfiles(fils2, "r")
        time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)
        timeind:=ind(date.le.yearend2.and.date.ge.yearstart2)
        print(date(ind(date.eq.yearend2)))
        delete(date)     
        sst_ssp := lonFlip(f2[:]->tos(timeind,{-10:10},:))
        delete(timeind)
    end if
;;;;    

        ; sst_his =(/rmMonAnnCycTLL(sst_his)/)
        ; ssta_his = (/dtrend_msg_n(ispan(1,dimsizes(sst_his&time),1),ssta_his,False,True,0)/)
        ; copy_VarCoords(sst_his, ssta_his)
        ; ; printVarSummary(ssta_his)
        ; sst_ssp =(/rmMonAnnCycTLL(sst_ssp)/)

        ; ssta_ssp = (/dtrend_msg_n(ispan(1,dimsizes(sst_ssp&time),1),ssta_ssp,False,True,0)/)
        ; copy_VarCoords(sst_ssp, ssta_ssp)
        ; ; printVarSummary(ssta_ssp)
        sst =  array_append_record(sst_his,sst_ssp, 0)
        ; sst = (/rmMonAnnCycTLL(sst)/)
   




;;;;计算标准偏差(5S-5N)

        ; ssta_annual = month_to_annual(sst, 1)
        ; std = rmstd(ssta_annual,50)
        std = rmstd(sst(:,:,{120:290}),50)
        std_all(k,:,:,:) = std
  delete(std)

end do 

print("************"+"写入nc文件"+"*************")
std_all!0  = "model_name"


system("/bin/rm -f annual_runmvstd_1900-2099-rmanucycl-eachperiod.nc")
out = addfile("annual_runmvstd_1900-2099-rmanucycl-eachperiod.nc", "c")
out->std_all = std_all

end
