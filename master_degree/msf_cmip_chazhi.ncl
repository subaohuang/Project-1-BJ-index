;;;
begin
;画时间序列图的时候，把BCC-CSM1-1给去掉
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models=(/"BNU-ESM"/)
dim=dimsizes(models)
msf_total=new((/dim,12,144,408/), double)


yearstart=197901
yearend=200812  
year1=1979
year2=2008

year=year2-year1+1
locationxi_rea=new((/dim,year/), float)
;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")

 diri1="/mnt/f/share/CMIP5/msf_history/"
 diri2="/mnt/f/share/CMIP5/msf_rcp45/"
 fils1=systemfunc("ls "+diri1+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")
 fils2=systemfunc("ls "+diri2+"msf_time_"+model+"_rcp45_r1i1p1_200601-201212.nc")
 

    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    yearstart1=197901
    yearend1=200512

    yearstart2=200401   ;;;因为存在200512开始的
    yearend2=201212

    year1=1979
    year2=2012

;;;********************************************* 
    f  =addfiles(fils, "r")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f[:]->time, 1)

    date1=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)


    timeind=ind(date.ge.197901.and.date.le.201212)
    time=f[:]->time(timeind)

    timeind1=ind(date1.le.yearend1.and.date1.ge.yearstart1)
    timeind2=ind(date2.le.yearend2.and.date2.ge.yearstart2)
    
    delete([/date1 ,date2/])

    year=year2-year1+1

    msf_time1=f1->msf_time(:,:,timeind1)
    delete(f1)
    printVarSummary(msf_time1)
    dim_data=dimsizes(msf_time1)
    msf_time_r=new((/dim_data(0),dim_data(1),year*12/),typeof(msf_time1))
    msf_time2=f2->msf_time(:,:,timeind2)
    delete(f2)

    printVarSummary(msf_time2)

    msf_time_r(:,:,0:dimsizes(timeind1)-1)=msf_time1
    msf_time_r(:,:,dimsizes(timeind1):)=msf_time2

    msf_total(i,:,:,:)=msf_time_r
    printVarSummary(msf_total)
    delete(msf_time_r)
    delete(timeind)
    delete(date)
    delete(timeind1)
    delete(timeind2)
    end do

    msf_cmip=dim_avg_n_Wrap(msf_total, 0)
    fi_time_new=msf_cmip
    delete(msf_cmip)

    nlev=dimsizes(fi_time_new&level)
    ntime=dimsizes(fi_time_new&time)
    nlon=dimsizes(fi_time_new&lon)
 
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diri="/mnt/f/share/CMIP5/msf_history/"
    filo1="msf_time_cmip_global.nc"
    system("/bin/rm -f "+diri + filo1+"")
    fout=addfile(diri + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量1950-200512，取5到-5平均的结果"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(fi_time_new&level), (/"level"/))
    filevardef(fout, "lon", typeof(fi_time_new&lon),(/"lon"/))
    filevardef(fout, "time", typeof(fi_time_new&time),(/"time"/))

    filevardef(fout, "msf_time", typeof(fi_time_new),(/"level","lon","time"/))

    filevarattdef(fout, "lon", fi_time_new&lon)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "msf_time", fi_time_new)
    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_time_new&level/)
    fout->lon=(/fi_time_new&lon/)
    fout->time=(/fi_time_new&time/)
    fout->msf_time=(/fi_time_new/)

    delete(time)
    delete(fi_time_new)
end