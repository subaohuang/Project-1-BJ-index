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
msf_total=new((/dim,12,144,84/), double)


yearstart=200601
yearend=201212  
year1=2006
year2=2012

year=year2-year1+1
locationxi_rea=new((/dim,year/), float)
;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")

 diri="/mnt/f/share/CMIP5/msf_rcp45/"
 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_rcp45_r1i1p1_200601-201212.nc")


;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    time=f1->time(timeind)
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind=get1Dindex(f1->level,z0)
   
    msf_time_r= f1->msf_time(levind,:,timeind)
   
    delete(f1)
    msf_total(i,:,:,:)=msf_time_r
    printVarSummary(msf_total)
    delete(msf_time_r)
    delete(timeind)
    delete(date)
    end do
    msf_amip=dim_avg_n_Wrap(msf_total, 0)
    fi_time_new=msf_amip
    delete(msf_amip)

    nlev=dimsizes(fi_time_new&level)
    ntime=dimsizes(fi_time_new&time)
    nlon=dimsizes(fi_time_new&lon)
 
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diri="/mnt/f/share/CMIP5/msf_rcp45/"
    filo1="msf_time_rcp45_global.nc"
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