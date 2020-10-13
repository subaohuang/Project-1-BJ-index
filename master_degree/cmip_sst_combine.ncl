
;;;********************************************************************
;;;主程序
;;;******************************************************************
begin

; models=(/"CSIRO-Mk3-6-0","INMCM4","GFDL-CM3","IPSL-CM5A-MR"/)
; dong_model=(/"MIROC5","MPI-ESM-MR"/)

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

;;;全球
lat_sst=fspan(89.5, -89.5, 180)
lon_sst=fspan(0.5, 359.5, 360)
;;;

sst_ximodel=new((/dimsizes(models),408,dimsizes(lat_sst),dimsizes(lon_sst)/),float)
dim=dimsizes(models)

diri_output="/mnt/f/share/CMIP5/tos_chazhi/"

 do j=0,dim-1
    
    model=models(j)
    yearstart=197901
    yearend=201212
    
    year1=1979
    year2=2012
    
    print("******************"+model+"*********************")
    fils1=systemfunc("ls /mnt/i/CMIP5/history/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_historical_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/rcp45/tos_new/"+model+"/r1i1p1/tos_Omon_"+model+"_rcp45_r1i1p1_"+"*"+".nc")

    printVarSummary(fils1)
    printVarSummary(fils2)

    printVarSummary(fils1)
    printVarSummary(fils2)
    print(fils1)
    print(fils2)
    
    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "join")
    delete(fils1)
    delete(fils2)

    vNames= getfilevarnames(f1[0]) ; get variable names of all groups on file
    print(vNames)

    if (any(vNames.eq."LONN179_180")) then 
    time=f1[:]->TIME
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    sst_r=lonFlip(f1[:]->TOS(timeind,:,:))
    sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r)
    dNames = getvardims(f1[:]->TOS)
      delete(f1)
    print(dNames)
    sst_r!0="time"
    sst_r!1="lat"
    sst_r!2="lon"
    end if
    delete(vNames)
    delete(dNames)
    sst_ximodel(j,:,:,:)=sst_r
    delete([/sst_r/])
    end do
    
    sst_xi_mean=dim_avg_n(sst_ximodel,0)
    printVarSummary(sst_ximodel)
    copy_VarMeta(sst_ximodel(1,:,:,:), sst_xi_mean)

   sst_xi_mean!0="time"
   sst_xi_mean!1="lat"
   sst_xi_mean!2="lon"
   ntime=dimsizes(sst_xi_mean&time)
   nlat=dimsizes(sst_xi_mean&lat)
   nlon=dimsizes(sst_xi_mean&lon)

   
   printVarSummary(sst_xi_mean)

    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="sst_cmip_r1i1p1_197901-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "sst_xi_mean插值后的结果，取1000hpa"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(sst_xi_mean&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(sst_xi_mean&lon),(/"lon"/))
    filevardef(fout, "time", typeof(sst_xi_mean&time),(/"time"/))
    filevardef(fout, "sst", typeof(sst_xi_mean),(/"time","lat","lon"/))



    filevarattdef(fout, "lat", sst_xi_mean&lat)
    filevarattdef(fout, "lon", sst_xi_mean&lon)
    filevarattdef(fout, "time", sst_xi_mean&time)
    
    filevarattdef(fout, "sst", sst_xi_mean)

    
    setfileoption(fout,"DefineMode",False)

    fout->lat=(/sst_xi_mean&lat/)
    fout->lon=(/sst_xi_mean&lon/)
    fout->time=(/sst_xi_mean&time/)
    fout->sst=(/sst_xi_mean/)

    delete(sst_xi_mean)
    delete(ntime)
    delete(nlon)
    delete(nlat)
    end 



    