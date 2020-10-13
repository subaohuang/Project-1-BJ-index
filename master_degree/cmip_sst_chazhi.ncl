
;********************************************************
function  sst_cold_annual(sst,timeind,year1,year2)
begin 

    ;;;先得到秋季和冬季的SONDJF海温
    sst_son=month_to_season(sst(time|:,lat|:,lon|:), "SON")
    sst_djf=month_to_season(sst(time|:,lat|:,lon|:), "DJF")
    printVarSummary(sst_son)

    sst_son&time=ispan(year1,year2,1)
    sst_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195lonmin1,195lonmin2的平均
    ;;注意点
    sst_djf(0:dimsizes(sst_djf&time)-2,:,:)=(/sst_djf(1:,:,:)/)
    sst_djf(dimsizes(sst_djf&time)-1,:,:)=(/sst(dimsizes(timeind)-1,:,:)/)
    
    ;;;
    sst_mean_annual=new((/dimsizes(timeind)/12,dimsizes(sst&lat),dimsizes(sst&lon)/),typeof(sst))
    printVarSummary(sst_mean_annual)
    ; sst_mean_annual=month_to_annual(sst_mean, 1)
     ; sst_mean_annual=(/sst_son(:,:,:)+sst_djf(:,:,:)/)*0.5
     sst_mean_annual=(/sst_djf(:,:,:)/)*1.0   ;取冬季

    delete(sst_son)
    delete(sst_djf)

    sst_mean_annual!0="time"
    sst_mean_annual!1="lat"
    sst_mean_annual!2="lon"
    sst_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    sst_mean_annual&lon=sst&lon
    sst_mean_annual&lat=sst&lat
    printVarSummary(sst_mean_annual)
return sst_mean_annual
end

begin

; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

models=(/"CanESM2", "CESM1-BGC", "CESM1-CAM5","CMCC-CMS","FIO-ESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H",\
    "GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC",\
    "HadGEM2-ES","MIROC-ESM-CHEM","NorESM1-ME"/)
;;;
printVarSummary(models)

yearstart=197901
yearend=201212  
year1=1979
year2=2012
year=year2-year1+1
dim=dimsizes(models)

rc_sst_total=new(dim, float)
rc_sst_yield=new(dim, float)

sst_rea=new((/dim,year/), float)

p  = 0.975
;;;循环
do i=0,dimsizes(models)-1

    yearstart=197901
    yearend=201212  
    year1=1979
    year2=2012
    year=year2-year1+1

    xi_model=models(i)
    diri="/mnt/f/share/CMIP5/sst_cha/"

    ; if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC" then
    ; yearstart=197901
    ; yearend=200412
    ; year1=1979
    ; year2=2004
    ; end if
   print("******************"+xi_model+"*********************")
    fils1=systemfunc("ls /mnt/i/CMIP5/history/tos_new/"+xi_model+"/r1i1p1/tos_Omon_"+xi_model+"_historical_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/rcp45/tos_new/"+xi_model+"/r1i1p1/tos_Omon_"+xi_model+"_rcp45_r1i1p1_"+"*"+".nc")

    printVarSummary(fils1)
    printVarSummary(fils2)
    print(fils1)
    print(fils2)
    
    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2

    f1=addfiles(fils, "r")
    delete(fils)
    vNames= getfilevarnames(f1[0]) ; get variable names of all groups on file
    print(vNames)

    
    if (any(vNames.eq."LONN179_180")) then 
    time=f1[:]->TIME
    date=cd_calendar(time, 1)

    print("******************日期*********************")

    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    if xi_model.eq."FGOALS-s2" then
    sst_r=f1[:]->TOS(timeind,:,:)
    else
    sst_r=lonFlip(f1[:]->TOS(timeind,:,:))
    end if    
    
     if xi_model.eq."MIROC5" then
      sst_r=where(sst_r.eq.0, sst_r@_FillValue, sst_r)
    end if
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
    
    ;;;插值 sst 

    xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)

     print("*****************开始水平插值************************")  
     sst=linint2_Wrap(sst_r&lon, sst_r&lat, sst_r, \
                      True, xo, yo, 0)


     sst!1="lat"
     sst!2="lon"
     
     sst&lat@units="degrees_north"
     sst&lon@units="degrees_east"

     delete(sst_r)
     printVarSummary(sst)

    ntime=dimsizes(sst&time)
    nlon=dimsizes(sst&lon)
    nlat=dimsizes(sst&lat)


;***************************************************
;写入nc文件
;*************************************************
    diro="/mnt/f/share/CMIP5/tos_chazhi"
    filo1="tos_Omon_"+xi_model+"_chazhi_r1i1p1_197901-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "tos_new 插值的结果"    
    fAtt@creation_date = "2018.12.7"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
   
    filevardef(fout, "time", typeof(sst&time),(/"time"/))
    filevardef(fout, "lat", typeof(sst&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(sst&lon),(/"lon"/))
    filevardef(fout, "sst", typeof(sst),(/"time","lat","lon"/))
 
     filevarattdef(fout, "time", sst&time)
     filevarattdef(fout, "lat", sst&lat)
     filevarattdef(fout, "lon", sst&lon)
    setfileoption(fout,"DefineMode",False)

    fout->time=(/sst&time/)
    fout->lat=(/sst&lat/)
    fout->lon=(/sst&lon/)
    fout->sst=(/sst/)
   delete(sst)
   delete(filo1)
   delete(ntime)
   delete(nlat)
   delete(nlon) 
end do 
end 