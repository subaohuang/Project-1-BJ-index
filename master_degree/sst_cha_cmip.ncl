
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

; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FGOALS-g2","FGOALS-s2","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

printVarSummary(models)

yearstart=197901
yearend=200512  
year1=1979
year2=2005
year=year2-year1+1
dim=dimsizes(models)

rc_sst_total=new(dim, float)
rc_sst_yield=new(dim, float)

sst_rea=new((/dim,year/), float)

p  = 0.975
;;;循环
do i=0,dimsizes(models)-1

    yearstart=197901
    yearend=200512  
    year1=1979
    year2=2005
    year=year2-year1+1

    xi_model=models(i)
    diri="/mnt/f/share/CMIP5/sst_cha/"

    if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC" then
    yearstart=197901
    yearend=200412
    year1=1979
    year2=2004
    end if

    fils1=systemfunc("ls /mnt/i/CMIP5/history/tos_new/"+xi_model+"/r1i1p1/tos_Omon_"+xi_model+"_historical_r1i1p1_"+"*"+".nc")
    printVarSummary(fils1)
    print(fils1)
  
    print("******************"+xi_model+"*********************")
    f1=addfiles(fils1, "join")
    delete(fils1)
    vNames= getfilevarnames(f1[:]) ; get variable names of all groups on file
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
    delete(sst_r&time)
    sst_r&time=date(timeind)
    if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC"  then
    sst=sst_cold_annual(sst_r,timeind,1979,2004)
    else 
    sst=sst_cold_annual(sst_r,timeind,1979,2005)
    end if
    
    delete(time)
    delete(date)
    printVarSummary(sst)
    delete(timeind)
;;;计算定义的区域的海温差
    sst=(/month_to_annual(sst_r, 1)/)
    delete(sst_r)
   ; sst_cha= dim_avg_n_Wrap(sst(:,{-5:5},{120:170}),(/1,2/)) - dim_avg_n_Wrap(sst(:,{-5:5},{210:270}),(/1,2/))
    sst_cha= dim_avg_n_Wrap(sst(:,{-5:5},{80:160}),(/1,2/)) - dim_avg_n_Wrap(sst(:,{-5:5},{200:280}),(/1,2/))
    sst_cha!0="time"
    sst_cha&time=ispan(year1,year2,1)
    printVarSummary(sst_cha)

    rc_sst=regCoef_n(ispan(0,dimsizes(sst_cha&time)-1,1),sst_cha , 0, 0)
    rc_sst=(/rc_sst*10/)

    ;;;;
    df   = rc_sst@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
    rc_sst_total(i)=rc_sst
    rc_sst_yield(i)=rc_sst@rstd*T*10.0
    if xi_model.eq."FGOALS-s2".or.xi_model.eq."HadGEM2-CC"  then
    sst_rea(i,0:year-2)=sst_cha
    else 
    sst_rea(i,:)=sst_cha    
    end if
    sst_rea!0="model"
    sst_rea!1="time"
    sst_rea&model=models

    year2=2005;

sst_rea&time=ispan(year1,year2,1)
printVarSummary(sst_rea)
delete(sst_cha)
delete(sst)

end do


nmodel=dimsizes(sst_rea&model)
ntime=dimsizes(sst_rea&time)

opt=diri+"his_sstca_trend.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",rc_sst_total))

opt1=diri+"his_sstca_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_sst_yield))  

;***************************************************
;写入nc文件
;*************************************************
    diro=diri
    filo1="sstca_cmip_r1i1p1_djf.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "location,msf djf平均后计算得到的位置结果"    
    fAtt@creation_date = "2018.4.26"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"model","time"/)
    dimsize=(/nmodel,ntime/)
    dimunlim=(/False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "model","integer",(/"model"/))
    filevardef(fout, "time", typeof(sst_rea&time),(/"time"/))

    filevardef(fout, "sst_cha", typeof(sst_rea),(/"model","time"/))

    
    filevarattdef(fout, "sst_cha", sst_rea)    
    setfileoption(fout,"DefineMode",False)

    fout->model=(/ispan(1,dim,1)/)
    fout->time=(/sst_rea&time/)
    fout->sst_cha=(/sst_rea/)


end