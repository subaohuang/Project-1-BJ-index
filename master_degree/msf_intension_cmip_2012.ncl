
function level_mean(data)
local data,msf_mean
begin
plev=data&level
psfc=conform_dims((/dimsizes(data&time),dimsizes(data&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
msfdp=data(level|:,time|:,lon|:)*dp
copy_VarCoords(data(level|:,time|:,lon|:), msfdp)
msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(data(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end 
;********************************************************
;draw_msf_hov
begin

;;有海温的几个模式的width情况
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FGOALS-g2","FGOALS-s2","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; ;;FGOALS-g2和s2没有海温数据l;;同时添加了models=(/"CanESM2",GISS-E2-R","MIROC-ESM-CHEM","NorESM1-ME"/)
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CanESM2","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-R",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
;;共33个模式
; ;FGOALS-g2和s2没有海温数据l;;同时添加了models=(/"CanESM2",GISS-E2-R","MIROC-ESM-CHEM","NorESM1-ME"/)
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CCSM4","CanESM2","CESM1-BGC","CESM1-CAM5","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M",\
;     "BNU-ESM","CanESM2", "CESM1-BGC", "CESM1-CAM5","CMCC-CMS","FIO-ESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H",\
;     "GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC",\
;     "HadGEM2-ES","MIROC-ESM-CHEM","NorESM1-ME"/)
; models=(/"cmip","BNU-ESM","CanESM2", "CESM1-BGC", "CESM1-CAM5","CMCC-CMS","FIO-ESM","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H",\
;     "GISS-E2-H-CC","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC",\
;     "HadGEM2-ES","MIROC-ESM-CHEM","NorESM1-ME"/)
; ;;;

;;;具有AMIP实验和HIS实验的两类模式

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
        "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CMCC-CM","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
; ;;


; models=(/"GFDL-ESM2M"/)
dim=dimsizes(models)

rc_inten_total =new(dim, float)
rc_inten_yield = new(dim, float)

rc_inten_id_total=rc_inten_total
rc_inten_id_yield=rc_inten_yield

rc_inten_alt_total=rc_inten_total
rc_inten_alt_yield=rc_inten_yield

rc_inten_yield_90=rc_inten_yield
rc_inten_id_yield_90 = rc_inten_yield
rc_inten_alt_yield_90 = rc_inten_yield


yearstart=197901
yearend=200812  
year1=1979
year2=2008
; p  = 0.975
p=0.95

year=year2-year1+1
inten=new((/dim,year/), float)
inten_alt=new((/dim,year/), float)
inten_id=new((/dim,year/), float)

;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("**********model*************="+model)
 ; model="ACCESS1-0"
 diri="/mnt/f/share/CMIP5/msf_history/"
 diri1="/mnt/f/share/CMIP5/msf_rcp45/"

 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")
 fils2=systemfunc("ls "+diri1+"msf_time_"+model+"_rcp45_r1i1p1_200601-201212.nc")
  
  fils=new(2,string)
  fils(0)=fils1
  fils(1)=fils2

yearstart1=197901
yearend1=200512

yearstart2=200401   ;;;因为存在200512开始的
yearend2=200812

year1=1979
year2=2008

;;;********************************************* 
    f=addfiles(fils, "r")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f[:]->time, 1)

    date1=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)


    timeind=ind(date.ge.197901.and.date.le.200812)
    timeind1=ind(date1.le.yearend1.and.date1.ge.yearstart1)
    timeind2=ind(date2.le.yearend2.and.date2.ge.yearstart2)
    
    delete([/date1 ,date2/])
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,150.0,100.0/)
    levind1=get1Dindex(f1->level,z0)
 
    levind2=get1Dindex(f2->level,z0)

    year=year2-year1+1
 
    ; msf_time=f1->msf_time({850:500},:,timeind)
    msf_time1=f1->msf_time(levind1,:,timeind1)
    delete(f1)
    delete(levind1)
   
    printVarSummary(msf_time1)
    dim_data=dimsizes(msf_time1)
    msf_time=new((/dim_data(0),dim_data(1),year*12/),typeof(msf_time1))


    msf_time2=f2->msf_time(levind2,:,timeind2)
    delete(f2)
    delete(levind2)

    printVarSummary(msf_time2)

    msf_time(:,:,0:dimsizes(timeind1)-1)=msf_time1
    msf_time(:,:,dimsizes(timeind1):)=msf_time2



    delete(msf_time1)
    delete(msf_time2)
    delete(timeind1)
    delete(timeind2)
    printVarSummary(msf_time)

    delete(f)
    msf_time=(/msf_time/10^9/)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    print(msf_time&time)

    delete(date)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)
 
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ;;;
  
    ; msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
     ; msf_mean_annual=month_to_annual(msf_mean, opt)
     msf_annual=msf_djf
    msf_annual=(/month_to_annual(msf_mean(time|:,level|:,lon|:), 1)/)

    delete(msf_mean)
    printVarSummary(msf_annual)

    msf_mean_annual=level_mean(msf_annual(:,:,:))  ;time,lon
    ; msf_mean_annual=level_mean(msf_annual(:,{700:300},:))  ;time,lon
    delete(msf_annual)
    delete(msf_son)
    delete(msf_djf)
   
    printVarSummary(msf_mean_annual)
 

    msf_inten_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{150:210}), 1)  ;;150:210
      ; msf_inten_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{120:180}), 1)  ;;150:210
    msf_inten_id_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{70:100}), 1)
    msf_inten_alt_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{320:340}), 1)

    delete(msf_mean_annual)
    ;;;calculate trend

    msf_inten_rc= regCoef_n(ispan(1, dimsizes(msf_inten_aunnal&time), 1), msf_inten_aunnal,0,0)
    printVarSummary(msf_inten_rc)

    msf_inten_id_rc= regCoef_n(ispan(1, dimsizes(msf_inten_aunnal&time), 1), msf_inten_id_aunnal,0,0)
    printVarSummary(msf_inten_rc)
    msf_inten_alt_rc= regCoef_n(ispan(1, dimsizes(msf_inten_aunnal&time), 1), msf_inten_alt_aunnal,0,0)
    printVarSummary(msf_inten_rc)

    df   = msf_inten_rc@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
   
    rc_inten_total(i)=(/tofloat(msf_inten_rc)/)*10.0
    rc_inten_yield(i)=tofloat(msf_inten_rc@rstd)*10.0
    rc_inten_yield_90(i)=tofloat(msf_inten_rc@rstd)*T*10.0
    
    rc_inten_id_total(i)=(/tofloat(msf_inten_id_rc)/)*10.0
     rc_inten_id_yield(i)=tofloat(msf_inten_id_rc@rstd)*10.0
    rc_inten_id_yield_90(i)=tofloat(msf_inten_id_rc@rstd*T)*10.0

    rc_inten_alt_total(i)=(/tofloat(msf_inten_alt_rc)/)*10.0
    rc_inten_alt_yield(i)=tofloat(msf_inten_alt_rc@rstd)*10.0
    rc_inten_alt_yield_90(i)=tofloat(msf_inten_alt_rc@rstd*T)*10.0
    ;;;calculate trend
    inten(i,:)=tofloat(msf_inten_aunnal)
    inten_id(i,:)=tofloat(msf_inten_id_aunnal)
    inten_alt(i,:)=tofloat(msf_inten_alt_aunnal)

    delete([/msf_inten_id_aunnal,msf_inten_alt_aunnal,msf_inten_aunnal/])
    delete([/msf_inten_id_rc,msf_inten_alt_rc,msf_inten_rc/])
    delete(year)
    end do


inten!0="model"
inten!1="time"
inten&model=models
inten&time=ispan(year1,year2,1)


inten_id!0="model"
inten_id!1="time"
inten_id&model=models
inten_id&time=ispan(year1,year2,1)


inten_alt!0="model"
inten_alt!1="time"
inten_alt&model=models
inten_alt&time=ispan(year1,year2,1)


nmodel=dimsizes(inten&model)
ntime=dimsizes(inten&time)



opt=diri+"cmip_msf_inten_trend.ascii"    
asciiwrite(opt, sprintf ("%6.3f,",rc_inten_total))

opt1=diri+"cmip_msf_inten_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_yield))

opt=diri+"cmip_msf_inten_id_trend.ascii"    
asciiwrite(opt, sprintf ("%6.3f,",rc_inten_id_total))

opt1=diri+"cmip_msf_inten_id_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_id_yield))

opt=diri+"cmip_msf_inten_alt_trend.ascii"    
asciiwrite(opt, sprintf ("%6.3f,",rc_inten_alt_total))

opt1=diri+"cmip_msf_inten_alt_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_alt_yield))


opt1=diri+"cmip_msf_inten_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_yield_90))

opt1=diri+"cmip_msf_inten_id_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_id_yield_90))

opt1=diri+"cmip_msf_inten_alt_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_alt_yield_90))


;***************************************************
;写入nc文件
;*************************************************
    diro=diri
    filo1="intension_cmip_r1i1p1_djf.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "intensin,msf djf平均后计算得到的位置结果"    
    fAtt@creation_date = "2018.4.26"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"model","time"/)
    dimsize=(/nmodel,ntime/)
    dimunlim=(/False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "model","integer",(/"model"/))
    filevardef(fout, "time", typeof(inten&time),(/"time"/))

    filevardef(fout, "intension_pa", typeof(inten),(/"model","time"/))
    filevardef(fout, "intension_id", typeof(inten_id),(/"model","time"/))
      filevardef(fout, "intension_alt", typeof(inten_alt),(/"model","time"/))


    filevarattdef(fout, "intension_pa", inten)
    filevarattdef(fout, "intension_id", inten_id)
    filevarattdef(fout, "intension_alt", inten_alt)


    
    setfileoption(fout,"DefineMode",False)

    fout->model=(/ispan(1,dim,1)/)
    fout->time=(/inten&time/)

    fout->intension_pa=(/inten/)
    fout->intension_id=(/inten_id/)
    fout->intension_alt=(/inten_alt/)

end
