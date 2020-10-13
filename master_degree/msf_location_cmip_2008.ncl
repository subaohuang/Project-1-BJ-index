
;;;
undef("level_mean")
function level_mean(msf_time)   ;;;;select 850hpa-200hpa 
begin
; msf_time=(/msf_time/10^9/)
plev=msf_time&level
psfc=conform_dims((/dimsizes(msf_time&time),dimsizes(msf_time&lon)/), 1013.0, 0)
ptop=min(plev)
dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
; dp=(/50/)
msfdp=msf_time(level|:,time|:,lon|:)*dp
copy_VarCoords(msf_time(level|:,time|:,lon|:), msfdp)

msf_mean=dim_sum_n_Wrap(msfdp, 0)/dim_sum_n_Wrap(dp, 0)
copy_VarCoords(msf_time(time|:,lon|:,level|1), msf_mean)
return(msf_mean)
end
;;; 
undef("search_location")
function  search_location(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
msf_mean_annual_xi=msf_mean_annual(:,{120:250})
msf_mean_annual_dong=msf_mean_annual(:,{230:290})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+8).gt.0))  then 
         a=(msf_mean_annual_xi(t,i+1)-msf_mean_annual_xi(t,i))/(msf_mean_annual_xi&lon(i+1)-msf_mean_annual_xi&lon(i))
         locationxi(t)=tofloat(msf_mean_annual_xi&lon(i)- msf_mean_annual_xi(t,i)/a)
        ; locationxi(t)=(msf_mean_annual_xi&lon(i)+msf_mean_annual_xi&lon(i+1))/2
        end if    
        if ismissing(locationxi(t)).ne.True  then 
        break
        end if
     end do
end do

;;;judge locationdong
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_dong&lon)-2
        ; if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0).and.(msf_mean_annual_dong(t,i+3).lt.0))  then 
         if  (msf_mean_annual_dong(t,i+1).le.0).and.((msf_mean_annual_dong(t,i).ge.0))  then 
         a=(msf_mean_annual_dong(t,i+1)-msf_mean_annual_dong(t,i))/(msf_mean_annual_dong&lon(i+1)-msf_mean_annual_dong&lon(i))
         locationdong(t)=tofloat(msf_mean_annual_dong&lon(i)- msf_mean_annual_dong(t,i)/a)
        ; locationdong(t)=(msf_mean_annual_dong&lon(i)+msf_mean_annual_dong&lon(i+1))/2
        end if    
        if ismissing(locationdong(t)).ne.True  then 
        break
        end if
     end do
    if ismissing(locationdong(t)).eq.True  then 
   ; co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
   co_in=closest_val(0,msf_mean_annual_dong(t,:))
   locationdong(t)=msf_mean_annual_dong&lon(co_in)
   delete(co_in)
   end if

end do
locationxi=(/linmsg(locationxi, -1)/)
; locationdong=(/linmsg(locationdong, -1)/)
location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 
;********************************************************
;draw_msf_hov
begin

; models=(/"ACCESS1-0","bcc-csm1-1","BNU-ESM","CanESM2" ,"CCSM4","CMCC-CM","CSIRO-Mk3-6-0","FGOALS-s2","FGOALS-g2","GFDL-CM2p1","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;             "GISS-E2-H-CC","GISS-E2-R" ,"GISS-E2-R-CC","HadGEM2-AO","inmcm4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models=(/"ACCESS1-0","BCC-CSM1-1","BNU-ESM","CanESM2" ,"CCSM4","CMCC-CM","CSIRO-Mk3-6-0","FGOALS-s2","FGOALS-g2","GFDL-CM2p1","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;             "GISS-E2-H-CC","GISS-E2-R" ,"HadGEM2-AO","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanESM2" ,"CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","FGOALS-s2","FGOALS-g2","GFDL-CM2p1","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;             "GISS-E2-H-CC","GISS-E2-R" ,"HadGEM2-AO","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)


;;有海温的几个模式的width情况
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
        "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
         "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)

; models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
;         "CSIRO-Mk3-6-0","FIO-ESM","FGOALS-s2","GFDL-CM3","GFDL-ESM2M",\
;          "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-MR","NorESM1-M"/)


printVarSummary(models)


dim=dimsizes(models)
rc_xi_total=new(dim, float)
rc_xi_yield=new(dim, float)

rc_dong_total=new(dim, float)
rc_dong_yield=new(dim, float)

rc_width_total=new(dim, float)
rc_width_yield=new(dim, float)


yearstart=197901
yearend=200812  
year1=1979
year2=2008
year=year2-year1+1

locationxi_rea=new((/dim,year/), float)
locationdong_rea=new((/dim,year/), float)


p  = 0.95
; p=0.95
;;;循环

do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("**********model*************="+model)
 ; model="ACCESS1-0"
 diri="/mnt/f/share/CMIP5/msf_history/"
 diri1="/mnt/f/share/CMIP5/msf_rcp45/"

 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")
 fils2=systemfunc("ls "+diri1+"msf_time_"+model+"_rcp45_r1i1p1_200601-200812.nc")
  
  fils=new(2,string)
  fils(0)=fils1
  fils(1)=fils2

yearstart1=197901
yearend1=200512

yearstart2=200401   ;;;因为存在200512开始的
yearend2=200812

year1=1979
year2=2008

; if model.eq."HadGEM2-CC".or. model.eq."HadGEM2-ES" then
; yearstart=197901
; yearend=200412  
; year1=1979
; year2=2004
; end if


;;;********************************************* 
    f=addfiles(fils, "join")
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
    delete(date)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)
 
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_jja=month_to_season(msf_mean(time|:,level|:,lon|:), "JJA")
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")

    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ;;;
   
   
    msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
    printVarSummary(msf_mean_annual)
    
    ; msf_mean_annual=(/msf_jja(:,{500},:)/)*1.0
    msf_djf_mean=level_mean(msf_djf(:,{600:400},:))
    msf_mean_annual=(/msf_djf_mean(:,:)/)*1.0
    delete(msf_djf_mean)

    ; msf_mean_annual=(/msf_djf(:,{500},:)/)*1.0    ;;;取冬季的情况
       ; levind=ind(msf_mean&level.eq.500)

       ; msf_mean_annual=(/month_to_annual(msf_mean(time|:,level|levind,lon|:), 1)/)
    
    delete(msf_son)
    delete(msf_djf)
    delete(msf_jja)
    msf_mean_annual!0="time"
    msf_mean_annual!1="lon"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    msf_mean_annual&lon=msf_mean&lon
    
    delete(msf_mean)

    printVarSummary(msf_mean_annual)
    msf_mean_annual!0="time"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    delete(timeind)
    
    location=search_location(msf_mean_annual)
    delete(msf_mean_annual)
    locationxi=location[0]
    locationdong=location[1]
    
   printVarSummary(locationxi_rea)
   print("year="+year)
;   if model.eq."HadGEM2-CC".or.model.eq."HadGEM2-ES" then
;     locationxi_rea(i,0:year)=locationxi
;     locationdong_rea(i,0:year)=locationdong
;    else
;     locationxi_rea(i,:)=locationxi
;     locationdong_rea(i,:)=locationdong
; end if
    ;;;;

    rc_xi=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),locationxi , 0, 0)
    rc_xi=(/rc_xi*10/)
    ;;;;
    df   = rc_xi@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
    rc_xi_total(i)=rc_xi
    rc_xi_yield(i)=rc_xi@rstd*T*10.0
    ; rc_xi_yield(i)=rc_xi@rstd*T


    rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
    rc_dong=(/rc_dong*10/)
    
    rc_dong_total(i)=rc_dong
    rc_dong_yield(i)=rc_dong@rstd*T*10.0
    ; rc_dong_yield(i)=rc_dong@rstd*T

    width=locationdong - locationxi

    time=locationxi&time

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
    rc_width=(/rc_width*10/)
   
    rc_width_total(i)=rc_width
    rc_width_yield(i)=rc_width@rstd*T*10.0
    ; rc_width_yield(i)=rc_width@rstd*T
     delete(locationxi)
     delete(locationdong)
     delete(width)
     delete(time)
     delete(rc_width)
    end do

locationxi_rea!0="model"
locationxi_rea!1="time"
locationxi_rea&model=models
locationxi_rea&time=ispan(year1,year2,1)
copy_VarCoords(locationxi_rea, locationdong_rea)
printVarSummary(locationxi_rea)
printVarSummary(locationdong_rea)


nmodel=dimsizes(locationxi_rea&model)
ntime=dimsizes(locationxi_rea&time)

opt=diri+"his_msf_locationxi_trend.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",rc_xi_total))

opt1=diri+"his_msf_locationxi_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_xi_yield))  

opt2=diri+"his_msf_locationdong_trend.ascii"    
asciiwrite(opt2, sprintf ("%6.2f,",rc_dong_total))

opt3=diri+"his_msf_locationdong_yield.ascii"    
asciiwrite(opt3, sprintf ("%6.2f,",rc_dong_yield))  

opt4=diri+"his_msf_locationwidth_trend.ascii"    
asciiwrite(opt4, sprintf ("%6.2f,",rc_width_total))

opt5=diri+"his_msf_locationwidth_yield.ascii"    
asciiwrite(opt5, sprintf ("%6.2f,",rc_width_yield))  

opt6=diri+"his_msf_locationwidth_trend_sst.ascii"    
asciiwrite(opt4, sprintf ("%6.2f,",rc_width_total))




;***************************************************
;写入nc文件
;*************************************************
    diro=diri
    filo1="location_cmip_r1i1p1_djf.nc"
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
    filevardef(fout, "time", typeof(locationxi_rea&time),(/"time"/))

    filevardef(fout, "locationxi", typeof(locationxi_rea),(/"model","time"/))
    filevardef(fout, "locationdong", typeof(locationxi_rea),(/"model","time"/))

    
    filevarattdef(fout, "locationxi", locationxi_rea)
    filevarattdef(fout, "locationdong", locationdong_rea)

    
    setfileoption(fout,"DefineMode",False)

    fout->model=(/ispan(1,dim,1)/)
    fout->time=(/locationxi_rea&time/)
    fout->locationxi=(/locationxi_rea/)
    fout->locationdong=(/locationdong_rea/)

end