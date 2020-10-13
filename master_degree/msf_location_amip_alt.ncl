
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
undef("search_location_alt") ;;;altlantic
function  search_location_alt(msf_mean_annual)
local msf_mean_annual_xi , msf_mean_annual_dong
begin
locationxi=new(dimsizes(msf_mean_annual&time),typeof(msf_mean_annual&lon))
locationxi!0="time"
locationxi&time=ispan(0,dimsizes(msf_mean_annual&time)-1,1)
locationdong=locationxi
; msf_mean_annual_xi=msf_mean_annual(:,{-90:-30})
; msf_mean_annual_dong=msf_mean_annual(:,{-40:10})
msf_mean_annual_xi=msf_mean_annual(:,{270:330})   ;;;
msf_mean_annual_dong=msf_mean_annual(:,{320:370})
;;;judge locationxi
do t=0,dimsizes(msf_mean_annual&time)-1
    do i=0,dimsizes(msf_mean_annual_xi&lon)-9
        if  (msf_mean_annual_xi(t,i+1).ge.0).and.((msf_mean_annual_xi(t,i).le.0).and.(msf_mean_annual_xi(t,i+2).gt.0))  then 
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
   if ismissing(locationxi(t)).eq.True  then 
   co_in=closest_val(max(msf_mean_annual(t,:))/10,msf_mean_annual_dong(t,:))
   ; test=msf_mean_annual_xi(t,{-70:-30})
   ; co_in=closest_val(0,msf_mean_annual_xi(t,{-70:-30}))
     test=msf_mean_annual_xi(t,{290:330})
   co_in=closest_val(0,msf_mean_annual_xi(t,{290:330}))
   locationxi(t)=test&lon(co_in)
   delete(co_in)
  end if

end do
; locationxi=(/linmsg(locationxi, -1)/)
locationdong=(/linmsg(locationdong, -1)/)
 location = [/ locationxi, locationdong /]   ; return two variables as a type 'list'
return location
end 


;********************************************************
;draw_msf_hov
begin

;画时间序列图的时候，把BCC-CSM1-1给去掉
; models=(/"amip","ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
;     "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)

models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanAM4","CCSM4","CMCC-CM","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FGOALS-s2","GFDL-CM3","GFDL-HIRAM-C180","HadGEM2-A","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-AGCM3-2H","MRI-CGCM3","NorESM1-M"/)

dim=dimsizes(models)
print("dim="+dim)

rc_xi_total=new(dim, float)
rc_xi_yield=new(dim, float)

rc_dong_total=new(dim, float)
rc_dong_yield=new(dim, float)

rc_width_total=new(dim, float)
rc_width_yield=new(dim, float)

rc_xi_yield_90=rc_xi_yield
rc_dong_yield_90=rc_dong_yield
rc_width_yield_90=rc_width_yield

yearstart=197901
yearend=200812  
year1=1979
year2=2008
; p  = 0.975
p=0.95
year=year2-year1+1

locationxi_rea=new((/dim,year/), float)
locationdong_rea=new((/dim,year/), float)
;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 ; model="ACCESS1-0"
 diri="/mnt/f/share/CMIP5/msf_AMIP/"

 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_amip_r1i1p1_197901-200812.nc")

;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    msf_time_r=f1->msf_time(:,:,timeind)
    print("*************"+models(i)+"*************")
    msf_time_f=lonPivot((msf_time_r(time|:,level|:,lon|:)),30)
    delete(msf_time_r)
    ; if models(i).eq."MERRA" then
    ;     msf_time_f&lon=fspan(-180,177.5,143)
    ;  end if
    msf_time=msf_time_f(level|:,lon|:,time|:)
    delete(msf_time_f)
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
    msf_mam=month_to_season(msf_mean(time|:,level|:,lon|:), "MAM")
    printVarSummary(msf_son)
    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均

    ; ;;注意点
    ; msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    ; msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ; ;;;
   
    msf_djf=(/(msf_djf + msf_mam)/2.0/)
    ; msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
    ; printVarSummary(msf_mean_annual)
    
    ; msf_mean_annual=(/msf_jja(:,{500},:)/)*1.0
    ; levind=ind(msf_mean&level.eq.500)

    ; msf_djf_mean=level_mean(msf_djf(:,:,:))
    ; msf_djf_mean=level_mean(msf_djf(:,{600:400},:))
    ; msf_mean_annual=(/msf_djf_mean(:,:)/)*1.0
    ; delete(msf_djf_mean)
    msf_mean_annual=(/msf_djf(:,{500},:)/)*1.0
    
    ; msf_mean_annual=(/month_to_annual(msf_mean(time|:,level|levind,lon|:), 1)/)



    msf_mean_annual!0="time"
    msf_mean_annual!1="lon"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    msf_mean_annual&lon=msf_mean&lon
    delete(msf_mean)
    printVarSummary(msf_mean_annual)
    msf_mean_annual!0="time"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    
    
    location=search_location_alt(msf_mean_annual)
    locationxi=location[0]
    locationdong=location[1]

    locationxi_rea(i,:)=locationxi
    locationdong_rea(i,:)=locationdong
    ;;;;
     

    rc_xi=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),locationxi , 0, 0)

    rc_xi=(/rc_xi*10/)
    ;;;;
    df   = rc_xi@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
    rc_xi_total(i)=rc_xi
    rc_xi_yield_90(i)=rc_xi@rstd*T*10.0
    rc_xi_yield(i)=rc_xi@rstd*10


    rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
    rc_dong=(/rc_dong*10/)
    
    rc_dong_total(i)=rc_dong
    rc_dong_yield_90(i)=rc_dong@rstd*T*10.0
    rc_dong_yield(i)=rc_dong@rstd*10

    width=locationdong - locationxi

    time=locationxi&time

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)

   
    rc_width_total(i)=rc_width*10.0
    rc_width_yield_90(i)=rc_width@rstd*T*10.0
    rc_width_yield(i)=rc_width@rstd*10.0
    delete(msf_son)
    delete(msf_djf)
    delete(msf_jja)
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

opt=diri+"msf_locationxi_alt_trend.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",rc_xi_total))

opt1=diri+"msf_locationxi_alt_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_xi_yield))  

opt2=diri+"msf_locationdong_alt_trend.ascii"    
asciiwrite(opt2, sprintf ("%6.2f,",rc_dong_total))

opt3=diri+"msf_locationdong_alt_yield.ascii"    
asciiwrite(opt3, sprintf ("%6.2f,",rc_dong_yield))  

opt4=diri+"msf_locationwidth_alt_trend.ascii"    
asciiwrite(opt4, sprintf ("%6.2f,",rc_width_total))

opt5=diri+"msf_locationwidth_alt_yield.ascii"    
asciiwrite(opt5, sprintf ("%6.2f,",rc_width_yield))  

;;;
opt1=diri+"msf_locationxi_alt_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_xi_yield_90))  

opt1=diri+"msf_locationdong_alt_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_dong_yield_90))  

opt1=diri+"msf_locationwidth_alt_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_width_yield_90))  



;***************************************************
;写入nc文件
;*************************************************
    diro=diri
    filo1="location_alt_amip_r1i1p1_djf.nc"
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