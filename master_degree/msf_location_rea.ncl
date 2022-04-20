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

; models=(/"mre","20cr","era","jra55","MERRA","ncep1","ncep2"/)
models=(/"20cr","era","jra55","MERRA","ncep1","ncep2","mre"/)

dim=dimsizes(models)
rc_xi_total=new(dimsizes(models), float)
rc_xi_yield=new(dimsizes(models), float)

rc_dong_total=new(dimsizes(models), float)
rc_dong_yield=new(dimsizes(models), float)

rc_width_total=new(dimsizes(models), float)
rc_width_yield=new(dimsizes(models), float)

prob_xi=rc_xi_total
prob_dong=rc_dong_total
prob_width=rc_width_total

rc_xi_yield_90=rc_xi_yield
rc_dong_yield_90=rc_dong_yield
rc_width_yield_90=rc_width_yield



yearstart=197901
yearend=201212  
year1=1979
year2=2012
year=year2-year1+1

locationxi_rea=new((/dim,year/), float)
locationdong_rea=new((/dim,year/), float)

; p  = 0.975
p=0.95
;;;循环
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)

 diri="/mnt/f/share/rea_msf/"

 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_global.nc")

;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    year=dimsizes(timeind)/12-1
    ; timeind=ind(date.ge.197901.and.date.le.201212)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    ; msf_time=f1->msf_time({500},:,timeind)
    msf_time=f1->msf_time(:,:,timeind)
    levind=ind(msf_time&level.eq.500)
    delete(f1)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    printVarSummary(msf_time)
    msf_time=(/msf_time/10^9/)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    ; msf_mean=msf_time(time|:,lon|:)
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
    msf_mam&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    

    ; ;注意点
    ; msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    ; msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ; ;;
    msf_djf=(/(msf_djf + msf_mam)/2.0/)


    msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
    printVarSummary(msf_mean_annual)

    ; msf_mean_annual=month_to_annual(msf_mean, 1)
     msf_mean_annual=(/msf_djf(:,{500},:)/)*1.0
     ; msf_mean_annual=(/msf_jja(:,{500},:)/)*1.0
    ; msf_mean_annual=(/month_to_annual(msf_mean(time|:,level|levind,lon|:), 1)/)
     
    ; msf_djf_mean=level_mean(msf_djf(:,{600:400},:))
    ; msf_mean_annual=(/msf_djf_mean(:,:)/)*1.0
    ; delete(msf_djf_mean)


    delete(msf_son)
    delete(msf_djf)
    delete(msf_jja)
    delete(msf_mam)
    msf_mean_annual!0="time"
    msf_mean_annual!1="lon"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    msf_mean_annual&lon=msf_mean&lon

    printVarSummary(msf_mean_annual)
    msf_mean_annual!0="time"
    msf_mean_annual&time=ispan(0,dimsizes(timeind)/12-1,1)
    delete(msf_mean)
    
    location=search_location(msf_mean_annual)
    locationxi= dble2flt(location[0])
    locationdong=dble2flt(location[1])

    locationxi_rea(i,:)=locationxi
    locationdong_rea(i,:)=locationdong

    ;;;;
    delete(msf_mean_annual)

    rc_xi=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),locationxi , 0, 0)
    rc_xi=(/rc_xi*10/)

    ;;;;
    df   = rc_xi@nptxy-2


    T = cdft_t(p, df)    ; T =  2.12
    rc_xi_total(i)=rc_xi
    rc_xi_yield_90(i)=rc_xi@rstd*T*10.0
    rc_xi_yield(i)=rc_xi@rstd*10.0
    prob_xi(i) = ( 1-betainc(df/(df+rc_xi@tval^2), df/2.0, 0.5) )
    delete(rc_xi)
    

    rc_dong=regCoef_n(ispan(0,dimsizes(locationdong)-1,1),locationdong, 0, 0)
    rc_dong=(/rc_dong*10/)
    
    rc_dong_total(i)=rc_dong
    rc_dong_yield_90(i)=rc_dong@rstd*T*10.0
    rc_dong_yield(i)=rc_dong@rstd*10.0
    prob_dong(i) = (1 - betainc(df/(df+rc_dong@tval^2), df/2.0, 0.5) )
    ; rc_dong_yield(i)=rc_dong@rstd*T
    delete(rc_dong)
    width=locationdong - locationxi

    time=locationxi&time

    rc_width=regCoef_n(ispan(0,dimsizes(locationxi)-1,1),width , 0, 0)
  
   
    rc_width_total(i)=rc_width*10
    rc_width_yield(i)=rc_width@rstd*10.0
    rc_width_yield_90(i)=rc_width@rstd*T*10.0

    prob_width(i) = (1 - betainc(df/(df+rc_width@tval^2), df/2.0, 0.5) )
    ; rc_width_yield(i)=rc_width@rstd*T
    delete(rc_width)
    delete(locationdong)
    delete(locationxi)
    print(i)
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

opt=diri+"rea_msf_locationxi_trend.ascii"    
asciiwrite(opt, sprintf ("%6.2f,",rc_xi_total))

opt1=diri+"rea_msf_locationxi_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_xi_yield))  

opt2=diri+"rea_msf_locationdong_trend.ascii"    
asciiwrite(opt2, sprintf ("%6.2f,",rc_dong_total))

opt3=diri+"rea_msf_locationdong_yield.ascii"    
asciiwrite(opt3, sprintf ("%6.2f,",rc_dong_yield))  

opt4=diri+"rea_msf_locationwidth_trend.ascii"    
asciiwrite(opt4, sprintf ("%6.2f,",rc_width_total))

opt5=diri+"rea_msf_locationwidth_yield.ascii"    
asciiwrite(opt5, sprintf ("%6.2f,",rc_width_yield))  


opt6=diri+"rea_msf_locationxi_prob.ascii"    
asciiwrite(opt6, sprintf ("%6.2f,",prob_xi))

opt7=diri+"rea_msf_locationdong_prob.ascii"    
asciiwrite(opt7, sprintf ("%6.2f,",prob_dong))

opt8=diri+"rea_msf_locationwidth_prob.ascii"    
asciiwrite(opt8, sprintf ("%6.2f,",prob_width))
;;;;

opt1=diri+"rea_msf_locationxi_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_xi_yield_90))  

opt1=diri+"rea_msf_locationdong_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_dong_yield_90))  

opt1=diri+"rea_msf_locationwidth_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.2f,",rc_width_yield_90))  


;***************************************************
;写入nc文件
;*************************************************
    diro=diri
    filo1="location_rea_r1i1p1_djf.nc"
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