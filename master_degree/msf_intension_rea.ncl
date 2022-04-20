
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
;画时间序列图的时候，把BCC-CSM1-1给去掉
models=(/"ens","20cr","era","jra55","MERRA","ncep1","ncep2"/)
; models=(/"ens","era","ncep1"/)

dim=dimsizes(models)

rc_inten_total =new(dim, float)
rc_inten_yield = new(dim, float)
rc_inten_y = new(dim, float)
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


;;;id
rc_inten_id_total =rc_inten_total
rc_inten_id_yield = rc_inten_yield
;;;

;;;alt
rc_inten_alt_total =rc_inten_total
rc_inten_alt_yield = rc_inten_yield
;;;
;;;
rc_inten_yield_90=rc_inten_yield
rc_inten_id_yield_90 = rc_inten_yield
rc_inten_alt_yield_90 = rc_inten_yield

prob=rc_inten_yield
prob_id=rc_inten_yield
prob_alt=rc_inten_yield

;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")
  diri="/mnt/f/share/rea_msf/"

 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_global.nc")

;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    ; timeind=ind(date.ge.197901.and.date.le.200812)
    ; msf_time=f1->msf_time({850:500},:,timeind)
    msf_time_r= f1->msf_time(:,:,timeind)
    delete(f1)

    msf_time_r=(/msf_time_r/10^9/)
    msf_time=tofloat(msf_time_r)
    copy_VarMeta(msf_time_r, msf_time)
    printVarSummary(msf_time)
    delete(msf_time_r)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    delete(date)
    ; msf_mean=level_mean(msf_time(:,:,:))
    ; msf_mean=dim_avg_n_Wrap(msf_time(time|:,level|:,lon|:),1)
    msf_mean=msf_time
    printVarSummary(msf_mean)
    delete(msf_time)
 
 
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_mean(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_mean(time|:,level|:,lon|:), "DJF")
    msf_jja=month_to_season(msf_mean(time|:,level|:,lon|:), "JJA")
    msf_mam=month_to_season(msf_mean(time|:,level|:,lon|:), "MAM")
  
    printVarSummary(msf_son)

    msf_son&time=ispan(year1,year2,1)
    msf_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    msf_mam&time=ispan(year1,year2,1)
    
    ; ;注意点
    ; msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    ; msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_mean(:,:,{yearend})/)
    ; ;;
    ; msf_mean_annual=new((/dimsizes(timeind)/12,dimsizes(msf_mean&lon)/),typeof(msf_mean))
    ; msf_mean_annual=(/msf_son(:,{500},:)+msf_djf(:,{500},:)/)*0.5

     msf_annual=msf_jja
    msf_annual=(/(msf_djf + msf_mam)/2.0/)
    ; msf_annual=(/(msf_son + msf_jja)/2.0/)

    msf_annual=(/month_to_annual(msf_mean(time|:,level|:,lon|:), 1)/)
      delete(msf_mean)
    

    printVarSummary(msf_annual)
    msf_mean_annual=level_mean(msf_annual(:,:,:))  ;time,lon
    ; msf_mean_annual=level_mean(msf_annual(:,{600:400},:))  ;time,lon
    delete(msf_annual)
    delete(msf_mam)
    delete(msf_son)
    delete(msf_djf)
    delete(msf_jja)
   
   printVarSummary(msf_mean_annual)
 

    msf_inten_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{150:210}), 1)
    ; msf_inten_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{120:180}), 1)
    ; msf_inten_id_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{70:100}), 1)
    msf_inten_id_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{60:120}), 1)
    msf_inten_alt_aunnal=dim_avg_n_Wrap(msf_mean_annual(:,{320:340}), 1)
    delete(msf_mean_annual)
    ;;;calculate trend

    msf_inten_rc= regline(ispan(1, dimsizes(msf_inten_aunnal&time), 1), msf_inten_aunnal)
    printVarSummary(msf_inten_rc)
    df   = msf_inten_rc@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
    rc_inten_total(i)=(/msf_inten_rc/)*10.0
    rc_inten_yield_90(i)=msf_inten_rc@rstd*T*10.0
    rc_inten_yield(i)=msf_inten_rc@rstd*10.0
    rc_inten_y(i)=msf_inten_rc@yintercept

    prob(i) = ( 1-betainc(df/(df+msf_inten_rc@tval^2), df/2.0, 0.5) )


     msf_inten_id_rc= regline(ispan(1, dimsizes(msf_inten_id_aunnal&time), 1), msf_inten_id_aunnal)
    df   = msf_inten_id_rc@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
    rc_inten_id_total(i)=(/msf_inten_id_rc/)*10.0
    rc_inten_id_yield_90(i)=msf_inten_id_rc@rstd*T*10.0
    rc_inten_id_yield(i)=msf_inten_id_rc@rstd*10.0
    prob_id(i) = ( 1-betainc(df/(df+msf_inten_id_rc@tval^2), df/2.0, 0.5) )

    msf_inten_alt_rc= regline(ispan(1, dimsizes(msf_inten_alt_aunnal&time), 1), msf_inten_alt_aunnal)
    df   = msf_inten_alt_rc@nptxy-2
    T = cdft_t(p, df)    ; T =  2.12
    rc_inten_alt_total(i)=(/msf_inten_alt_rc/)*10.0
    rc_inten_alt_yield_90(i)=msf_inten_alt_rc@rstd*T*10.0
    rc_inten_alt_yield(i)=msf_inten_alt_rc@rstd*10.0
    prob_alt(i) = ( 1-betainc(df/(df+msf_inten_alt_rc@tval^2), df/2.0, 0.5) )
    ;;;
     inten(i,:)=msf_inten_aunnal
    inten_id(i,:)=msf_inten_id_aunnal
    inten_alt(i,:)=msf_inten_alt_aunnal

    delete([/msf_inten_id_aunnal,msf_inten_alt_aunnal,msf_inten_aunnal/])

    end do

opt=diri+"rea_inten_trend.ascii"    
asciiwrite(opt, sprintf ("%6.3f,",rc_inten_total))

opt1=diri+"rea_inten_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_yield))

opt=diri+"rea_inten_id_trend.ascii"    
asciiwrite(opt, sprintf ("%6.3f,",rc_inten_id_total))

opt1=diri+"rea_inten_id_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_id_yield))

opt=diri+"rea_inten_alt_trend.ascii"    
asciiwrite(opt, sprintf ("%6.3f,",rc_inten_alt_total))

opt1=diri+"rea_inten_alt_yield.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_alt_yield))

;;;;
opt1=diri+"rea_inten_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_yield_90))

opt1=diri+"rea_inten_id_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_id_yield_90))

opt1=diri+"rea_inten_alt_yield_90.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",rc_inten_alt_yield_90))
;;;
opt1=diri+"rea_inten_prob.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",prob))

opt1=diri+"rea_inten_id_prob.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",prob_id))

opt1=diri+"rea_inten_alt_prob.ascii"    
asciiwrite(opt1, sprintf ("%6.3f,",prob_alt))



;;；画图


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



;***************************************************
;写入nc文件,年的文件
;*************************************************
    diro=diri
    filo1="intension_rea_r1i1p1_djf.nc"
    ; filo1="intension_rea_r1i1p1_jja.nc"
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
