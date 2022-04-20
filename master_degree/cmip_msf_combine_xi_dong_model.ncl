
begin

; FGOALS-g2和s2没有海温数据l;;同时添加了models=(/"CanESM2",GISS-E2-R","MIROC-ESM-CHEM","NorESM1-ME"/)
models=(/"ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","CCSM4","CanESM2","CESM1-BGC","CMCC-CM","CMCC-CMS","CNRM-CM5",\
        "CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-R",\
         "HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4","IPSL-CM5A-LR",\
    "IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC5","MIROC-ESM-CHEM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)
;;;共30个模式


 diri1="/mnt/f/share/CMIP5/msf_history/"
 diri2="/mnt/f/share/CMIP5/msf_rcp45/"

 diri_output="/mnt/f/share/CMIP5/msf_history/"
 dim=dimsizes(models)
 print("dim="+dim)
 width_trend=asciiread(diri1+"cmip_msf_inten_trend.ascii" ,dim, "float")
 sst_trend=asciiread("/mnt/f/share/CMIP5/sst_cha/his_sstca_trend.ascii" ,dim, "float")

 sst_poid=models(ind(sst_trend.ge.0.05))
 width_poid=models(ind(width_trend.ge.0.1))
 sst_neid=models(ind(sst_trend.le.-0.1))
 width_neid=models(ind(width_trend.le.-0.1))

 ; xi_model =sst_poid(get1Dindex(sst_poid,width_poid))
 ; dong_model =sst_neid(get1Dindex(sst_neid,width_neid))
 
 xi_model= width_poid
 dong_model=width_neid

 msf_time_xi=new((/dimsizes(xi_model),12,144,408/),float)
 msf_time_dong=new((/dimsizes(dong_model),12,144,408/),float)


do j=0,1

	if j.eq.0 then 
	model=xi_model
	else
    model=dong_model
	end if

do  i=0,dimsizes(model)-1
    print("*************************"+model(i)+"***************************")
	fils1=diri1+"msf_time_"+model(i)+"_history_r1i1p1_197901-200512.nc"
	fils2= diri2+"msf_time_"+model(i)+"_rcp45_r1i1p1_200601-201212.nc"

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
    f=addfiles(fils, "r")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f[:]->time,1)
    date1=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)
    

    timeind=ind(date.ge.197901.and.date.le.201212)
    timeind1=ind(date1.le.yearend1.and.date1.ge.yearstart1)
    timeind2=ind(date2.le.yearend2.and.date2.ge.yearstart2)


    delete([/date1 ,date2, date/])
   
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

    delete(msf_time1)
    delete(msf_time2)
    delete(timeind1)
    delete(timeind2)
	printVarSummary(msf_time_r)
    msf_time=tofloat(msf_time_r)
    copy_VarMeta(msf_time_r, msf_time)
    delete(msf_time_r)
    printVarSummary(msf_time)
 
	if j.eq.0 then
	msf_time_xi(i,:,:,:)=msf_time(:,:,timeind) 
    else 
	msf_time_dong(i,:,:,:)=msf_time(:,:,timeind)
    end if
	delete(msf_time)
    delete(timeind)
	end do
    
    if j.eq.0 then
	msf_xi_model=dim_avg_n_Wrap(msf_time_xi, 0)
	printVarSummary(msf_xi_model)
    else
	msf_dong_model=dim_avg_n_Wrap(msf_time_dong, 0)
	printVarSummary(msf_dong_model)
    end if
    
    
	time=msf_xi_model&time
	level=msf_xi_model&level
	lon=msf_xi_model&lon

    ntime=dimsizes(time)
    nlev=dimsizes(level)
    nlon=dimsizes(lon)
       
   print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    if j.eq.0  then
    filo1="msf_time_xi_cmip_r1i1p1_197901-201212.nc"
    else
    filo1="msf_time_dong_cmip_r1i1p1_197901-201212.nc"
    end if
    
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
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
    filevardef(fout, "level",typeof(level), (/"level"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(time),(/"time"/))
  
    filevardef(fout, "msf_time", typeof(msf_xi_model),(/"level","lon","time"/))
 

    filevarattdef(fout, "level", level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "time", time)
    filevarattdef(fout, "msf_time",msf_xi_model)
    setfileoption(fout,"DefineMode",False)

    fout->level=(/level/)
    fout->lon=(/lon/)
    fout->time=(/time/)
    if j.eq.0 then
    fout->msf_time=(/msf_xi_model/)
    else
    fout->msf_time=(/msf_dong_model/)
    end if
    delete(time)
    delete(level)
    delete(lon)
    delete(fout)
    delete(model)
end do
end 