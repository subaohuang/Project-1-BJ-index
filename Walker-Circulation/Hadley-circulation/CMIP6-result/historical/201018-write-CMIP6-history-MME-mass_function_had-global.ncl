function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end


begin
;读入数据

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;     "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;     "GFDL-CM4", "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         "HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM"/)
model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-ESM1",\
              "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
   , "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
  "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
        ,"HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
        "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM","SAM0-UNICON","UKESM1-0-LL"/)

; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(1,graphic)
iplot= plot
plot1 = plot


a=6371000
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180

diri_output="/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/"
wks = gsn_open_wks("pdf", "MME的CMIP6模式东太平洋哈德来环流的趋势")

msf_total = new((/dim_model,1980,11,72,144/),float)
do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart = 185001
    yearend = 201412
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/1850-1949/msf_"+model+"_history_r1i1p1f1_185001-194912.nc",\
               "/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/1950-2014/msf_"+model+"_history_r1i1p1f1_195001-201412.nc"/)

    f1 := addfiles(fils1, "r")
    
    time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    msf_r := f1[:]->msf(timeind,:,:,:)
    msf_r = (msf_r/10^10)
    printVarSummary(msf_r)

    msf_total(k,:,:,:,:) = msf_r
     print("********************************finish calculate"+" "+model+"************************************")

end  do 


    msf = dim_avg_n_Wrap(msf_total, 0)

    printVarSummary(msf)
    ntime = dimsizes(msf&time)
    nlev = dimsizes(msf&level)
    nlat = dimsizes(msf&lat)
    nlon = dimsizes(msf&lon)

  print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="msf_mme_history_r1i1p1f1_185001-201412.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "CMIP6 historical run Mass stream Function for Hadley"    
    fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time",typeof(msf&time), (/"time"/))
    filevardef(fout, "level",typeof(msf&level), (/"level"/))
    filevardef(fout, "lat", typeof(msf&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(msf&lon),(/"lon"/))

   
    filevardef(fout, "msf", typeof(msf),(/"time","level","lat","lon"/))

   
    filevarattdef(fout, "level", msf&level)
    filevarattdef(fout, "lat", msf&lat)
    filevarattdef(fout, "lon", msf&lon)
    filevarattdef(fout, "time", msf&time)
    filevarattdef(fout, "msf", msf)
 
    
    setfileoption(fout,"DefineMode",False)


    fout->time=(/msf&time/)
    fout->level=(/msf&level/)
    fout->lat=(/msf&lat/)
    fout->lon=(/msf&lon/)
    fout->msf = (/msf/)
    delete(msf)
 print("********************************finish calculate"+" "+model+"************************************")
end 