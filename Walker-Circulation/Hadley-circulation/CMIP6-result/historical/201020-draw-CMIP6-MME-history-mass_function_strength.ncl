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
wks = gsn_open_wks("pdf", "MME的CMIP6模式东太平洋哈德来环流强度的时间序列")



;;;;
   
    yearstart = 190001
    yearend = 201412
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/msf_mme_history_r1i1p1f1_185001-201412.nc"/)

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    msf = f1->msf(timeind,:,:,:)

    printVarSummary(msf)

  
    msf_mean = dim_avg_n_Wrap(msf(:,:,:,{240:270}), (/3/))



    msf_mean_annual = month_to_annual(msf_mean, 1)




    ;;;;;区域平均
    
    msf_mean_annual_vermean = dim_avg_n_Wrap(msf_mean_annual(:,:,{-10:10}), (/1,2/))
   
    msf_mean_annual_vermean = (/-1*msf_mean_annual_vermean/)

   






    ; ;;;计算趋势
    ; ; rc = msf_mean_annual(1,:,:)
    ; rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual(:,:,{-60:60}),0,0)

    ; copy_VarCoords(msf_mean_annual(1,:,{-60:60}), rc)

    ; rc  = (/rc*dim(0)/)


    ; prob := rc_test(rc)

    ; copy_VarCoords(msf_mean_annual(1,:,{-60:60}), prob)
    ; prob&level@units = "hPa"
    
    res  = True
    res@vpHeightF = 0.3
    res@vpWidthF = 0.7
    res@tiYAxisString = "strenth"
    x = ispan(1900,2014,1)



    plot = gsn_csm_xy(wks, x, msf_mean_annual_vermean, res)
    ;;;  
end 



