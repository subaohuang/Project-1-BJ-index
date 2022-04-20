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

model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-ESM1",\
              "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
   , "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
  "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
        ,"HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
        "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM","SAM0-UNICON","UKESM1-0-LL"/)

; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(1,graphic)


rc_all_v = new((/dim_model,115*12/),float)


do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start read"+" "+model+"**********************************************")
    yearstart = 190001
    yearend = 201412

    fils3 := systemfunc("ls /home/yangsong3/data-observation/cmip6/"+model+"/historical/uv_fixed/va_Amon_"+"*")
    ; print(fils1)

    f3 = addfiles(fils3, "r")

    time :=f3[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)


    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    time_r := time(timeind)

    vwnd := lonFlip(f3[:]->va(timeind,0,:,:))  ;;;
    delete(f3)
    


   ;;;;计算趋势
   
    
    vep = dim_avg_n_Wrap(vwnd(:,{-10:10},{150:270}), (/1,2/))
    rc_all_v(k,:) = vep
end do 
printVarSummary(rc_all_v)

rc_all_v!0 = "model_name"
rc_all_v&model_name = ispan(1,30,1)

rc_all_v!1 = "time"
rc_all_v&time = vep&time


outpath = "/home/ys17-19/lsh/Project/Walker-Circulation/Hadley-circulation/CMIP6-result/historical/"
outfile = "Vep_timeseries_1900-2014.nc"

system("/bin/rm -f "+outpath+outfile)
ncout = addfile(outpath+outfile, "c")

ncout->vep_all = rc_all_v



end 


 


  