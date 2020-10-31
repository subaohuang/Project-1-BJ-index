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

rc_all = new((/dim_model,72,144/),float)
rc_all_u = new((/dim_model,72,144/),float)
rc_all_v = new((/dim_model,72,144/),float)


wks = gsn_open_wks("pdf","201026-MME的CESM模式的SST趋势的水平分布")         ; send graphics to PNG file
do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start read"+" "+model+"**********************************************")
    yearstart = 190001
    yearend = 201412


    fils1 := systemfunc("ls /home/yangsong3/data-observation/cmip6/"+model+"/historical/tos_fixed/tos_Omon_"+"*")
    fils2 := systemfunc("ls /home/yangsong3/data-observation/cmip6/"+model+"/historical/uv_fixed/ua_Amon_"+"*")
    fils3 := systemfunc("ls /home/yangsong3/data-observation/cmip6/"+model+"/historical/uv_fixed/va_Amon_"+"*")
    ; print(fils1)

   if dimsizes(fils1).eq.1 then 

        f1 := addfile(fils1, "r")

        time:=f1->time    
        date:=cd_calendar(time, 1)
        delete(time)

        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        delete(date)     
        sst = lonFlip(f1->tos(timeind,:,:))
    else 
        f1 := addfiles(fils1, "r")

        time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
        date:=cd_calendar(time, 1)
        delete(time)
        timeind:=ind(date.le.yearend.and.date.ge.yearstart)
        delete(date)     
        sst = lonFlip(f1[:]->tos(timeind,:,:))
    end if

    if k.eq.12 then 
    sst =  where(sst.ge.40.0 , sst@_FillValue , sst)
    end if
        delete(timeind)
        ssta :=(/rmMonAnnCycTLL(sst)/)
        copy_VarCoords(sst, ssta)
        printVarSummary(ssta)


;;;;


    f2 = addfiles(fils2, "r")
    f3 = addfiles(fils3, "r")

    time=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)


    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    time_r = time(timeind)

    ; z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,150.0,100.0/)*100.0

    ; levind=get1Dindex(plev,z0)
  

    uwnd := lonFlip(f2[:]->ua(timeind,0,:,:))  ;;;195801-201212
    vwnd := lonFlip(f3[:]->va(timeind,0,:,:))  ;;;
    
    delete(f1)
    delete(f2)
    delete(f3)
    

    ;;;计算趋势
    ssta_annual = month_to_annual(ssta, 1)
    dim = dimsizes(ssta_annual)
    rc := regCoef_n(ispan(1,dim(0),1),ssta_annual,0,0)
    copy_VarCoords(ssta_annual(1,:,:), rc)

    rc  = (/rc*dim(0)/)
    rc_all(k,:,:) = rc
    ;;;;

    ;;;;计算趋势
    uwnd_annual = month_to_annual(uwnd, 1)
    dim = dimsizes(uwnd_annual)
    rc_uwnd := regCoef_n(ispan(1,dim(0),1),uwnd_annual,0,0)
    copy_VarCoords(uwnd_annual(1,:,:), rc_uwnd)

    rc_uwnd  = (/rc_uwnd*dim(0)/)
    rc_all_u(k,:,:) = rc_uwnd
 
   ;;;;计算趋势
    vwnd_annual = month_to_annual(vwnd, 1)
    dim = dimsizes(vwnd_annual)
    rc_vwnd := regCoef_n(ispan(1,dim(0),1),vwnd_annual,0,0)
    copy_VarCoords(vwnd_annual(1,:,:), rc_vwnd)

    rc_vwnd  = (/rc_vwnd*dim(0)/)
    rc_all_v(k,:,:) = rc_vwnd
;;;;;
end do 
printVarSummary(rc_all_u)
rc_all!0 = "model_name"
rc_all&model_name = ispan(1,30,1)

rc_all!1 = "lat"
rc_all&lat = uwnd_annual&lat

rc_all!2 = "lon"
rc_all&lon = uwnd_annual&lon

rc_all_u!0 = "model_name"
rc_all_u&model_name = ispan(1,30,1)

rc_all_u!1 = "lat"
rc_all_u&lat = uwnd_annual&lat

rc_all_u!2 = "lon"
rc_all_u&lon = uwnd_annual&lon


rc_all_v!0 = "model_name"
rc_all_v&model_name = ispan(1,30,1)

rc_all_v!1 = "lat"
rc_all_v&lat = uwnd_annual&lat

rc_all_v!2 = "lon"
rc_all_v&lon = uwnd_annual&lon

outpath = "/home/ys17-19/lsh/Project/Walker-Circulation/Hadley-circulation/CMIP6-result/historical/"
outfile = "sst_u_v_1000hpa_trend_1900-2014.nc"

system("/bin/rm -f "+outpath+outfile)
ncout = addfile(outpath+outfile, "c")

ncout->rc_sst = rc_all
ncout->rc_u = rc_all_u
ncout->rc_v = rc_all_v
; ncout->model_name = model_name


end 


 


  