  begin
;;读入数据

    
   ; model=(/"INMCM4","IPSL-CM5A-MR",\
   ; "MIROC5","MPI-ESM-MR"/)
   model=(/"CSIRO-Mk3-6-0"/)
   
   do i=0,dimsizes(model)-1
    diri_output="/mnt/i/CMIP5/RCP45/uv/"
 
    yearstart=200601
    yearend=201212
     
    
    if model(i).eq."FGOALS-s2".or.model(i).eq."HadGEM2-CC".or.model(i).eq."HadGEM2-ES" then

    yearstart=197901
    yearend=200412

    end if

    fils1=systemfunc("ls /mnt/i/CMIP5/RCP45/ua/"+model(i)+"/r1i1p1/ua_Amon_"+model(i)+"_rcp45_r1i1p1_"+"*"+".nc")
    fils2=systemfunc("ls /mnt/i/CMIP5/RCP45/va/"+model(i)+"/r1i1p1/va_Amon_"+model(i)+"_rcp45_r1i1p1_"+"*"+".nc")

    print(fils1)
    print(fils2)
    
    f1=addfiles(fils1, "join")
    f2=addfiles(fils2, "join")
    
    delete(fils1)
    delete(fils2)

    time=f1[:]->time    
    printVarSummary(time)
    date=cd_calendar(time, 1)
  
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    delete(date)
    print("读入风场")
      if model(i).eq."CSIRO-Mk3-6-0" then 
    z0=(/925.0/)*100.0
    else
    z0=(/1000.0/)*100.0
end if
    levind=get1Dindex(f1[0]->plev,z0)

    uwnd_r=f1[:]->ua(timeind,levind,:,:)  ;;;195801-201212
    vwnd_r=f2[:]->va(timeind,levind,:,:)  ;;;
   delete(timeind)
    delete(levind)
    delete(f1)
    delete(f2)

    xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)

     print("*****************开始水平插值************************")  
    uwnd=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, \
                      True, xo, yo, 0)

    vwnd=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, \
                      True, xo, yo, 0)
     ; uwnd!0="time"
     uwnd!1="lat"
     uwnd!2="lon"
     
     ; vwnd!0="time"
     vwnd!1="lat"
     vwnd!2="lon"
   
 
    uwnd&lat@units="degrees_north"
    uwnd&lon@units="degrees_east"

    vwnd&lat@units="degrees_north"
    vwnd&lon@units="degrees_east"

    delete(uwnd_r)
    delete(vwnd_r)

    printVarSummary(uwnd)
    printVarSummary(vwnd)

    lat=uwnd&lat
    lon=uwnd&lon
  
   ;;;interp missing value
     print("****************interp missing value*********************")
       
    ;******************************************************************
    ; use poisson_grid_fill to fill missing value
    ;******************************************************************
    nscan = 1500           ; usually *much* fewer
    eps   = 0.001            ; variable depended
    gtype = False            ; "gendat" does not generate cyclic fields
    guess = 0                ; use zonal means
    relc  = 0.6              ; standard relaxation coef
    opt   = 0
    
    poisson_grid_fill(uwnd, gtype, guess, nscan, eps, relc, opt)
    
    poisson_grid_fill(vwnd, gtype, guess, nscan, eps, relc, opt)



    ntime=dimsizes(uwnd&time)
    nlat=dimsizes(uwnd&lat)
    nlon=dimsizes(uwnd&lon)
  
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="uv_"+model(i)+"_rcp45_r1i1p1_200601-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "u,v插值后的结果，取1000hpa"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(uwnd&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(uwnd&lon),(/"lon"/))
    filevardef(fout, "time", typeof(uwnd&time),(/"time"/))

    filevardef(fout, "uwnd", typeof(uwnd),(/"time","lat","lon"/))
    filevardef(fout, "vwnd", typeof(vwnd),(/"time","lat","lon"/))


    filevarattdef(fout, "lat", uwnd&lat)
    filevarattdef(fout, "lon", uwnd&lon)
    filevarattdef(fout, "time", time)
    
    filevarattdef(fout, "uwnd", uwnd)
    filevarattdef(fout, "vwnd", vwnd)
    
    setfileoption(fout,"DefineMode",False)

    fout->lat=(/uwnd&lat/)
    fout->lon=(/uwnd&lon/)
    fout->time=(/uwnd&time/)
    fout->uwnd=(/uwnd/)
    fout->vwnd=(/vwnd/)
    delete(uwnd)
    delete(vwnd)
    delete(ntime)
    delete(nlon)
    delete(nlat)
    delete(time)
 end do

end 



