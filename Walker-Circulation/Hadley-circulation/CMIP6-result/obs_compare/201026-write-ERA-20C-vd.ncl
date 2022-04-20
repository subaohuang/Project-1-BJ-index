  begin
;;读入数据


;;;;
    model = "ERA-20C"
    diri_output="/home/yangsong3/data-observation/linshh/obs-derived/vd/"
    
    yearstart = 190001
    yearend = 201012

    fils1 := systemfunc("ls /home/yangsong3/data-observation/ERA-20C/pressure/uwnd/uwnd.*")
    fils2 := systemfunc("ls /home/yangsong3/data-observation/ERA-20C/pressure/vwnd/vwnd.*")

    f1=addfiles(fils1, "r")
    f2=addfiles(fils2, "r")

    date=cd_calendar(f1[:]->time, 1)
    timeind=ind(date.ge.yearstart.and.date.le.yearend)
    time = f1[:]->time(timeind)
    

    lat = f1[1]->latitude
    lon = f1[1]->longitude
    plev = f1[1]->level
  
    ;;;region 
      lat1 = -90
      lat2 = 90
      lon1 = 0
      lon2 = 360
    ;;;

       
    z0= (/1000.0,   950.0,   900.0,   850.0,   800.0,   750.0,\
    700.0,   650.0,   600.0,   550.0,   500.0,   450.0,\
    400.0,   350.0,   300.0,   250.0,   200.0,   150.0,\
    100.0,   70.00,   50.00,   30.00,   20.00,   10.00/)

    levind=get1Dindex(plev,z0)


    uwnd = new((/dimsizes(time), dimsizes(z0), dimsizes(lat),dimsizes(lon)/),float)
    vwnd = new((/dimsizes(time), dimsizes(z0), dimsizes(lat),dimsizes(lon)/),float)


    do i = 0,dimsizes(timeind)/12-1
      uwnd(i*12:(i+1)*12-1,:,:,:) = short2flt(f1[i]->u(:,levind, {lat1:lat2},{lon1:lon2}))
      vwnd(i*12:(i+1)*12-1,:,:,:) = short2flt(f2[i]->v(:,levind, {lat1:lat2},{lon1:lon2}))
    end do


    delete(f1)
    delete(f2)

      
    print("*********************开始计算辐散风******************************")  
    vrdvx=uv2vrdvF(uwnd,vwnd)
    ; vrdvx=uv2dv_cfd(uwnd, vwnd, uwnd&lat, vwnd&lon, 0)
    uvd=dv2uvF(vrdvx(1,:,:,:,:))
    delete(vrdvx)
    ; ud=uvd(0,:,:,:,:)
    vd=uvd(1,:,:,:,:)  ;;;time lev,lat,lon
    delete(uvd)
    vd!0="time"
    vd!1="level"
    vd!2="lat"
    vd!3="lon"
    copy_VarCoords(vwnd, vd)
    copy_VarAtts(vwnd, vd)
    printVarSummary(vd)
    delete(uwnd)
    delete(vwnd)

    dim = dimsizes(vd)
    ntime = dim(0)
    nlev = dim(1)
    nlat = dim(2)
    nlon = dim(3)
    
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="vd_"+model+"_190001-201012.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "ERA-20C Meridional Divergence wind"    
    fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time",typeof(vd&time), (/"time"/))
    filevardef(fout, "level",typeof(vd&level), (/"level"/))
    filevardef(fout, "lat", typeof(vd&latitude),(/"lat"/))
    filevardef(fout, "lon", typeof(vd&longitude),(/"lon"/))

   
    filevardef(fout, "vd", typeof(vd),(/"time","level","lat","lon"/))

   
    filevarattdef(fout, "level", vd&level)
    filevarattdef(fout, "lat", vd&latitude)
    filevarattdef(fout, "lon", vd&longitude)
    filevarattdef(fout, "time", time)
    filevarattdef(fout, "vd", vd)
 
    
    setfileoption(fout,"DefineMode",False)


    fout->time=(/time/)
    fout->level=(/vd&level/)
    fout->lat=(/vd&latitude/)
    fout->lon=(/vd&longitude/)
    fout->vd = (/vd/)
    delete(vd)
   delete(time)

 

end 



