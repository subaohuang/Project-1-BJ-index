  begin
;;读入数据


;;;;
    model = "20CR"
    diri_output="/home/yangsong3/data-observation/linshh/obs-derived/vd/"
    
    yearstart = 187101
    yearend = 201212

    fils1 := systemfunc("ls /home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc")
    fils2 := systemfunc("ls /home/yangsong3/data-observation/20thCRV2/vwnd.mon.mean.nc")

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    time=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)


    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    time_r = time(timeind)
    print(date(timeind))

    print("读入风场")
    plev = f1->level
    
    print(plev)
    exit()

    z0=(/1000,925,850,700,600,500,400,300,250,150,100/)*1.0

    levind=get1Dindex(plev,z0)

    print(levind)
    exit()
  

    uwnd := f1->uwnd(timeind,levind,:,:)  ;;;195801-201212
    vwnd := f2->vwnd(timeind,levind,:,:)  ;;;
    

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
    filo1="vd_"+model+"_187101-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "20CR Meridional Divergence wind"    
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
    filevardef(fout, "lat", typeof(vd&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(vd&lon),(/"lon"/))

   
    filevardef(fout, "vd", typeof(vd),(/"time","level","lat","lon"/))

   
    filevarattdef(fout, "level", vd&level)
    filevarattdef(fout, "lat", vd&lat)
    filevarattdef(fout, "lon", vd&lon)
    filevarattdef(fout, "time", time_r)
    filevarattdef(fout, "vd", vd)
 
    
    setfileoption(fout,"DefineMode",False)


    fout->time=(/time_r/)
    fout->level=(/vd&level/)
    fout->lat=(/vd&lat/)
    fout->lon=(/vd&lon/)
    fout->vd = (/vd/)
    delete(vd)
   delete(time_r)
   delete(time)
 

end 



