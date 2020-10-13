 begin
  ;数据从187101到201212
  ;精度为2度乘2度。
    diri="/mnt/f/share/MRE/"
    fils1=systemfunc("ls "+diri+"mre_u_v.nc")
    
    f1=addfile(fils1, "r")

    a=6371000   
    g=9.8
    pi=get_pi("double")
    dfi=(pi*10)/180
    time=f1->time
    level=f1->level
    lat=f1->lat
    lon=f1->lon


    date=cd_calendar(time, 1)
    
    tind1=ind(date.eq.195801)
    tind2=ind(date.eq.201212)
    time_wr=time(tind1:tind2)

    latind1=ind(lat.eq.60)
    latind2=ind(lat.eq.-60)
    lonind1=ind(lon.eq.0)
    lonind2=ind(lon.eq.357.5)
    levind1=ind(level.eq.100)
    levind2=ind(level.eq.1000)



    uwnd=f1->uwnd(tind1:tind2,:,::-1,:)
    vwnd=f1->vwnd(tind1:tind2,:,::-1,:)

    vrdvx=uv2vrdvF(uwnd,vwnd)
    uvd=dv2uvF(vrdvx(1,:,:,:,:))
    delete(vrdvx)
    ud=uvd(0,:,:,:,:)   ;;;这里为了得到经向的辐散风场
    vd=uvd(1,:,:,:,:)   ;;;这里为了得到经向的辐散风场
    delete(uvd)
    ud!0="time"
    ud!1="level"
    ud!2="lat"
    ud!3="lon"
    copy_VarCoords(ud, vd)
    printVarSummary(vd)
    ud_new=ud(level|:,lat|latind1:latind2,lon|lonind1:lonind2,time|:)  
    vd_new=vd(level|:,lat|latind1:latind2,lon|lonind1:lonind2,time|:)  
    delete(ud)
    delete(vd)
    
    ud_new&lon=lon(lonind1:lonind2)
    ud_new&lat=lat(latind1:latind2)
    ud_new&time= time_wr
    ud_new&level=level

    vd_new&lon=lon(lonind1:lonind2)
    vd_new&lat=lat(latind1:latind2)
    vd_new&time= time_wr
    vd_new&level=level


    ntime=dimsizes(ud_new&time)
    nlev=dimsizes(ud_new&level)
    nlon=dimsizes(ud_new&lon)
    nlat=dimsizes(ud_new&lat)

    ;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro=diri
    filo1="ud_vd_mre_global.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "ud的量,,4个再分析数据的结果，有时间变量1958-2012"    
    fAtt@creation_date = "2017.12.2"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lat","lon","time"/)
    dimsize=(/nlev,nlat,nlon,ntime/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    filevardef(fout, "level",typeof(level), (/"level"/))
    filevardef(fout, "lat", typeof(lat),(/"lat"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(time),(/"time"/))
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "ud_new", typeof(ud_new),(/"level","lat","lon","time"/))
    filevardef(fout, "vd_new", typeof(vd_new),(/"level","lat","lon","time"/))

    filevarattdef(fout, "level", level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "lat", lat)
    filevarattdef(fout, "time", time_wr)

    fout->level=(/ud_new&level/)
    fout->lon=(/ud_new&lon/)
    fout->lat=(/ud_new&lat/)
    fout->time=(/time_wr/)
    fout->ud_new=(/ud_new/)
    fout->vd_new=(/vd_new/)

    end