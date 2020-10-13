
begin
    fils1="/mnt/f/share/JRA-55/omega.jra55.mon.mean.nc"
    fils2="/mnt/f/share/ncep/20CR/uwnd.10m.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f2->time, 1)

    timeind=ind(date.le.201212.and.date.ge.197901)
    time=f2->time(timeind)
    printVarSummary(time)

    level=f1->level
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind=get1Dindex(f1->level,z0)
      w_r=f1->omega(252:659,levind,{-5:5},:)
    if isatt(w_r, "add_offset") then 
      w=w_r@scale_factor*w_r + w_r@add_offset
    else 
      w=w_r
    end if
    copy_VarCoords(w_r, w)
    copy_VarAtts(w_r, w)
    w!0="time"
    w!1="level"
    w!2="lat"
    w!3="lon"
    delete(w&time)
    w&time=time
    w&time@units=time@units
    printVarSummary(w)
   

    delete(f1)
    delete(fils1)
    delete(levind)
    delete(level)

;;读入T
    fils1="/mnt/f/share/JRA-55/air.jra55.mon.mean.nc"
    f1=addfile(fils1, "r")
    level=f1->level
    lat=f1->lat
    lon=f1->lon
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    levind=get1Dindex(f1->level,z0)
    t_r=f1->air(252:659,levind,{-5:5},:)

    if isatt(t_r, "add_offset") then 
    t=t_r@scale_factor*t_r + t_r@add_offset
    else 
    t=t_r
    end if
    copy_VarCoords(t_r, t)
    copy_VarAtts(t_r, t)
    t!0="time"
    t!1="level"
    t!2="lat"
    t!3="lon"
    delete(t&time)
    t&time=date(timeind)
    printVarSummary(t)

 p=z0*100.0

 pt=pot_temp(p, t, 1, False)
 printVarSummary(pt)
 delete(t)
 printMinMax(pt, 1)                                        
 dtdp=center_finite_diff_n(pt, p,False , 0, 1)                                                         ;计算稳定度
 copy_VarCoords(pt, dtdp)
 printMinMax(w, 1)
 printVarSummary(dtdp)
 Q=new((/dimsizes(pt&time),dimsizes(pt&level),dimsizes(pt&lat),dimsizes(pt&lon)/),float)
 do i=0,dimsizes(p)-1
 Q(:,i,:,:)=((p(i)/100000.0)^0.286)*w(:,i,:,:)*dtdp(:,i,:,:)*3600.0*24.0
 end do
 delete(w)
 copy_VarCoords(pt, Q)
 Q_time=dim_avg_n_Wrap(Q(time|:,lat|:,level|:,lon|:), 1)
 Q_time_new=Q_time(level|:,lon|:,time|:)
nlev=dimsizes(Q_time_new&level)
nlon=dimsizes(Q_time_new&lon)
ntime=dimsizes(Q_time_new&time)

;****************************************************************************
;随时间变化的量
    diro="/mnt/f/share/JRA-55/"
    filo1="Q_time_JRA55.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "Q的值，取5到-5平均了"    
    fAtt@creation_date = "2017.8.25"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(level), (/"level"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(time),(/"time"/))
    filevardef(fout, "Q_time", typeof(Q_time_new),(/"level","lon","time"/))    
 
    
  
    filevarattdef(fout, "level", Q_time_new&level)
    filevarattdef(fout, "time", time)
    filevarattdef(fout, "lon", Q_time_new&lon)

    setfileoption(fout,"DefineMode",False)

    fout->level=(/Q_time_new&level/)
    fout->lon=(/Q_time_new&lon/)
    fout->time=(/time/)
    fout->Q_time=(/Q_time_new/)

 end