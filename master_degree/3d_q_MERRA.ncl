
begin
    ;;;读入MERRA的数据
    m=""+ispan(1979,2012,1)+""
    fils=new((/34,12/),string)
    do i=0,33
    if i.le.13 then
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA/MERRA100.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
    end if
    if i.ge.14.and.i.le.21 then 
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA/MERRA200.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
    end if

    if i.ge.22 then 
    fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA/MERRA300.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
    end if

    fils(i,:)= fils3
    end do
    delete(fils3)
    fils_new=onedtond(fils, 34*12)
    fileh=new(408,string)
    fileh(0:407)=fils_new

    f1=addfiles(fileh, "join")
    printVarSummary(f1)
    time=f1[:]->TIME
    lon=f1[1]->lon
    lat=f1[1]->lat
    date=cd_calendar(f1[:]->TIME, 1)
    level=f1[1]->Height
    timeind=ind(date.le.201212.and.date.ge.197901)
    levind=ind(level.le.1000.and.level.ge.100)
    print("levelind="+levind)
    w=f1[:]->OMEGA(timeind,levind,{5:-5},:)
    ; w=(/w*86400.0/)

    w!0="time"
    w!1="level"
    w!2="lat"
    w!3="lon"
    delete(f1)
    delete(fils)
    delete(date)
    ;;;

    fils=new((/34,12/),string)
    do i=0,33
  if i.le.13 then
  fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_air/MERRA100.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
  end if
  if i.ge.14.and.i.le.21 then 
  fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_air/MERRA200.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
  end if
 
  if i.ge.22 then 
  fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA_air/MERRA300.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
  end if
    fils(i,:)= fils3
    end do
    delete(fils3)
    fils_new=onedtond(fils, 34*12)
    fileh=new(408,string)
    fileh(0:407)=fils_new

    f1=addfiles(fileh, "join")
    printVarSummary(f1)
    date=cd_calendar(f1[:]->TIME, 1)
    timeind=ind(date.le.201212.and.date.ge.197901)
    t=f1[:]->T(timeind,levind,{5:-5},:)
 
    t!0="time"
    t!1="level"
    t!2="lat"
    t!3="lon"
    printVarSummary(t)

  p=level(levind)*100.0

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
 Q(:,i,:,:)=tofloat(((p(i)/100000.0)^0.286)*w(:,i,:,:)*dtdp(:,i,:,:)*3600.0*24.0)
 end do
 delete(w)
 copy_VarCoords(pt, Q)
 Q_time=dim_avg_n_Wrap(Q(time|:,lat|:,level|:,lon|:), 1)
 Q_time_new=Q_time(level|:,lon|:,time|:)
nlev=dimsizes(Q_time_new&level)
nlon=dimsizes(Q_time_new&lon)
ntime=dimsizes(Q_time_new&time)

;****************************************************************************
diri="/mnt/f/share/MERRA/"
;随时间变化的量
    diro=diri
    filo1="Q_time_MERRA.nc"
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
    filevarattdef(fout, "time", Q_time_new&time)
    filevarattdef(fout, "lon", Q_time_new&lon)

    setfileoption(fout,"DefineMode",False)

    fout->level=(/Q_time_new&level/)
    fout->lon=(/Q_time_new&lon/)
    fout->time=(/Q_time_new&time/)
    fout->Q_time=(/Q_time_new/)

 end