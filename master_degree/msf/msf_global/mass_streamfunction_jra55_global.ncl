  begin
;;读入数据
  diri="/mnt/f/share/JRA-55/"
  fils4=systemfunc("ls /mnt/f/share/ncep/ncep1/uwnd.mon.mean.nc")
  f4=addfile(fils4, "r")
  time=f4->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
  data=cd_calendar(time, 1)
  tind1=ind(data.eq.195801)
  tind2=ind(data.eq.201212)
  time_wr=time(tind1:tind2)
  delete(f4)

  fils1=systemfunc("ls "+diri+"uwnd.jra55.mon.mean.nc")
  fils2=systemfunc("ls "+diri+"vwnd.jra55.mon.mean.nc")
        
  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")


  print("读入风场")
  z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)
  levind=get1Dindex(f1->level,z0)

  
  uwnd_r=f1->uwnd(0:659,levind,::-1,:)  ;;;195801-201212
  vwnd_r=f2->vwnd(0:659,levind,::-1,:)  ;;;
  print(uwnd_r&level)
   
  level=uwnd_r&level   ;level反着取
  print("level="+level)
  levelw=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
  levindw=get1Dindex(f1->level,levelw)
  delete(f1)
  delete(f2)
   xo=fspan(0,357.5,144)
   yo=fspan(-90, 90, 73)
 
  
   print("开始水平插值")  
   uwnd=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, \
                          True, xo, yo, 0)
  
   vwnd=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, \
                          True, xo, yo, 0)
    
   uwnd!2="lat"
   uwnd!3="lon"
   
   vwnd!2="lat"
   vwnd!3="lon"

   delete(uwnd_r)
   delete(vwnd_r)
   
   printVarSummary(uwnd)
   printVarSummary(vwnd)
    
    lat=uwnd&lat
    lon=uwnd&lon
     
    lon@units="degree_east"


    latind1=ind(uwnd&lat.eq.5)
    latind2=ind(uwnd&lat.eq.-5)

    lonind1=ind(uwnd&lon.eq.0)
    lonind2=ind(uwnd&lon.eq.357.5)



   print("开始计算辐散风")  
    vrdvx=uv2vrdvF(uwnd,vwnd)
    uvd=dv2uvF(vrdvx(1,:,:,:,:))
    delete(vrdvx)
    ud=uvd(0,:,:,:,:)
    delete(uvd)
    ud!0="time"
    ud!1="level"
    ud!2="lat"
    ud!3="lon"
    copy_VarCoords(uwnd, ud)
    delete(uwnd)
    delete(vwnd)

    ud_new=ud(time|:,lat|latind2:latind1,lon|lonind1:lonind2,level|:)  
    delete(ud)

    fils3=systemfunc("ls "+diri+"omega.jra55.mon.mean.nc")
    f3=addfile(fils3, "r")
  
    w_r=f3->omega(0:659,levindw ,{5:-5},:)
    delete(f3)
    
    print("开始垂直速度的水平插值")  
     w=linint1_n_Wrap(w_r&lon,w_r,True, xo,0, 3) 
     w!2="lat"
     w!3="lon"
     delete(w_r)

    w_new=w(time|:,lat|:,level|:,lon|:) 
    delete(w)
    w_avg=dim_avg_n_Wrap(w_new, (/0,1/))
   
    w_time=dim_avg_n_Wrap(w_new, 1)
    w_time_t=w_time(level|:,lon|:,time|:)
    delete(w_time)
    w_time_new=w_time_t*(-50)
    copy_VarMeta(w_time_t, w_time_new)
    delete(w_time_t)
    w_avg_new=w_avg*(-50)
    copy_VarCoords(w_avg, w_avg_new)
    delete(w_avg)
    
    ud_avg=dim_avg_n_Wrap(ud_new,(/0,1/))
    ud_avg_new=ud_avg(level|0:dimsizes(levelw)-1,lon|:)  
    
    delete(ud_avg)

    ud_time=dim_avg_n_Wrap(ud_new(time|:,lat|:,level|0:dimsizes(levelw)-1,lon|:), 1)
    
    ud_time_new=ud_time(level|:,lon|:,time|:)
    delete(ud_time)
    


    a=6371000
    g=9.8
    pi=get_pi("double")
    dfi=(pi*10)/180
   print("开始积分")
    pbot=levelw
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(levelw),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    do i=0,(dimsizes(levelw)-1)
    vint(i,:,:,:)=vibeta(ud_new&level, ud_new(:,:,:,:), 1, pres_new,pbot(i),0)
    end do
    fi=((a*dfi)/g)*vint*100
    delete(vint)

    fi_mean=dim_avg_n_Wrap(fi, (/1,2/))
     
    fi_mean!0="level"
    fi_mean!1="lon"
    fi_mean&level=pbot
    fi_mean&level@units="hpa"
    fi_mean&lon=ud_new&lon
    fi_time=dim_avg_n_Wrap(fi, 2)
    delete(fi)

    fi_time!0="level"
    fi_time!1="time"
    fi_time!2="lon"
    fi_time&level=pbot
    fi_time&level@units="hpa"
    fi_time&lon=ud_new&lon
    fi_time&time=ud_new&time
    fi_time_new=fi_time(level|:,lon|:,time|:)
    
    ntime=dimsizes(fi_time_new&time)
    nlev=dimsizes(fi_time_new&level)
    nlon=dimsizes(fi_time_new&lon)
 
  print("开始写nc文件")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro=diri
    filo1="msf_time_jra55_global.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量1958-2012，取5到-5平均的结果"    
    fAtt@creation_date = "2017.12.13"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(fi_time_new&level), (/"level"/))
    filevardef(fout, "lon", typeof(fi_time_new&lon),(/"lon"/))
    filevardef(fout, "time", typeof(time_wr),(/"time"/))
    filevardef(fout, "msf", typeof(fi_mean),(/"level","lon"/))    
    filevardef(fout, "ud", typeof(ud_avg_new),(/"level","lon"/))
    filevardef(fout, "w", typeof(w_avg_new),(/"level","lon"/))
    filevardef(fout, "msf_time", typeof(fi_time_new),(/"level","lon","time"/))
    filevardef(fout, "ud_time", typeof(ud_time_new),(/"level","lon","time"/))
    filevardef(fout, "w_time", typeof(w_time_new),(/"level","lon","time"/))

    filevarattdef(fout, "level", fi_time_new&level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "time", time_wr)
    filevarattdef(fout, "ud", ud_avg_new)
    filevarattdef(fout, "w", w_avg_new)
    filevarattdef(fout, "msf_time", fi_time_new)
    filevarattdef(fout, "ud_time", ud_time_new)
    filevarattdef(fout, "w_time", w_time_new)
 
   
    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_mean&level/)
    fout->lon=(/fi_mean&lon/)
    fout->time=(/time_wr/)
    fout->ud=(/ud_avg_new/)
    fout->w=(/w_avg_new/)
    fout->msf=(/fi_mean/)
    fout->msf_time=(/fi_time_new/)
    fout->ud_time=(/ud_time_new/)
    fout->w_time=(/w_time_new/)

end 



