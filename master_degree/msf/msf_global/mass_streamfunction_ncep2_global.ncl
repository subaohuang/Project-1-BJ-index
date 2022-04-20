  begin
  ;数据从187101到201212
  ;精度为2度乘2度。
  diri="/mnt/f/share/ncep/ncep2/"
    fils1=systemfunc("ls "+diri+"uwnd.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"vwnd.mon.mean.nc")
    fils4=systemfunc("ls "+diri+"omega.mon.mean.nc")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f4=addfile(fils4, "r")
    a=6371000   
    g=9.8
    pi=get_pi("double")
    ; dfi=(pi*20)/180
    dfi=(pi*10)/180

    time=f1->time
    level=f1->level
    lat=f1->lat
    lon=f1->lon
    levelw=f4->level

    date=cd_calendar(time, 1)
    
    tind1=ind(date.eq.197901)
    tind2=ind(date.eq.201212)
    time_wr=time(tind1:tind2)
    latind1=ind(lat.eq.5)
    latind2=ind(lat.eq.-5)
    lonind1=ind(lon.eq.0)
    lonind2=ind(lon.eq.357.5)
    levind1=ind(level.eq.100)
    levind2=ind(level.eq.1000)

    levindw1=ind(levelw.eq.100)
    levindw2=ind(levelw.eq.1000)
  
    uwnd_r=f1->uwnd(tind1:tind2,:,::-1,:)
    vwnd_r=f2->vwnd(tind1:tind2,:,::-1,:)
    

    uwnd=dble2flt(uwnd_r*uwnd_r@scale_factor+uwnd_r@add_offset)
    vwnd=dble2flt(vwnd_r*vwnd_r@scale_factor+vwnd_r@add_offset)

    copy_VarCoords(uwnd_r, uwnd)
    copy_VarCoords(vwnd_r, vwnd)
    delete(uwnd_r)
    delete(vwnd_r)

    w_r=f4->omega(tind1:tind2,levindw2:levindw1,latind2:latind1,lonind1:lonind2)
    
    w=dble2flt(w_r*w_r@scale_factor+w_r@add_offset)
     copy_VarCoords(w_r, w)
     delete(w_r)

    w_new=w(time|:,lat|:,level|:,lon|:) 
    w_avg=dim_avg_n_Wrap(w_new, (/0,1/))
   
    w_time=dim_avg_n_Wrap(w_new, 1)
    w_time_t=w_time(level|:,lon|:,time|:)
    delete(w_time)
    w_time_new=w_time_t*(-50)
    copy_VarMeta(w_time_t, w_time_new)
    delete(w_time_t)

    vrdvx=uv2vrdvF(uwnd,vwnd)
    uvd=dv2uvF(vrdvx(1,:,:,:,:))
    delete(vrdvx)
    ud=uvd(0,:,:,:,:)
    delete(uvd)
    ud!0="time"
    ud!1="level"
    ud!2="lat"
    ud!3="lon"
    ud_new=ud(time|:,lat|latind2:latind1,lon|lonind1:lonind2,level|:)  
    delete(ud)
    ud_new&lon=w&lon
    ud_new&lat=fspan(-5,5,dimsizes(w&lat))
    ud_new&time=ispan(1,dimsizes(w&time),1)
    ud_new&level=level
    
    ud_avg=dim_avg_n_Wrap(ud_new,(/0,1/))
    ud_avg_new=ud_avg(level|0:dimsizes(w&level)-1,lon|:)  

  

    ud_time=dim_avg_n_Wrap(ud_new(time|:,lat|:,level|0:dimsizes(w&level)-1,lon|:), 1)
    
    ud_time_new=ud_time(level|:,lon|:,time|:)
   
    w_avg_new=w_avg*(-50)
    copy_VarCoords(w_avg, w_avg_new)

   
    pbot=w_avg_new&level
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    do i=0,(dimsizes(w_avg_new&level)-1)
    vint(i,:,:,:)=vibeta(level, ud_new, 1, pres_new,pbot(i),0)
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
 
  
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro=diri
    filo1="msf_time_ncep2_global.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量1950-2012，取5到-5平均的结果"    
    fAtt@creation_date = "2017.8.7"
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
    filevardef(fout, "msf", typeof(fi_mean),(/"level","lon"/))    
    filevardef(fout, "ud", typeof(ud_avg_new),(/"level","lon"/))
    filevardef(fout, "w", typeof(w_avg_new),(/"level","lon"/))
    filevardef(fout, "msf_time", typeof(fi_time_new),(/"level","lon","time"/))
    filevardef(fout, "ud_time", typeof(ud_time_new),(/"level","lon","time"/))
    filevardef(fout, "w_time", typeof(w_time_new),(/"level","lon","time"/))

    filevarattdef(fout, "level", level)
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




