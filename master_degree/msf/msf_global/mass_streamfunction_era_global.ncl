  begin
;;读入数据
  diri="/mnt/f/share/era/data/data_series/"
  fils1=systemfunc("ls "+diri+"era_pre_global.nc")
  f1=addfile(fils1, "r")
 
  a=6371000
  g=9.8
  pi=get_pi("double")
  dfi=(pi*10)/180

  time=f1->time
  lat=f1->latitude
  lon=f1->longitude
  level=f1->level(::-1)   ;level反着取

  levelw=f1->level

  data=cd_calendar(time, 1)

  tind1=ind(data.eq.197901)
  tind2=ind(data.eq.201701)
   
  timew=time(tind1:tind2)

  latind1=ind(lat.eq.5)
  latind2=ind(lat.eq.-5)

  lonind1=ind(lon.eq.0)
  lonind2=ind(lon.eq.357.5)

  levind1=ind(level.eq.100)
  levind2=ind(level.eq.1000)

  levindw1=ind(levelw.eq.100)
  levindw2=ind(levelw.eq.1000)

  d=f1->d(tind1:tind2,::-1,::-1,:)
  omega=f1->w(tind1:tind2,levindw2:levindw1,latind2:latind1,lonind1:lonind2)
  printVarSummary(omega)
  div=dble2flt(d*d@scale_factor+d@add_offset)
  w=dble2flt(omega*omega@scale_factor+omega@add_offset)

  copy_VarMeta(omega, w)
  copy_VarMeta(d, div)
  delete(d)


  w_new=w(time|:,latitude|:,level|:,longitude|:) 
  delete(w)
  w_avg=dim_avg_n_Wrap(w_new, (/0,1/))

  w_time=dim_avg_n_Wrap(w_new, 1)
  w_time_t=w_time(level|:,longitude|:,time|:)
  w_time_new=w_time_t*(-50)
  copy_VarMeta(w_time_t, w_time_new)
  delete(w_time_t)
  delete(w_time)
  delete(w_new)

 ;;计算辐散风
    uvd=dv2uvF_Wrap(div)
    ud=uvd(0,:,:,:,:)
    delete(uvd)
    printVarSummary(ud)
    ud_new=ud(time|:,latitude|latind2:latind1,longitude|lonind1:lonind2,level|:)  
    delete(ud)
    ud_new&time=ispan(1,dimsizes(ud_new&time),1)
     

    ud_avg=dim_avg_n_Wrap(ud_new,(/0,1/)) 
    print(ud_avg&level) 
    ud_avg_new=ud_avg(level|levind2:levind1,longitude|:) 

    ud_time=dim_avg_n_Wrap(ud_new(time|:,latitude|:,level|levind2:levind1,longitude|:), 1)

    ud_time_new=ud_time(level|:,longitude|:,time|:)

    w_avg_new=w_avg*(-50)
    copy_VarCoords(w_avg, w_avg_new)

    pbot=w_avg_new&level
  
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&latitude),dimsizes(ud_new&longitude)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&latitude),dimsizes(ud_new&longitude)/),float)
     
     do i=0,dimsizes(w_avg_new&level)-1
    vint(i,:,:,:)=vibeta(level(:), ud_new, 1, pres_new,pbot(i),0)
     end do
   
    fi=((a*dfi)/g)*vint*100
 
    fi_mean=dim_avg_n_Wrap(fi, (/1,2/))
    fi_mean!0="level"
    fi_mean!1="lon"
    fi_mean&level=pbot
    fi_mean&level@units="hpa"
    fi_mean&lon=ud_new&longitude

    
    nlev=dimsizes(fi_mean&level)
    nlon=dimsizes(fi_mean&lon)
 
    fi_time=dim_avg_n_Wrap(fi, 2)
    fi_time!0="level"
    fi_time!1="time"
    fi_time!2="lon"
    fi_time&level=pbot
    fi_time&level@units="hpa"
    fi_time&lon=ud_new&longitude
    fi_time&time=ud_new&time
    fi_time_new=fi_time(level|:,lon|:,time|:)
    
    ntime=dimsizes(fi_time_new&time)
    nlev=dimsizes(fi_time_new&level)
    nlon=dimsizes(fi_time_new&lon)

;****************************************************************************
;随时间变化的量
    diro=diri
    filo1="msf_time_era_global.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量，取5到-5平均的结果,经度范围取的是全球"    
    fAtt@creation_date = "2017.10.30"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "level",typeof(level), (/"level"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(timew),(/"time"/))
    filevardef(fout, "msf", typeof(fi_mean),(/"level","lon"/))    
    filevardef(fout, "ud", typeof(ud_avg_new),(/"level","lon"/))
    filevardef(fout, "w", typeof(w_avg_new),(/"level","lon"/))
    filevardef(fout, "msf_time", typeof(fi_time_new),(/"level","lon","time"/))
    filevardef(fout, "ud_time", typeof(ud_time_new),(/"level","lon","time"/))
    filevardef(fout, "w_time", typeof(w_time_new),(/"level","lon","time"/))

    filevarattdef(fout, "level", level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "time", timew)
    filevarattdef(fout, "ud", ud_avg_new)
    filevarattdef(fout, "w", w_avg_new)
    filevarattdef(fout, "msf_time", fi_time_new)
    filevarattdef(fout, "ud_time", ud_time_new)
    filevarattdef(fout, "w_time", w_time_new)
    
 
   
    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_mean&level/)
    fout->lon=(/fi_mean&lon/)
    fout->time=(/timew/)
    fout->ud=(/ud_avg_new/)
    fout->w=(/w_avg_new/)
    fout->msf=(/fi_mean/)
    fout->msf_time=(/fi_time_new/)
    fout->ud_time=(/ud_time_new/)
    fout->w_time=(/w_time_new/)

   
end 


