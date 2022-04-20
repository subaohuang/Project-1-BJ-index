  begin
  ;数据从187101到201212
  ;精度为2度乘2度。
  diri="/mnt/f/share/era/era40/"
  diri_in="/mnt/f/share/era/data/data_series/"
    fils1=systemfunc("ls "+diri+"era_40_u_v.nc")
    fils2=systemfunc("ls "+diri+"era_40_w.nc")
    fils3=systemfunc("ls "+diri_in+"era_pre_global.nc")
    
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f3=addfile(fils3, "r")
    a=6371000   
    g=9.8
    pi=get_pi("double")
    ; dfi=(pi*20)/180
    dfi=(pi*10)/180

    time=f1->time
    
    lat=f1->latitude
    lon=f1->longitude
   
    
    date_in=cd_calendar(f3->time, 1)
    tind_in=ind(date_in.le.201212.and.date_in.ge.200201)

    date=cd_calendar(time, 1)
    tind1=ind(date.eq.195801)
    tind2=ind(date.eq.200112)

    time_wr=new(660,typeof(time))

    time_wr(0:527)=time(tind1:tind2)
    time_wr(528:659)=f3->time(tind_in)
    

    latind1=ind(lat.eq.5)
    latind2=ind(lat.eq.-5)
    lonind1=ind(lon.eq.0)
    lonind2=ind(lon.eq.357.5)
    levind1=ind(f1->level.eq.100)
    levind2=ind(f1->level.eq.1000)
    
    level=f1->level(levind2:levind1)   ;;从下往上取
    levelw=level
    uwnd_40_r=f1->u(tind1:tind2,{1000:100},::-1,:)
    vwnd_40_r=f1->v(tind1:tind2,{1000:100},::-1,:)
    
    w_40_r=f2->w(tind1:tind2,{1000:100},{5:-5},:)

    uwnd_40=dble2flt(uwnd_40_r*uwnd_40_r@scale_factor+uwnd_40_r@add_offset)
    vwnd_40=dble2flt(vwnd_40_r*vwnd_40_r@scale_factor+vwnd_40_r@add_offset)
    w_40=dble2flt(w_40_r*w_40_r@scale_factor+w_40_r@add_offset)
    
    copy_VarCoords(uwnd_40_r, uwnd_40)
    copy_VarCoords(vwnd_40_r, vwnd_40)
    copy_VarCoords(w_40_r, w_40)
    delete(uwnd_40_r)
    delete(vwnd_40_r)
    delete(w_40_r)

    ;;;读入era_in的数据
    levind=get1Dindex(f3->level,uwnd_40&level)    ;;;这里容易出现错误，请注意

    uwnd_in_r=f3->u(tind_in,levind,::-1,:)
    vwnd_in_r=f3->v(tind_in,levind,::-1,:)
    
    w_in_r=f3->w(tind_in,levind,{5:-5},:)

    uwnd_in=dble2flt(uwnd_in_r*uwnd_in_r@scale_factor+uwnd_in_r@add_offset)
    vwnd_in=dble2flt(vwnd_in_r*vwnd_in_r@scale_factor+vwnd_in_r@add_offset)
    w_in=dble2flt(w_in_r*w_in_r@scale_factor+w_in_r@add_offset)

    copy_VarCoords(uwnd_in_r, uwnd_in)
    copy_VarCoords(vwnd_in_r, vwnd_in)
    copy_VarCoords(w_in_r, w_in)

    delete(uwnd_in_r)
    delete(vwnd_in_r)
    delete(w_in_r)
 

    uwnd=new((/dimsizes(uwnd_40&time)+dimsizes(uwnd_in&time), \
        dimsizes(uwnd_40&level),dimsizes(uwnd_40&latitude),dimsizes(uwnd_40&longitude)/),typeof(uwnd_40))
    
    w=new((/dimsizes(uwnd_40&time)+dimsizes(uwnd_in&time), \
        dimsizes(uwnd_40&level),dimsizes(w_40&latitude),dimsizes(w_40&longitude)/),typeof(w_40))
    uwnd!0="time"
    uwnd!1="level"
    uwnd!2="lat"
    uwnd!3="lon"
    
    w!0="time"
   
    uwnd&time=ispan(1,(/dimsizes(uwnd_40&time)+dimsizes(uwnd_in&time)/),1)
    copy_VarCoords(uwnd_40(1,:,:,:), uwnd(1,:,:,:))
    
    vwnd=uwnd
    
    w&time=ispan(1,(/dimsizes(uwnd_40&time)+dimsizes(uwnd_in&time)/),1)
    
    copy_VarCoords(w_40(1,:,:,:), w(1,:,:,:))
    w!1="level"
    w!2="lat"
    w!3="lon"

    printVarSummary(uwnd_in)
    printVarSummary(uwnd_40)
 
    uwnd(0:dimsizes(uwnd_40&time)-1,:,:,:)=(/uwnd_40/)
    uwnd(dimsizes(uwnd_40&time):659,:,:,:)=(/uwnd_in/)

    vwnd(0:dimsizes(vwnd_40&time)-1,:,:,:)=(/vwnd_40/)
    vwnd(dimsizes(vwnd_40&time):659,:,:,:)=(/vwnd_in/)
   

    printVarSummary(w_40)
    printVarSummary(w_in)
    printVarSummary(w)
    w(0:dimsizes(w_40&time)-1,:,:,:)=(/w_40/)
    w(dimsizes(w_40&time):659,:,:,:)=(/w_in/)

  
    printVarSummary(uwnd)
    printVarSummary(vwnd)
    printVarSummary(w)
 
    delete(uwnd_40)
    delete(vwnd_40)
    delete(w_40)
    delete(w_in)


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
    ud_avg_new=ud_avg(level|0:dimsizes(levelw)-1,lon|:)  

  

    ud_time=dim_avg_n_Wrap(ud_new(time|:,lat|:,level|0:dimsizes(levelw)-1,lon|:), 1)
    
    ud_time_new=ud_time(level|:,lon|:,time|:)
   
    w_avg_new=w_avg*(-50)
    copy_VarCoords(w_avg, w_avg_new)

   
    pbot=w_avg_new&level
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    do i=0,(dimsizes(w_avg_new&level)-1)
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
 
  
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro=diri
    filo1="msf_time_era40_in_global.nc"
    system("/bin/rm -f "+diro + filo1+"")
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




