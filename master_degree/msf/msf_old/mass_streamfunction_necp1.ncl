  begin
  diri="/media/sf_share/ncep/ncep1/"
    fils1=systemfunc("ls "+diri+"uwnd.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"vwnd.mon.mean.nc")
    fils3=systemfunc("ls "+diri+"pres.mon.mean.nc")
    fils4=systemfunc("ls "+diri+"omega.mon.mean.nc")
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f3=addfile(fils3, "r")
    f4=addfile(fils4, "r")
   
    time=f1->time
    lat=f1->lat
    lon=f1->lon
    level=f1->level
    data=cd_calendar(time, 1)
  
    tind1=ind(data.eq.197901)
    tind2=ind(data.eq.201701)
    latind1=ind(lat.eq.5)
    latind2=ind(lat.eq.-5)
    lonind1=ind(lon.eq.100)
    lonind2=ind(lon.eq.290)


    uwnd=f1->uwnd(tind1:tind2,:,::-1,:)
    vwnd=f2->vwnd(tind1:tind2,:,::-1,:)
    pres=f3->pres(tind1:tind2,latind2:latind1,lonind1:lonind2)
   
    w=f4->omega(tind1:tind2,:,latind2:latind1,lonind1:lonind2)

    w_new=w(time|:,lat|:,level|:,lon|:) 
    w_avg=dim_avg_n_Wrap(w_new, (/0,1/))
   



    a=6371000
    g=9.8
    pi=get_pi("double")
    dfi=(pi*10)/180

    vrdvx=uv2vrdvF(uwnd,vwnd)
    uvd=dv2uvF(vrdvx(1,:,:,:,:))
    delete(vrdvx)
    ud=uvd(0,:,:,:,:)
    delete(uvd)

    ud!0="time"
    ud!1="level"
    ud!2="lat"
    ud!3="lon"
    ud_new=ud(time|:,lat|latind1:latind2,lon|lonind1:lonind2,level|:)  
    ud_new&lon=pres&lon
    ud_new&lat=fspan(-5,5,5)
    ud_new&time=ispan(1,457,1)
    ud_new&level=uwnd&level
    ud_avg=dim_avg_n_Wrap(ud_new,(/0,1/))   
    ud_avg_new=ud_avg(level|0:dimsizes(w_avg&level)-1,lon|:)

    ud_time=dim_avg_n_Wrap(ud_new(time|:,lat|:,level|0:dimsizes(w_avg&level)-1,lon|:), 1)
    
    ud_time_new=ud_time(level|:,lon|:,time|:)
   
    w_avg_new=w_avg*(-50)
    copy_VarCoords(w_avg, w_avg_new)
     
    pbot=w_avg_new&level
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)

     do i=0,dimsizes(w_avg_new&level)-1
    vint(i,:,:,:)=vibeta(level, ud_new, 1, pres_new,pbot(i),0)
     end do
    fi=((a*dfi)/g)*vint*100
    
    fi_mean=dim_avg_n_Wrap(fi, (/1,2/))
    fi_mean!0="level"
    fi_mean!1="lon"
    fi_mean&level=pbot
    fi_mean&level@units="hpa"
    fi_mean&lon=ud_new&lon
    wks=gsn_open_wks("png", ""+diri+"mass_stream")

    nlev=dimsizes(fi_mean&level)
    nlon=dimsizes(fi_mean&lon)
 
    fi_time=dim_avg_n_Wrap(fi, 2)
    fi_time!0="level"
    fi_time!1="time"
    fi_time!2="lon"
    fi_time&level=pbot
    fi_time&level@units="hpa"
    fi_time&lon=ud_new&lon
    fi_time&time=ud_new&time
    fi_time_new=fi_time(level|:,lon|:,time|:)
    
    ntime=dimsizes(fi_time_new&time)
    
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

   diro=diri
    filo1="msf_time_necp1.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量，取5到-5平均的结果"    
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
  

    filevarattdef(fout, "level", level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "ud", ud_avg_new)
    filevarattdef(fout, "w", w_avg_new)
    filevarattdef(fout, "msf_time", fi_time_new)
    filevarattdef(fout, "ud_time", ud_time_new)
    
 
   
    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_mean&level/)
    fout->lon=(/fi_mean&lon/)
    fout->time=(/fi_time_new&time/)
    fout->ud=(/ud_avg_new/)
    fout->w=(/w_avg_new/)
    fout->msf=(/fi_mean/)
    fout->msf_time=(/fi_time_new/)
    fout->ud_time=(/ud_time_new/)
   
end 


