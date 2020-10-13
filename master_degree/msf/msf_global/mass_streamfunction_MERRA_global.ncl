  begin
  diri="/mnt/f/share/MERRA/MERRA/"
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


  fils_new=onedtond(fils, 34*12)
  fileh=new(408,string)
  fileh(0:407)=fils_new

  f1=addfiles(fileh, "join")
  printVarSummary(f1)
  date=cd_calendar(f1[:]->TIME, 1)
  
 
  timeind=ind(date.le.201212.and.date.ge.197901)
 
    a=6371000   
    g=9.8
    pi=get_pi("double")
    ; dfi=(pi*20)/180
    dfi=(pi*10)/180

    time=f1[:]->TIME
    level=f1[0]->Height
    lat=f1[0]->lat
  
    date=cd_calendar(time, 1)
    
    tind1=ind(date.eq.197901)
    tind2=ind(date.eq.201212)
    time_wr=time(tind1:tind2)

    latind1=ind(lat.eq.5)
    latind2=ind(lat.eq.-5)

    levind1=ind(level.eq.100)
    levind2=ind(level.eq.1000)

    
  
    uwnd_in = lonPivot(f1[:]->U(tind1:tind2,:,::-1,:),0)
  

    vwnd_in = lonPivot(f1[:]->V(tind1:tind2,:,::-1,:),0)
    printVarSummary(uwnd_in)

    ;;;interp 180 degree data
    xo = fspan(0, 357.5, 144)

    uwnd = linint1_n_Wrap(uwnd_in&lon, uwnd_in, True, xo, 0, 3)

    vwnd = linint1_n_Wrap(vwnd_in&lon, vwnd_in, True, xo, 0, 3)

    ;;;

    uwnd!3 = "lon"
    uwnd&lon = xo
    uwnd&lon@units = "degree_east"

    vwnd!3 = "lon"
    vwnd&lon = xo
    vwnd&lon@units = "degree_east"

    printVarSummary(uwnd)
 


 ;;;interp missing value
     print("****************interp missing value*********************")
       
    ;******************************************************************
    ; use poisson_grid_fill to fill missing value
    ;******************************************************************
    nscan = 1500           ; usually *much* fewer
    eps   = 0.001            ; variable depended
    gtype = False            ; "gendat" does not generate cyclic fields
    guess = 0                ; use zonal means
    relc  = 0.6              ; standard relaxation coef
    opt   = 0
    
    poisson_grid_fill(uwnd, gtype, guess, nscan, eps, relc, opt)
    
    poisson_grid_fill(vwnd, gtype, guess, nscan, eps, relc, opt)
;;;---------------------------------------------------

    w_r_in = lonPivot(f1[:]->OMEGA(tind1:tind2,levind2:levind1,latind2:latind1,:),0)

    w_r = linint1_n_Wrap(w_r_in&lon, w_r_in, True, xo, 0, 3)
     

    w_r!3 = "lon"
    w_r&lon = xo
    w_r&lon@units = "degree_east"


    le=w_r&Height
    ti=w_r&TIME
    lon=w_r&lon
    w_r!0="time"
    w_r!1="level"
    w_r&time=ti
    w_r&level=le
  
    printVarSummary(w_r)
    
     w=w_r
     ; w=w_r*86400.0
    ; print(w)

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
    ud_new=ud(time|:,lat|latind2:latind1,lon|:,level|:)  
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
    w_avg_new=w_avg
    copy_VarCoords(w_avg, w_avg_new)
   
    pbot=w_avg_new&level
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),double)
    
    printVarSummary(ud_new)

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
    filo1="msf_time_MERRA_global.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
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




