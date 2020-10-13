  begin
;;读入数据
   
    ; model="ACCESS1-0"
    model="CanESM2"
    model="GFDL-CM3"
    diri_output="/mnt/f/share/CMIP5/msf_history_natureforce/"
   ; system("if ! test -d " + diri_output +" ; then mkdir " + diri_output + " ; fi")   ;;;make dir
   ; diri_output="/mnt/f/share/CMIP5/msf_history/"+model
   yearstart=195001
   ; yearend=201212
   yearend=200512
   ; fils1="/mnt/j/CMIP5/history/ua/"+model+"/r1i1p1/ua_Amon_"+model+"_historical_r1i1p1_"+yearstart+"-"+yearecnd+".nc"
   ; fils2="/mnt/j/CMIP5/history/va/"+model+"/r1i1p1/va_Amon_"+model+"_historical_r1i1p1_"+yearstart+"-"+yearend+".nc"
   
    ; fils1="/mnt/j/CMIP5/history_natureforce/ua/"+model+"/r1i1p1/ua_Amon_"+model+"_historicalNat_r1i1p1_"+"185001-201212"+".nc"
    ; fils2="/mnt/j/CMIP5/history_natureforce/va/"+model+"/r1i1p1/va_Amon_"+model+"_historicalNat_r1i1p1_"+"185001-201212"+".nc"
   
    fils1=systemfunc("ls /mnt/j/CMIP5/history_natureforce/ua/"+model+"/r1i1p1/ua_Amon_"+model+"_historicalNat_r1i1p1"+"*"+".nc")
    fils2=systemfunc("ls /mnt/j/CMIP5/history_natureforce/va/"+model+"/r1i1p1/va_Amon_"+model+"_historicalNat_r1i1p1"+"*"+".nc")
    print(fils1)
    print(fils2)
    

    f1=addfiles(fils1, "join")
    f2=addfiles(fils2, "join")


    time=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date=cd_calendar(time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    print("读入风场")
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,150.0,100.0/)*100.0
    levind=get1Dindex(f1[0]->plev,z0)

    uwnd_r=f1[:]->ua(timeind,levind,:,:)  ;;;195801-201212
    vwnd_r=f2[:]->va(timeind,levind,:,:)  ;;;
    uwnd_r!1="level"
    vwnd_r!1="level"
    uwnd_r&level=(/uwnd_r&level/100.0/)
    vwnd_r&level=(/vwnd_r&level/100.0/)

    print(uwnd_r&level)
    delete(f1)
    delete(f2)

    level=uwnd_r&level   ;level反着取
    levelw=level

    xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)

     print("*****************开始水平插值************************")  
    uwnd=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, \
                      True, xo, yo, 0)

    vwnd=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, \
                      True, xo, yo, 0)
     uwnd!2="lat"
     uwnd!3="lon"
     
     vwnd!2="lat"
     vwnd!3="lon"
   
   ; uwnd=uwnd_r
   ; vwnd=vwnd_r

   uwnd&lat@units="degrees_north"
   uwnd&lon@units="degrees_east"

   vwnd&lat@units="degrees_north"
   vwnd&lon@units="degrees_east"

   delete(uwnd_r)
   delete(vwnd_r)
   
   printVarSummary(uwnd)
   printVarSummary(vwnd)
    
    lat=uwnd&lat
    lon=uwnd&lon
     
    lon@units="degree_east"
    latind1=ind(uwnd&lat.eq.5)
    latind2=ind(uwnd&lat.eq.-5)

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
      
   print("*********************开始计算辐散风******************************")  
    vrdvx=uv2vrdvF(uwnd,vwnd)
    ; vrdvx=uv2dv_cfd(uwnd, vwnd, uwnd&lat, vwnd&lon, 0)
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

    ud_new=ud(time|:,lat|latind2:latind1,lon|:,level|:)  
    delete(ud)
    printVarSummary(ud_new)

    
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
    pbot=dble2flt(ud_new&level)
    pres_new=new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
    pres_new=1013.0
    vint=new((/dimsizes(ud_new&level),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
  
    do i=0,(dimsizes(ud_new&level)-1)
    vint(i,:,:,:)=dble2flt(vibeta(ud_new&level, ud_new(:,:,:,:), 1, pres_new,pbot(i),0))
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
    delete(fi_time)
    ntime=dimsizes(fi_time_new&time)
    nlev=dimsizes(fi_time_new&level)
    nlon=dimsizes(fi_time_new&lon)
  
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="msf_time_"+model+"_historyNat_r1i1p1_"+yearstart+"-"+yearend+".nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量1950-200512，取5到-5平均的结果"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(fi_time_new&level), (/"level"/))
    filevardef(fout, "lon", typeof(fi_time_new&lon),(/"lon"/))
    filevardef(fout, "time", typeof(fi_time_new&time),(/"time"/))

    filevardef(fout, "msf", typeof(fi_mean),(/"level","lon"/))    
    filevardef(fout, "ud", typeof(ud_avg_new),(/"level","lon"/))
    ; filevardef(fout, "w", typeof(w_avg_new),(/"level","lon"/))
    filevardef(fout, "msf_time", typeof(fi_time_new),(/"level","lon","time"/))
    filevardef(fout, "ud_time", typeof(ud_time_new),(/"level","lon","time"/))
    ; filevardef(fout, "w_time", typeof(w_time_new),(/"level","lon","time"/))

    filevarattdef(fout, "level", fi_time_new&level)
    filevarattdef(fout, "lon", fi_time_new&lon)
    filevarattdef(fout, "time", fi_time_new&time)
    filevarattdef(fout, "ud", ud_avg_new)
    ; filevarattdef(fout, "w", w_avg_new)
    filevarattdef(fout, "msf_time", fi_time_new)
    filevarattdef(fout, "ud_time", ud_time_new)
    ; filevarattdef(fout, "w_time", w_time_new)
 
   
    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_time_new&level/)
    fout->lon=(/fi_time_new&lon/)
    fout->time=(/fi_time_new&time/)
    fout->ud=(/ud_avg_new/)
    ; fout->w=(/w_avg_new/)
    fout->msf=(/fi_mean/)
    fout->msf_time=(/fi_time_new/)
    fout->ud_time=(/ud_time_new/)
    ; fout->w_time=(/w_time_new/)

end 



