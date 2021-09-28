begin

    ;读入数据
    a=6371000
    g=9.8
    pi=get_pi("float")
 
    diri_output="/home/yangsong3/data-observation/linshh/obs-derived/msf/"

    model = "ERA5"
;;;;
 
    print("*******************start reading"+" "+"ERA5"+"**********************************************")
    yearstart = 197901
    yearend = 202012
    fils1 := "/home/yangsong3/data-observation/ERA5-monthly/pressure/w_mon_1x1_1979_2020.nc"
    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
 
    lat = f1->latitude
    lon = f1->longitude
    lev = f1->level
    levwant= (/1000.0, 925.0, 850.0, 700.0, 600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)   ;level反着取

    levind1 = get1Dindex(lev,levwant)

    level  = f1->level(levind1)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    ; w = short2flt(f1->w(timeind,levind1,::-1,:))
    w = short2flt(f1->w(timeind,{500},::-1,:))

    printVarSummary(w)

    if any(getVarDimNames(w).eq."longitude") then 
     w!1 = "lat"
     w!2 = "lon"
    end if


    print("****************计算fi矢量******************************")
    uvd=dv2uvF_Wrap(w)
    delete(w)
 
    ud=uvd(0,:,:,:)
    vd=uvd(1,:,:,:)
    delete(uvd)
    printVarSummary(ud)
    printVarSummary(vd)
    
    print("**********计算local 沃克和哈得来*************************")
    w_zonal = ud              
    T_grad_lat = ud             
    gradsf (ud, w_zonal, T_grad_lat)  

    fi=cos((w_zonal&lat/180)*pi)

    printVarSummary(w_zonal)
    delete(T_grad_lat)
   
    T_grad_lon = vd             
    w_meri = vd               
    gradsf (vd,T_grad_lon, w_meri)  
    printVarSummary(w_meri)
    delete(T_grad_lon)
   

    m_zonal=w_zonal
    m_meri=w_meri
  
    do i=0,dimsizes(w_zonal&lat)-1
        m_zonal(:,i,:)=(/-1*w_zonal(:,i,:)*fi(i)/g/)
        m_meri(:,i,:)=(/-1*w_meri(:,i,:)*fi(i)/g/)
    end do
   delete(w_zonal)
   delete(w_meri)

    ; m_zonal_mean = m_zonal(:,{500},:,:)
    ; m_meri_mean = m_meri(:,{500},:,:)
    
    printVarSummary(m_zonal)
 ;;;


;;;写nc文件********************************************************************

    dim=dimsizes(m_zonal)
    ntime=dim(0)
    ; nlev=dim(1)
    nlat=dim(1)
    nlon=dim(2)

;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

   
    filo1="m_zonal_meri_ERA5_500hpa.nc"
    system("/bin/rm -f "+diri_output + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diri_output + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "ERA5计算的local沃克和哈得来"    
    fAtt@creation_date = "2021.5.25"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(m_zonal&time),(/"time"/))
    ; filevardef(fout, "level",typeof(m_zonal&level), (/"level"/))
    filevardef(fout, "lat", typeof(m_zonal&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(m_zonal&lon),(/"lon"/))

    filevardef(fout, "m_zonal", typeof(m_zonal),(/"time","lat","lon"/))
    filevardef(fout, "m_meri", typeof(m_meri),(/"time","lat","lon"/))
    ; filevardef(fout, "fi_zonal", typeof(ud),(/"time","level","lat","lon"/))
    ; filevardef(fout, "fi_meri", typeof(vd),(/"time","level","lat","lon"/))  
 
  
    filevarattdef(fout, "time", m_zonal&time)   ;;;记得给时间附上时间单位
    ; filevarattdef(fout, "level", m_zonal&level)
    filevarattdef(fout, "lat", m_zonal&lat)
    filevarattdef(fout, "lon", m_zonal&lon)

    filevarattdef(fout, "m_zonal", m_zonal)
    filevarattdef(fout, "m_meri", m_meri)
    ; filevarattdef(fout, "fi_zonal", ud)
    ; filevarattdef(fout, "fi_meri", vd)    
 
   
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/m_zonal&time/)
    fout->lat=(/m_zonal&lat/)
    fout->lon=(/m_zonal&lon/)

    ; fout->m_zonal=(/m_zonal_mean/)
    ; fout->m_meri=(/m_meri_mean/)
    fout->m_zonal=(/m_zonal/)
    fout->m_meri=(/m_meri/)

    ; fout->fi_zonal=(/ud/)
    ; fout->fi_meri=(/vd/)

end 


