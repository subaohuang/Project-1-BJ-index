begin
    ; diri="/mnt/f/share/ncep/ncep1/"
    diri="/mnt/f/share/ncep/20CR/"
    diri_pi="/mnt/f/share/era/data/"
    a=6371000
    g=9.8
    fils1=systemfunc("ls "+diri+"omega.mon.mean.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.201212.and.date.ge.195801)
    w=f1->omega(timeind,{1000:100},::-1,:)
    delete(f1)
    print("计算fi矢量")
    uvd=dv2uvF_Wrap(w)
    delete(w)
 
    ud=uvd(0,:,:,:,:)
    vd=uvd(1,:,:,:,:)
    delete(uvd)
    printVarSummary(ud)
    printVarSummary(vd)
    
    print("**********计算local 沃克和哈得来*************************")
    w_zonal = ud              
    T_grad_lat = ud             
    gradsf (ud, w_zonal, T_grad_lat)  

    pi=get_pi("float")
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
    m_zonal(:,:,i,:)=(/-1*w_zonal(:,:,i,:)*fi(i)/g/)
    m_meri(:,:,i,:)=(/-1*w_meri(:,:,i,:)*fi(i)/g/)
   end do
   delete(w_zonal)
   delete(w_meri)

    m_zonal_mean=dim_avg_n_Wrap(m_zonal(:,{500},:,:), 0)
    m_meri_mean=dim_avg_n_Wrap(m_meri(:,{500},:,:), 0)
    
    printVarSummary(m_zonal)

 ;;;


;;;写nc文件********************************************************************

dim=dimsizes(m_zonal)
ntime=dim(0)
nlev=dim(1)
nlat=dim(2)
nlon=dim(3)

;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro=diri
    filo1="20cr_m_zonal_meri_new.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "20cr计算的local沃克和哈得来"    
    fAtt@creation_date = "2017.12.28"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(m_zonal&time),(/"time"/))
    filevardef(fout, "level",typeof(m_zonal&level), (/"level"/))
    filevardef(fout, "lat", typeof(m_zonal&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(m_zonal&lon),(/"lon"/))

    filevardef(fout, "m_zonal", typeof(m_zonal),(/"time","level","lat","lon"/))
    filevardef(fout, "m_meri", typeof(m_meri),(/"time","level","lat","lon"/))
    filevardef(fout, "fi_zonal", typeof(ud),(/"time","level","lat","lon"/))
    filevardef(fout, "fi_meri", typeof(vd),(/"time","level","lat","lon"/))  
 
  
    filevarattdef(fout, "time", m_zonal&time)   ;;;记得给时间附上时间单位
    filevarattdef(fout, "level", m_zonal&level)
    filevarattdef(fout, "lat", m_zonal&lat)
    filevarattdef(fout, "lon", m_zonal&lon)

    filevarattdef(fout, "m_zonal", m_zonal)
    filevarattdef(fout, "m_meri", m_meri)
    filevarattdef(fout, "fi_zonal", ud)
    filevarattdef(fout, "fi_meri", vd)    
 
   
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/m_zonal&time/)
    fout->level=(/m_zonal&level/)
    fout->lat=(/m_zonal&lat/)
    fout->lon=(/m_zonal&lon/)

    fout->m_zonal=(/m_zonal/)
    fout->m_meri=(/m_meri/)
    ; fout->fi_zonal=(/ud/)
    ; fout->fi_meri=(/vd/)

end 


