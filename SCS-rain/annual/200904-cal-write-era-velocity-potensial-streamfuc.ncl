
;;********************************************************
;;计算速度势函数以及辐散风场
;; ERA-interim:u,v，
;; used ncl function: uv2dv_cfd,
;；；grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2020.09.04
;**********************************************************
begin
    
;;;time and region
 timestart = 197901
 timeend = 201612
;;;


;;;;read data
    diri="/home/yangsong3/data-observation/ERA-interim/monthly/"
    fils1=systemfunc("ls "+diri+"uwnd.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"vwnd.mon.mean.nc")

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    date2=cd_calendar(f2->time, 1)

    timeind=ind(date.le.timeend.and.date.ge.timestart)
    timeind2=ind(date2.le.timeend.and.date2.ge.timestart)        


    u = short2flt(f1->u(timeind,{1000:100},::-1,:))
    v = short2flt(f2->v(timeind,{1000:100},::-1,:))
;;;;


;;;;calculate velocity potential
    print("计算速度势函数")

    vp_total = uv2sfvpF(u,v)
    sf = vp_total(0,:,:,:,:)
    vp = vp_total(1,:,:,:,:)
    
    copy_VarCoords(u, sf)
    copy_VarCoords(u, vp)
  
  printVarSummary(sf)
  printVarSummary(vp)


;;;写nc文件********************************************************************
    print("****写入nc文件********")
    dim=dimsizes(sf)
    ntime=dim(0)
    nlev=dim(1)
    nlat=dim(2)
    nlon=dim(3)

;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro="/home/ys17-19/lsh/data/era/"
    filo1="era_vp_sf.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "era计算的流函数,速度势,辐散风"    
    fAtt@creation_date = "2020.09.05"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(sf&time),(/"time"/))
    filevardef(fout, "level",typeof(sf&level), (/"level"/))
    filevardef(fout, "lat", typeof(sf&latitude),(/"lat"/))
    filevardef(fout, "lon", typeof(sf&longitude),(/"lon"/))


    filevardef(fout, "sf", typeof(sf),(/"time","level","lat","lon"/))
    filevardef(fout, "vp", typeof(vp),(/"time","level","lat","lon"/))  
 
  
    filevarattdef(fout, "time", sf&time)   ;;;记得给时间附上时间单位
    filevarattdef(fout, "level", sf&level)
    filevarattdef(fout, "lat", sf&latitude)
    filevarattdef(fout, "lon", sf&longitude)

    filevarattdef(fout, "sf", sf)
    filevarattdef(fout, "vp", vp)    
 
   
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/sf&time/)
    fout->level=(/sf&level/)
    fout->lat=(/sf&latitude/)
    fout->lon=(/sf&longitude/)

    fout->sf=(/sf/)
    fout->vp=(/vp/)


end 


