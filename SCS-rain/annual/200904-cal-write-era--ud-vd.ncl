
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
   
    print("计算辐散风")

    dv=uv2dvF(u,v)  ;;;divergence
    uvd=dv2uvF(dv)

    ud=uvd(0,:,:,:,:)
    vd=uvd(1,:,:,:,:)
    delete(uvd)
    
    copy_VarCoords(u, ud)
    copy_VarCoords(v, vd)
    

    printVarSummary(ud)
    printVarSummary(vd)


;;;写nc文件********************************************************************
    print("****写入nc文件********")
    dim=dimsizes(ud)
    ntime=dim(0)
    nlev=dim(1)
    nlat=dim(2)
    nlon=dim(3)

;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro="/home/ys17-19/lsh/data/era/"
    filo1="era_ud_vd.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "era计算的辐散风"    
    fAtt@creation_date = "2020.09.05"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(ud&time),(/"time"/))
    filevardef(fout, "level",typeof(ud&level), (/"level"/))
    filevardef(fout, "lat", typeof(ud&latitude),(/"lat"/))
    filevardef(fout, "lon", typeof(ud&longitude),(/"lon"/))

    filevardef(fout, "ud", typeof(ud),(/"time","level","lat","lon"/))
    filevardef(fout, "vd", typeof(vd),(/"time","level","lat","lon"/))
  
 
  
    filevarattdef(fout, "time", ud&time)   ;;;记得给时间附上时间单位
    filevarattdef(fout, "level", ud&level)
    filevarattdef(fout, "lat", ud&latitude)
    filevarattdef(fout, "lon", ud&longitude)

    filevarattdef(fout, "ud", ud)
    filevarattdef(fout, "vd", vd)
   
   
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/ud&time/)
    fout->level=(/ud&level/)
    fout->lat=(/ud&latitude/)
    fout->lon=(/ud&longitude/)


    fout->ud=(/ud/)
    fout->vd=(/vd/)

end 


