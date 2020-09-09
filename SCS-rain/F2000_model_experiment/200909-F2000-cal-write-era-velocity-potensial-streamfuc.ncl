
;;********************************************************
;;计算速度势函数以及辐散风场
;; ERA-interim:u,v，
;; used ncl function: uv2dv_cfd,
;；；grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2020.09.04
;**********************************************************

undef("cal_vp_sf")
function cal_vp_sf(infilepath,outfilepath,timestart,timeend)
begin
    fils1 = infilepath
    f1 = addfile(fils1, "r")

;;;;read data
  
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.le.timeend.and.date.ge.timestart)
      
    u = f1->U(timeind,{1000:100},:,:)
    v = f1->V(timeind,{1000:100},:,:)

   
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
    system("/bin/rm -f "+outfilepath+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(outfilepath,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "era计算的流函数,速度势,辐散风"    
    fAtt@creation_date = "2020.09.05"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lev","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(sf&time),(/"time"/))
    filevardef(fout, "lev",typeof(sf&lev), (/"lev"/))
    filevardef(fout, "lat", typeof(sf&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(sf&lon),(/"lon"/))


    filevardef(fout, "sf", typeof(sf),(/"time","lev","lat","lon"/))
    filevardef(fout, "vp", typeof(vp),(/"time","lev","lat","lon"/))  
 
  
    filevarattdef(fout, "time", sf&time)   ;;;记得给时间附上时间单位
    filevarattdef(fout, "lev", sf&lev)
    filevarattdef(fout, "lat", sf&lat)
    filevarattdef(fout, "lon", sf&lon)

    filevarattdef(fout, "sf", sf)
    filevarattdef(fout, "vp", vp)    
 
   
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/sf&time/)
    fout->lev=(/sf&lev/)
    fout->lat=(/sf&lat/)
    fout->lon=(/sf&lon/)

    fout->sf=(/sf/)
    fout->vp=(/vp/)
;;;
  
    k= "处理完成"
    print("处理完成")
    return(k)
end 
;;***************************************************************************
;;main program 
;;*************************************************************************

begin

  m = cal_vp_sf(infilepath,outfilepath,101,4012)
  
end 


