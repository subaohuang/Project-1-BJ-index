
;;********************************************************
;;计算速度势函数以及辐散风场
;; ERA-interim:u,v，
;; used ncl function: uv2dv_cfd,
;；；grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2020.09.04
;**********************************************************
;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
begin
  ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"uwnd","vwnd"/)
    filein = path_wind + var_name +".mon.mean.nc"
    
    outfilepath = "/home/ys17-19/lsh/Project/SCS-rain/annual/era_interium_sf_vp.nc"
 
    ;;;region
    lat1 = -5
    lat2 = 45
    lon1 = 90
    lon2 = 170


    ;;;time span
    startmon =197901
    endmon = 201612

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1



    fils1 = filein(0)
    fils2 = filein(1)
    f1 = addfile(fils1, "r")
    f2 = addfile(fils2, "r")

;;;;read data
  
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)
      
    u = short2flt(f1->u(timeind,{1000:100},:,:))
    v = short2flt(f2->v(timeind,{1000:100},:,:))

   
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
;;;;

;;;计算旋度风



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
    filevardef(fout, "lev",typeof(sf&level), (/"lev"/))
    filevardef(fout, "lat", typeof(sf&latitude),(/"lat"/))
    filevardef(fout, "lon", typeof(sf&longitude),(/"lon"/))


    filevardef(fout, "sf", typeof(sf),(/"time","lev","lat","lon"/))
    filevardef(fout, "vp", typeof(vp),(/"time","lev","lat","lon"/))  
 
  
    filevarattdef(fout, "time", sf&time)   ;;;记得给时间附上时间单位
    filevarattdef(fout, "lev", sf&level)
    filevarattdef(fout, "lat", sf&latitude)
    filevarattdef(fout, "lon", sf&longitude)

    filevarattdef(fout, "sf", sf)
    filevarattdef(fout, "vp", vp)    
 
   
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/sf&time/)
    fout->lev=(/sf&level/)
    fout->lat=(/sf&latitude/)
    fout->lon=(/sf&longitude/)

    fout->sf=(/sf/)
    fout->vp=(/vp/)
;;;
  
    k= "处理完成"
    print("处理完成")

end 
;;***************************************************************************
;;main program 
;;*************************************************************************



