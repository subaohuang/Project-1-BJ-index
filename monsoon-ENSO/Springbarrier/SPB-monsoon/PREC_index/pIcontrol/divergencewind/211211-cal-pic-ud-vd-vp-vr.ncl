;;;************************************************************
;;
;; linshuheng 2021.03.04
;; 计算质量流函数
;;;*********************************************************
begin

model = "CESM2"
;;;region and time period 
    timestart= 101
    timeend = 120012
  
;;;;
 
  a=6371000
  g=9.8
  pi=get_pi("double")
  dfi=(pi*10)/180

;;;读入数据

  diri="/home/ys17-23/CMIP6/2.5x2.5/CESM2/piControl/"

  fils1 = systemfunc("ls "+diri+"ua_Amon_CESM2_piControl_r1i1p1f1_gn_"+"*")
  fils2 = systemfunc("ls "+diri+"va_Amon_CESM2_piControl_r1i1p1f1_gn_"+"*")
 

  f1 = addfiles(fils1, "r")
  f2 = addfiles(fils2, "r")


  time = f1[:]->time
  lat = f1[0]->lat
  lon = lonFlip(f1[0]->lon)
  lev = f1[0]->plev

  ; levwant= (/1000.0, 925.0, 850.0, 700.0, 600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)   ;plev反着取
  levwant= (/92500.0,85000.0,20000.0/)   ;plev反着取
  levind1 = get1Dindex(lev,levwant)


  plev2  = f1[0]->plev(levind1)

  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)
  time_r = time(timeind)
  
;;;;

;;;;
   ;;; -------------Calcaulate divergence wind -------------------------------

      uwnd = lonFlip(f1[:]->ua(timeind,levind1,:,:))  ;;;variable must have ascending lat and lon order coordinates
      vwnd = lonFlip(f2[:]->va(timeind,levind1,:,:))   ;;;variable must have ascending lat and lon order coordinates
      printVarSummary(uwnd)
      ;;;


    ;;;interp missing value
     print("****************interp missing value*********************")
       
    ;******************************************************************
    ; use poisson_grid_fill to fill missing value
    ;******************************************************************
    nscan = 500           ; usually *much* fewer
    eps   = 0.001            ; variable depended
    gtype = False            ; "gendat" does not generate cyclic fields
    guess = 0                ; use zonal means
    relc  = 0.6              ; standard relaxation coef
    opt   = 0
    
    poisson_grid_fill(uwnd, gtype, guess, nscan, eps, relc, opt)
    
    poisson_grid_fill(vwnd, gtype, guess, nscan, eps, relc, opt)
  ;;;;   
  
      vrdvx = uv2vrdvF(uwnd,vwnd) 
    
      vor =   vrdvx(0,:,:,:,:)  ;;; vorticity  
      dv =   vrdvx(1,:,:,:,:)  ;;; vorticity  

      copy_VarCoords(uwnd, vor)
      copy_VarCoords(uwnd, dv)

      print("开始计算辐散风")
      uvd = dv2uvF_Wrap(dv)
      delete(dv)
      vd = uvd(1,:,:,:,:)
      ud = uvd(0,:,:,:,:)

      delete(uvd)

      print("结束计算辐散风")

      print("开始计算旋度风")
      uvr = vr2uvF_Wrap(vor)
      delete(vor)
      vr = uvr(1,:,:,:,:)
      ur = uvr(0,:,:,:,:)

      delete(uvr)

      print("结束计算旋转风")


      print("计算势函数&流函数")

      sfvp = uv2sfvpF(uwnd,vwnd)

      sf = sfvp(0,:,:,:,:)
      vp = sfvp(1,:,:,:,:)

      copy_VarCoords(uwnd, vp)
      copy_VarAtts(uwnd, vp)

      copy_VarCoords(uwnd, sf)
      copy_VarAtts(uwnd, sf)
      delete(sfvp)
      print("结束计算势函数")


     
     printVarSummary(ud)
     printVarSummary(vd)

     printVarSummary(ur)
     printVarSummary(vr)

     printVarSummary(vp)
     printVarSummary(sf)
   
       
    dim = dimsizes(vd)
    ntime = dim(0)
    nlev = dim(1)
    nlat = dim(2)
    nlon = dim(3)
      
   print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro="/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/divergencewind/"
    filo1="vpsf_uvd_uvr_"+model+"_1-1200.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "VPSF& Rotated_divergence wind 850 hpa and 200 hpa"    
    fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","plev","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time",typeof(vd&time), (/"time"/))
    filevardef(fout, "plev",typeof(vd&plev), (/"plev"/))
    filevardef(fout, "lat", typeof(vd&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(vd&lon),(/"lon"/))

   
    filevardef(fout, "vd", typeof(vd),(/"time","plev","lat","lon"/))
    filevardef(fout, "ud", typeof(ud),(/"time","plev","lat","lon"/))

    filevardef(fout, "vr", typeof(vr),(/"time","plev","lat","lon"/))
    filevardef(fout, "ur", typeof(ur),(/"time","plev","lat","lon"/))

    filevardef(fout, "vp", typeof(vp),(/"time","plev","lat","lon"/))
    filevardef(fout, "sf", typeof(sf),(/"time","plev","lat","lon"/))

    filevarattdef(fout, "plev", vd&plev)
    filevarattdef(fout, "lat", vd&lat)
    filevarattdef(fout, "lon", vd&lon)
    filevarattdef(fout, "time", time_r)

    filevarattdef(fout, "vd", vd)
    filevarattdef(fout, "vp", vp)
    filevarattdef(fout, "ud", ud)

    filevarattdef(fout, "vr", vr)
    filevarattdef(fout, "sf", sf)
    filevarattdef(fout, "ur", ur)

    setfileoption(fout,"DefineMode",False)


    fout->time=(/time_r/)
    fout->plev=(/vd&plev/)
    fout->lat=(/vd&lat/)
    fout->lon=(/vd&lon/)

    fout->vd = (/vd/)
    fout->ud = (/ud/)
    fout->vp = (/vp/)

    fout->vr = (/vr/)
    fout->ur = (/ur/)
    fout->sf = (/sf/)

   delete(time_r)
   delete(time)

 print("********************************finish calculate"+" "+model+"************************************")

print("creat successfully")

end 




 