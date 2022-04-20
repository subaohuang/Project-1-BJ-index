;;;************************************************************
;;
;; linshuheng 2021.03.04
;; 计算质量流函数
;;;*********************************************************
begin

model = "20CR"
;;;region and time period 
   timestart  = 187101
   timeend  =201212
;;;;
 
  a=6371000
  g=9.8
  pi=get_pi("double")
  dfi=(pi*10)/180

;;;读入数据
  diri = "/home/ys17-23/Extension/20thCRV2/"

  fils1 = systemfunc("ls "+diri+"uwnd.mon.mean.nc")
  fils2 = systemfunc("ls "+diri+"vwnd.mon.mean.nc")
 

  f1 = addfile(fils1, "r")
  f2 = addfile(fils2, "r")


  time = f1->time
  lat = f1->lat
  lon = f1->lon
  lev = f1->level

  ; levwant= (/1000.0, 925.0, 850.0, 700.0, 600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)   ;level反着取
  levwant= (/950.0,850.0,200.0/)   ;level反着取
  levind1 = get1Dindex(lev,levwant)
  level  = f1->level(levind1)

  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)
  time_r = time(timeind)

   ;;; -------------Calcaulate divergence wind -------------------------------

      uwnd = f1->uwnd(timeind,levind1,::-1,:)  ;;;variable must have ascending lat and lon order coordinates
      vwnd = f2->vwnd(timeind,levind1,::-1,:)   ;;;variable must have ascending lat and lon order coordinates
      printVarSummary(uwnd)
      ;;;
      ; print("开始计算辐散风")
      ; vrdvx = uv2vrdvF(uwnd,vwnd)
      ; uvd = dv2uvF(vrdvx(1,:,:,:,:))
      ; delete(vrdvx)
      ; vd = uvd(1,:,:,:,:)
      ; print("结束计算辐散风")

      print("开始计算辐散风")
      dv =   uv2dvF(uwnd,vwnd)
      uvd = dv2uvF(dv)
      delete(dv)
      vd = uvd(1,:,:,:,:)
      ud = uvd(0,:,:,:,:)

      copy_VarCoords(vwnd, vd)
      copy_VarAtts(vwnd, vd)

      copy_VarCoords(uwnd, ud)
      copy_VarAtts(uwnd, ud)
      delete(uvd)

      print("结束计算辐散风")



      print("计算势函数")

      sfvp = uv2sfvpF(uwnd,vwnd)
      vp = sfvp(1,:,:,:,:)
      
      copy_VarCoords(uwnd, vp)
      copy_VarAtts(uwnd, vp)
      delete(sfvp)

      print("结束计算势函数")
     
     printVarSummary(ud)
     printVarSummary(vd)
     printVarSummary(vp)
      ; vd!0="time"
      ; vd!1="level"
      ; vd!2="lat"
      ; vd!3="lon"
      ; printVarSummary(vd)
      ; delete(uwnd)
      ; delete(vwnd)
       
      dim = dimsizes(vd)
      ntime = dim(0)
      nlev = dim(1)
      nlat = dim(2)
      nlon = dim(3)
      
   print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro="/home/ys17-23/Extension/linshh/obs-derived/"
    filo1="vp_ud_vd_"+model+"_187101-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "Potensial Velocity& divergence wind 850 hpa and 200 hpa"    
    fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "time",typeof(vd&time), (/"time"/))
    filevardef(fout, "level",typeof(vd&level), (/"level"/))
    filevardef(fout, "lat", typeof(vd&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(vd&lon),(/"lon"/))

   
    filevardef(fout, "vd", typeof(vd),(/"time","level","lat","lon"/))
    filevardef(fout, "ud", typeof(ud),(/"time","level","lat","lon"/))
    filevardef(fout, "vp", typeof(vp),(/"time","level","lat","lon"/))

    filevarattdef(fout, "level", vd&level)
    filevarattdef(fout, "lat", vd&lat)
    filevarattdef(fout, "lon", vd&lon)
    filevarattdef(fout, "time", time_r)
    filevarattdef(fout, "vd", vd)
    filevarattdef(fout, "vp", vp)
    filevarattdef(fout, "ud", ud)

    setfileoption(fout,"DefineMode",False)


    fout->time=(/time_r/)
    fout->level=(/vd&level/)
    fout->lat=(/vd&lat/)
    fout->lon=(/vd&lon/)
    fout->vd = (/vd/)
    fout->ud = (/ud/)
    fout->vp = (/vp/)
   delete(vd)
   delete(time_r)
   delete(time)
 print("********************************finish calculate"+" "+model+"************************************")

print("creat successfully")

end 




 