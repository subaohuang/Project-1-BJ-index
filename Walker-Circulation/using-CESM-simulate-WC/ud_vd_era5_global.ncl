;;;************************************************************
;;
;; linshuheng 2021.03.04
;; 计算质量流函数
;;;*********************************************************
begin

model = "ERA5"
;;;region and time period 
   timestart  = 197901
   timeend  =202012
;;;;
 
  a=6371000
  g=9.8
  pi=get_pi("double")
  dfi=(pi*10)/180

;;;读入数据
  diri="/home/yangsong3/data-observation/ERA5-monthly/pressure/"
  fils1 = systemfunc("ls "+diri+"u_mon_1x1_1979_2020.nc")
  fils2 = systemfunc("ls "+diri+"v_mon_1x1_1979_2020.nc")
 
  

  f1 = addfile(fils1, "r")
  f2 = addfile(fils2, "r")

 

  time = f1->time
  lat = f1->latitude
  lon = f1->longitude
  lev = f1->level
  levwant= (/1000.0, 925.0, 850.0, 700.0, 600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)   ;level反着取
  
  levind1 = get1Dindex(lev,levwant)
  
  level  = f1->level(levind1)


  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)
  time_r = time(timeind)



   ;;; -------------Calcaulate divergence wind -------------------------------

      uwnd = short2flt(f1->u(timeind,levind1,::-1,:))   ;;;variable must have ascending lat and lon order coordinates
      vwnd = short2flt(f2->v(timeind,levind1,::-1,:))   ;;;variable must have ascending lat and lon order coordinates
      printVarSummary(uwnd)
      ;;;
      print("开始计算辐散风")
      vrdvx = uv2vrdvF(uwnd,vwnd)
      uvd = dv2uvF(vrdvx(1,:,:,:,:))
      delete(vrdvx)
      vd = uvd(1,:,:,:,:)
      print("结束计算辐散风")

      vd=uvd(1,:,:,:,:)  ;;;time lev,lat,lon
      delete(uvd)
 
      copy_VarCoords(vwnd, vd)
      copy_VarAtts(vwnd, vd)

      vd!0="time"
      vd!1="level"
      vd!2="lat"
      vd!3="lon"
      printVarSummary(vd)
      delete(uwnd)
      delete(vwnd)
       
      dim = dimsizes(vd)
      ntime = dim(0)
      nlev = dim(1)
      nlat = dim(2)
      nlon = dim(3)
       


   print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro="/home/yangsong3/data-observation/linshh/data/rea_msf/"
    filo1="vd_"+model+"_197901-202012.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "CMIP6 historical run Meridional Divergence wind"    
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
    filevarattdef(fout, "level", vd&level)
    filevarattdef(fout, "lat", vd&lat)
    filevarattdef(fout, "lon", vd&lon)
    filevarattdef(fout, "time", time_r)
    filevarattdef(fout, "vd", vd)
 
    
    setfileoption(fout,"DefineMode",False)


    fout->time=(/time_r/)
    fout->level=(/vd&level/)
    fout->lat=(/vd&lat/)
    fout->lon=(/vd&lon/)
    fout->vd = (/vd/)
    delete(vd)
   delete(time_r)
   delete(time)
 print("********************************finish calculate"+" "+model+"************************************")

print("vd creat successfully")

end 




 