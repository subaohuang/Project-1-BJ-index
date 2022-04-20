;;*************************************************************
;;;;MRE_DATA  combine: U,V
;**************************************************************
begin
  diri="/mnt/f/share/"
;**************************************************************
;read data
;*************************************************************
print("ncep1")
;;NCEP1 U,V 
    fils1="/mnt/f/share/ncep/ncep1/uwnd.mon.mean.nc"
    fils2="/mnt/f/share/ncep/ncep1/vwnd.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    level=f1->level
    lat=f1->lat
    lon=f1->lon
    time=f1->time(timeind)
    uwnd_necp1=f1->uwnd(timeind,:,::-1,:)   ;;;插值前先逆着取
    vwnd_necp1=f2->vwnd(timeind,:,::-1,:)   ;;;
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
    delete(fils2)
;;;

print("20cr")
;;;20CR
    fils1="/mnt/f/share/ncep/20CR/uwnd.mon.mean.nc"
    fils2="/mnt/f/share/ncep/20CR/vwnd.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)

    uwnd_r=f1->uwnd(timeind,:,::-1,:)
    vwnd_r=f2->vwnd(timeind,:,::-1,:)
    delete(f1)
    delete(f2)
    delete(fils1)
    delete(fils2)
    ;;;interpolation
    xo=lon
    yo=lat(::-1)

    print("开始水平插值")  
    uwnd_20cr_shuiping=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, True, xo, yo, 0)
    vwnd_20cr_shuiping=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, True, xo, yo, 0)
    delete(uwnd_r)
    delete(vwnd_r)
    delete(xo)
    delete(yo)
    uwnd_20cr_shuiping!2="lat"
    uwnd_20cr_shuiping!3="lon"  
    vwnd_20cr_shuiping!2="lat"
    vwnd_20cr_shuiping!3="lon"
    printVarSummary(uwnd_20cr_shuiping)
    ;;;
   
    print("开始垂直插值")
    
   
    uwnd_20cr=int2p_n_Wrap(uwnd_20cr_shuiping&level, uwnd_20cr_shuiping, level, 1, 1)
    vwnd_20cr=int2p_n_Wrap(vwnd_20cr_shuiping&level, vwnd_20cr_shuiping, level, 1, 1)

    delete(uwnd_20cr_shuiping)
    delete(vwnd_20cr_shuiping)

    uwnd_20cr!1="level"
    uwnd_20cr!1="level"  
    vwnd_20cr!1="level"
    vwnd_20cr!1="level"
    
    printVarSummary(uwnd_20cr)
    delete(date)
    delete(timeind)


;;;
    print("jra55")
    ;;;JRA-55
    fils1="/mnt/f/share/JRA-55/uwnd.jra55.mon.mean.nc"
    fils2="/mnt/f/share/JRA-55/vwnd.jra55.mon.mean.nc"
        
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    levind=get1Dindex(f1->level,level)
    uwnd_r=f1->uwnd(0:659,levind,::-1,:)  ;;;195801-201212
    vwnd_r=f2->vwnd(0:659,levind,::-1,:)  ;;;
    delete(f1)
    delete(f2)

    xo=fspan(0,357.5,144)
    yo=fspan(-90, 90, 73)

   print("开始水平插值")  
   uwnd_jra55=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, True, xo, yo, 0)
   vwnd_jra55=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, True, xo, yo, 0)
   uwnd_jra55!2="lat"
   uwnd_jra55!3="lon"
   vwnd_jra55!2="lat"
   vwnd_jra55!3="lon"

   delete(uwnd_r)
   delete(vwnd_r)
   
   printVarSummary(uwnd_jra55)
   printVarSummary(vwnd_jra55)
   delete(fils1)
   delete(fils2)
   delete(levind)
;;;

    print("era40")
    ;;; era40&interium

    fils1="/mnt/f/share/era/era40/era_40_u_v.nc"
    f1=addfile(fils1, "r")
    levind=get1Dindex(f1->level,level)
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.200112)

    uwnd_40_r=f1->u(timeind,levind,::-1,:)
    vwnd_40_r=f1->v(timeind,levind,::-1,:)
     
    uwnd_40=dble2flt(uwnd_40_r*uwnd_40_r@scale_factor+uwnd_40_r@add_offset)
    vwnd_40=dble2flt(vwnd_40_r*vwnd_40_r@scale_factor+vwnd_40_r@add_offset) 
    copy_VarCoords(uwnd_40_r, uwnd_40)
    copy_VarCoords(vwnd_40_r, vwnd_40)
    delete(uwnd_40_r)
    delete(vwnd_40_r)
   
    delete(fils1)
    delete(f1)
    delete(levind)
    delete(timeind)
    delete(date)
    ;;;读入era_in的数据
    fils2= "/mnt/f/share/era/data/data_series/era_pre_global.nc"
    f2=addfile(fils2, "r")
    levind=get1Dindex(f2->level,level)    ;;;这里容易出现错误，请注意
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.200201.and.date.le.201212) 
    uwnd_in_r=f2->u(timeind,levind,::-1,:)
    vwnd_in_r=f2->v(timeind,levind,::-1,:)
    
    uwnd_in=dble2flt(uwnd_in_r*uwnd_in_r@scale_factor+uwnd_in_r@add_offset)
    vwnd_in=dble2flt(vwnd_in_r*vwnd_in_r@scale_factor+vwnd_in_r@add_offset)
    copy_VarCoords(uwnd_in_r, uwnd_in)
    copy_VarCoords(vwnd_in_r, vwnd_in)
   
    uwnd_era=new((/dimsizes(uwnd_40&time)+dimsizes(uwnd_in&time), \
        dimsizes(uwnd_40&level),dimsizes(uwnd_40&latitude),dimsizes(uwnd_40&longitude)/),typeof(uwnd_40))
    uwnd_era!0="time"
    uwnd_era!1="level"
    uwnd_era!2="lat"
    uwnd_era!3="lon"
    
    vwnd_era=uwnd_era

    uwnd_era(0:dimsizes(uwnd_40&time)-1,:,:,:)=(/uwnd_40/)
    uwnd_era(dimsizes(uwnd_40&time):659,:,:,:)=(/uwnd_in/)

    vwnd_era(0:dimsizes(vwnd_40&time)-1,:,:,:)=(/vwnd_40/)
    vwnd_era(dimsizes(vwnd_40&time):659,:,:,:)=(/vwnd_in/)
    
    delete(uwnd_40)
    delete(uwnd_in)
    delete(vwnd_40)
    delete(vwnd_in)

    delete(fils2)
    delete(levind)
    delete(date)
    delete(timeind)
;;;

;;;
    print("ncep2")
    ;;NCEP1 U,V 
    fils1="/mnt/f/share/ncep/ncep1/uwnd.mon.mean.nc"
    fils2="/mnt/f/share/ncep/ncep1/vwnd.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    level=f1->level
    lat=f1->lat
    lon=f1->lon
    time=f1->time(timeind)
    uwnd_necp1=f1->uwnd(timeind,:,::-1,:)   ;;;插值前先逆着取
    vwnd_necp1=f2->vwnd(timeind,:,::-1,:)   ;;;
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
    delete(fils2)
;;;



;***********************************************************************
;data_mean
;************************************************************************
    uwnd=uwnd_necp1
    vwnd=vwnd_necp1

    uwnd=(/(uwnd_necp1+uwnd_20cr+uwnd_era+uwnd_jra55)/4.0/)
    vwnd=(/(vwnd_necp1+vwnd_20cr+vwnd_era+vwnd_jra55)/4.0/)
   delete(uwnd_necp1)
   delete(uwnd_20cr)
   delete(uwnd_era)
   delete(uwnd_jra55)

   delete(vwnd_necp1)
   delete(vwnd_20cr)
   delete(vwnd_era)
   delete(vwnd_jra55)

  print("write nc")
  ntime=dimsizes(time)
  nlev=dimsizes(level)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
;***********************************************************************
;write nc
;***********************************************************************

    diro=diri
    filo1="mre_u_v.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "1958-2012,era40,ncep1,jra-55以及20cr的风场平均结果，精度为2.5乘2.5"    
    fAtt@creation_date = "2017.12.21"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    dimsize=(/ntime,nlev,nlat,nlon/)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    filevardef(fout, "level",typeof(level), (/"level"/))
    filevardef(fout, "lat", typeof(lat),(/"lat"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(time),(/"time"/))
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "uwnd", typeof(uwnd),(/"time","level","lat","lon"/))
    filevardef(fout, "vwnd", typeof(vwnd),(/"time","level","lat","lon"/))
   
    filevarattdef(fout, "level", level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "lat", lat)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "uwnd", uwnd)
    filevarattdef(fout, "vwnd", vwnd)

    fout->level=(/level/)
    fout->lon=(/lon/)
    fout->lat=(/lat/)
    fout->time=(/time/)
    fout->uwnd=(/uwnd(:,:,::-1,:)/)
    fout->vwnd=(/vwnd(:,:,::-1,:)/)


end 
   
