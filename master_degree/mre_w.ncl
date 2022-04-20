;;*************************************************************
;;;;MRE_DATA  combine: w
;**************************************************************
begin
  diri="/mnt/f/share/"
;**************************************************************
;read data
;*************************************************************
print("ncep1")
;;NCEP1 w
    fils1="/mnt/f/share/ncep/ncep1/omega.mon.mean.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    level=f1->level
    lat=f1->lat
    lon=f1->lon
    time=f1->time(timeind)
    w_necp1=f1->omega(timeind,:,::-1,:)   ;;;插值前先逆着取
    
    print("ncep1的垂直速度取值")
    printMinMax(w_necp1, 1)
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
;;;

print("20cr")
;;;20CR
    fils1="/mnt/f/share/ncep/20CR/omega.mon.mean.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    w_r=f1->omega(timeind,:,::-1,:)

    delete(f1)
    delete(fils1)
    ;;;interpolation
    xo=lon
    yo=lat(::-1)

    print("开始水平插值")  
    w_20cr_shuiping=linint2_Wrap(w_r&lon, w_r&lat, w_r, True, xo, yo, 0)
    delete(w_r)
    delete(xo)
    delete(yo)
    w_20cr_shuiping!2="lat"
    w_20cr_shuiping!3="lon"  
   
    ;;;
   
    print("开始垂直插值")
    w_20cr=int2p_n_Wrap(w_20cr_shuiping&level, w_20cr_shuiping, level, 1, 1)
    delete(w_20cr_shuiping)
    w_20cr!1="level"
    w_20cr!1="level"  
 
    printVarSummary(w_20cr)
    delete(date)
    delete(timeind)
    

    print("*********************20cr的垂直速度取值**************************")
    printMinMax(w_20cr, 1)

;;;
    print("jra55")
    ;;;JRA-55
    fils1="/mnt/f/share/JRA-55/omega.jra55.mon.mean.nc"
   
    f1=addfile(fils1, "r")
    levind=get1Dindex(f1->level,level)
    w_r=f1->omega(0:659,levind,::-1,:)  ;;;195801-201212
   
    delete(f1)
   

    xo=fspan(0,357.5,144)
    yo=fspan(-90, 90, 73)

    print("开始水平插值")  
    w_jra55=linint2_Wrap(w_r&lon, w_r&lat, w_r, True, xo, yo, 0)

    w_jra55!2="lat"
    w_jra55!3="lon"

    delete(w_r)

   
   printVarSummary(w_jra55)
   delete(fils1)
   delete(levind)

   
    print("**************jra55的垂直速度取值**********")
    printMinMax(w_jra55, 1)

;;;

    print("****************era40*****************")
    ;;; era40&interium

    fils1="/mnt/f/share/era/era40/era_40_w.nc"
    f1=addfile(fils1, "r")
    levind=get1Dindex(f1->level,level)
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.200112)

    w_40_r=f1->w(timeind,levind,::-1,:)
    w_40=dble2flt(w_40_r*w_40_r@scale_factor+w_40_r@add_offset)
    copy_VarCoords(w_40_r,w_40)
    delete(w_40_r)

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
    w_in_r=f2->w(timeind,levind,::-1,:)
   
    w_in=dble2flt(w_in_r*w_in_r@scale_factor+w_in_r@add_offset)
    copy_VarCoords(w_in_r, w_in)

    w_era=new((/dimsizes(w_40&time)+dimsizes(w_in&time), \
        dimsizes(w_40&level),dimsizes(w_40&latitude),dimsizes(w_40&longitude)/),typeof(w_40))
    w_era!0="time"
    w_era!1="level"
    w_era!2="lat"
    w_era!3="lon"
    

    w_era(0:dimsizes(w_40&time)-1,:,:,:)=(/w_40/)
    w_era(dimsizes(w_40&time):659,:,:,:)=(/w_in/)


    delete(w_40)
    delete(w_in)


    delete(fils2)
    delete(levind)
    delete(date)
    delete(timeind)
    print("****************era垂直速度的取值**************")
    printMinMax(w_era, 1)
;;;

;***********************************************************************
;data_mean
;************************************************************************
    w=w_necp1
    w=(/(w_necp1+w_20cr+w_era+w_jra55)/4.0/)
 
   delete(w_necp1)
   delete(w_20cr)
   delete(w_era)
   delete(w_jra55)
     print("**************平均垂直速度的取值**************")
    printMinMax(w, 1)

  print("write nc")
  ntime=dimsizes(time)
  nlev=dimsizes(level)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
;***********************************************************************
;write nc
;***********************************************************************

    diro=diri
    filo1="mre_w.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "1958-2012,era40,ncep1,jra-55以及20cr的w风场平均结果，精度为2.5乘2.5"    
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
    filevardef(fout, "w", typeof(w),(/"time","level","lat","lon"/))

    filevarattdef(fout, "level", level)
    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "lat", lat)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "w", w)
 
    fout->level=(/level/)
    fout->lon=(/lon/)
    fout->lat=(/lat/)
    fout->time=(/time/)
    fout->w=(/w(:,:,::-1,:)/)
 


end 
   
