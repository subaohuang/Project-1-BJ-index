;;*************************************************************
;;;;MRE_DATA  combine: w
;**************************************************************
begin
  diri="/mnt/f/share/"
;**************************************************************
;read data
;*************************************************************
xo=fspan(0,357.5,144)
yo=fspan(-87.5, 87.5, 71)
print("ncep1")
;;NCEP1 w
    fils1="/mnt/f/share/ncep/ncep1/omega.mon.mean.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    level=f1->level
    time=f1->time(timeind)
    w_r=f1->omega(timeind,:,::-1,:)   ;;;插值前先逆着取
    
    print("开始水平插值")  
    w_necp1=linint2_Wrap(w_r&lon, w_r&lat, w_r, True, xo, yo, 0)
    copy_VarCoords(w_r, w_necp1)
    delete(w_r)
    w_necp1!2="lat"
    w_necp1!3="lon"
    printVarSummary(w_necp1)

    lat=w_necp1&lat(::-1)
    lon=w_necp1&lon

    lat@units="degree_north"
    lon@units="degree_east"

    print("ncep1的垂直速度取值")
    printMinMax(w_necp1, 1)
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
;;;
print("*********************ncep2*****************************")
    ;;NCEP1 U,V 
    fils1="/mnt/f/share/ncep/ncep2/omega.mon.mean.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    time=f1->time(timeind)
    levind=get1Dindex(f1->level,level)
    w_r=f1->omega(timeind,levind,::-1,:)   ;;;插值前先逆着取
    w_in=tofloat(w_r*w_r@scale_factor+w_r@add_offset)
    copy_VarCoords(w_r, w_in)
    print("开始水平插值")  
    w_necp2=linint2_Wrap(w_in&lon, w_in&lat, w_in, True, xo, yo, 0)
    delete(w_in)

    w_necp2!2="lat"
    w_necp2!3="lon"

    delete(w_r)
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
    timeind=ind(date.ge.197901.and.date.le.201212)
    w_r=f1->omega(timeind,:,::-1,:)

    delete(f1)
    delete(fils1)
    ;;;interpolation

    print("开始水平插值")  
    w_20cr_shuiping=linint2_Wrap(w_r&lon, w_r&lat, w_r, True, xo, yo, 0)
    delete(w_r)
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

print("**************jra55**********")
    print("jra55")
    ;;;JRA-55
    fils1="/mnt/f/share/JRA-55/omega.jra55.mon.mean.nc"
   
    f1=addfile(fils1, "r")
    levind=get1Dindex(f1->level,level)
    w_r=f1->omega(252:659,levind,::-1,:)  ;;;195801-201212
   
    delete(f1)
   

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

    print("****************era*****************")
    ;;; era40&interium

    ;;;读入era_in的数据
    fils2= "/mnt/f/share/era/data/data_series/era_pre_global.nc"
    f2=addfile(fils2, "r")
    levind=get1Dindex(f2->level,level)    ;;;这里容易出现错误，请注意
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212) 
    w_in_r=f2->w(timeind,levind,::-1,:)
   
    w_in=dble2flt(w_in_r*w_in_r@scale_factor+w_in_r@add_offset)
    copy_VarCoords(w_in_r, w_in)
     delete(w_in_r)

    print("开始水平插值")  
    w_era=linint2_Wrap(w_in&longitude, w_in&latitude, w_in, True, xo, yo, 0)
     delete(w_in)

    w_era!0="time"
    w_era!1="level"
    w_era!2="lat"
    w_era!3="lon"
    
    printVarSummary(w_era)

    delete(fils2)
    delete(levind)
    delete(date)
    delete(timeind)
    print("****************era垂直速度的取值**************")
    printMinMax(w_era, 1)
;;;


;;; ;;;读入MERRA的数据
   print("**************************MERRA****************************")
   diri="/mnt/f/share/MERRA/MERRA/"
  m=""+ispan(1979,2012,1)+""

  fils=new((/34,12/),string)
  do i=0,33

  if i.le.13 then
  fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA/MERRA100.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
  end if
  if i.ge.14.and.i.le.21 then 
  fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA/MERRA200.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
  end if
  if i.ge.22 then 
  fils3=systemfunc("ls /mnt/f/share/MERRA/MERRA/MERRA300.prod.assim.instM_3d_asm_Cp."+m(i)+"*"+".SUB.nc4")
  end if
  fils(i,:)= fils3
  end do
  
  fils_new=onedtond(fils, 34*12)
  fileh=new(408,string)
  fileh(0:407)=fils_new
  delete(fils_new)

  f1=addfiles(fileh, "join")
  delete(fileh)
  printVarSummary(f1)
  date=cd_calendar(f1[:]->TIME, 1)
  timeind=ind(date.le.201212.and.date.ge.197901)
  levind=get1Dindex(f1[:]->Height,level)

    w_r=lonPivot(f1[:]->OMEGA(timeind,levind,::-1,:),0)
    printVarSummary(w_r)

    print("开始水平插值")  
    w_merra=linint2_Wrap(w_r&lon, w_r&lat, w_r, True, xo, yo, 0)
 
    w_merra!0="time"
    w_merra!1="level"
    w_merra!2="lat"
    w_merra!3="lon"
    printVarSummary(w_merra)
    delete(w_r)

    delete(f1)
    delete(timeind)
    delete(date)




;***********************************************************************
;data_mean
;************************************************************************
    w=w_necp1
    printVarSummary(w_necp1)
    printVarSummary(w_20cr)
    printVarSummary(w_era)
    printVarSummary(w_jra55)
    printVarSummary(w_necp2)
    
    w=(/(w_necp1+w_20cr+w_era+w_jra55+w_merra+w_necp2)/6.0/)
 
   delete(w_necp1)
   delete(w_20cr)
   delete(w_era)
   delete(w_jra55)
   delete(w_necp2)
   delete(w_merra)
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
    filo1="mre_w_all.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "1958-2012,era40,ncep1,jra-55以及20cr,merra,necp2的w风场平均结果，精度为2.5乘2.5"    
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
   
