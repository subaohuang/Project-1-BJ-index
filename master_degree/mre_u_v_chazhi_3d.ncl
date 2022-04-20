;;*************************************************************
;;;;MRE_DATA  combine: U,V
;**************************************************************
begin
  diri_out="/mnt/f/share/era/data/"
;**************************************************************
;read data
;*************************************************************

xo=fspan(0,357.5,144)
yo=fspan(-87.5, 87.5, 71)

print("************ncep1*********************")
;;NCEP1 U,V 
    fils1="/mnt/f/share/ncep/ncep1/uwnd.mon.mean.nc"
    fils2="/mnt/f/share/ncep/ncep1/vwnd.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    time=f1->time(timeind)
    uwnd_r=f1->uwnd(timeind,{1000},::-1,:)   ;;;插值前先逆着取
    vwnd_r=f2->vwnd(timeind,{1000},::-1,:)   ;;;

    print("开始水平插值")  
    uwnd_necp1=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, True, xo, yo, 0)
    vwnd_necp1=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, True, xo, yo, 0)
    uwnd_necp1!1="lat"
    uwnd_necp1!2="lon"
    vwnd_necp1!1="lat"
    vwnd_necp1!2="lon"

    lat=uwnd_necp1&lat(::-1)
    lon=uwnd_necp1&lon
    lat@units="degree north"
    lon@units="degree east"

    delete(uwnd_r)
    delete(vwnd_r)
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
    delete(fils2)
;;;

print("**************20cr******************************")
;;;20CR
    fils1="/mnt/f/share/ncep/20CR/uwnd.mon.mean.nc"
    fils2="/mnt/f/share/ncep/20CR/vwnd.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)

    uwnd_r=f1->uwnd(timeind,{1000},::-1,:)
    vwnd_r=f2->vwnd(timeind,{1000},::-1,:)
    delete(f1)
    delete(f2)
    delete(fils1)
    delete(fils2)
    ;;;interpolation

    print("开始水平插值")  
    uwnd_20cr=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, True, xo, yo, 0)
    vwnd_20cr=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, True, xo, yo, 0)
    delete(uwnd_r)
    delete(vwnd_r)

    uwnd_20cr!1="lat"
    uwnd_20cr!2="lon"  
    vwnd_20cr!1="lat"
    vwnd_20cr!2="lon"
    printVarSummary(uwnd_20cr)
    ;;;
   
    delete(date)
    delete(timeind)


;;;
    print("*******************jra55***********************************")
    ;;;JRA-55
    fils1="/mnt/f/share/JRA-55/uwnd.jra55.mon.mean.nc"
    fils2="/mnt/f/share/JRA-55/vwnd.jra55.mon.mean.nc"
        
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    levind=get1Dindex(f1->level,1000)
    uwnd_r=f1->uwnd(252:659,levind,::-1,:)  ;;;195801-201212
    vwnd_r=f2->vwnd(252:659,levind,::-1,:)  ;;;
    delete(f1)
    delete(f2)

   print("开始水平插值")  
   uwnd_jra55=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, True, xo, yo, 0)
   vwnd_jra55=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, True, xo, yo, 0)
   uwnd_jra55!1="lat"
   uwnd_jra55!2="lon"
   vwnd_jra55!1="lat"
   vwnd_jra55!2="lon"

   delete(uwnd_r)
   delete(vwnd_r)
   
   printVarSummary(uwnd_jra55)
   printVarSummary(vwnd_jra55)
   delete(fils1)
   delete(fils2)
   delete(levind)
;;;

    ;;;读入era_in的数据
    fils1= "/mnt/f/share/era/data/data_series/era_pre_global.nc"
    f1=addfile(fils1, "r")
    levind=get1Dindex(f1->level,1000)    ;;;这里容易出现错误，请注意
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212) 
    uwnd_r=f1->u(timeind,levind,::-1,:)
    vwnd_r=f1->v(timeind,levind,::-1,:)
    
    uwnd_in=dble2flt(uwnd_r*uwnd_r@scale_factor+uwnd_r@add_offset)
    vwnd_in=dble2flt(vwnd_r*vwnd_r@scale_factor+vwnd_r@add_offset)
    copy_VarCoords(uwnd_r, uwnd_in)
    copy_VarCoords(vwnd_r, vwnd_in)
    delete(uwnd_r)
    delete(vwnd_r)
    dnames = getVarDimNames(uwnd_in)
     
    if any(dnames.eq."latitude")
    uwnd_in!1="lat"
    uwnd_in!2="lon"
    vwnd_in!1="lat"
    vwnd_in!2="lon"
    end if

    print("开始水平插值")  
    uwnd_era=linint2_Wrap(uwnd_in&lon, uwnd_in&lat, uwnd_in, True, xo, yo, 0)
    vwnd_era=linint2_Wrap(vwnd_in&lon, vwnd_in&lat, vwnd_in, True, xo, yo, 0)
    uwnd_era!1="lat"
    uwnd_era!2="lon"
    vwnd_era!1="lat"
    vwnd_era!2="lon"
    delete(uwnd_in)
    delete(vwnd_in)

    delete(fils1)
    delete(levind)
    delete(date)
    delete(timeind)
;;;

;;;
    print("*********************ncep2*****************************")
    ;;NCEP1 U,V 
    fils1="/mnt/f/share/ncep/ncep2/uwnd.mon.mean.nc"
    fils2="/mnt/f/share/ncep/ncep2/vwnd.mon.mean.nc"
    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)

    time=f1->time(timeind)
    uwnd_r=f1->uwnd(timeind,{1000},::-1,:)   ;;;插值前先逆着取
    vwnd_r=f2->vwnd(timeind,{1000},::-1,:)   ;;;

    uwnd_n=uwnd_r*uwnd_r@scale_factor+uwnd_r@add_offset
    vwnd_n=vwnd_r*vwnd_r@scale_factor+vwnd_r@add_offset

    copy_VarCoords(uwnd_r, uwnd_n)
    copy_VarCoords(vwnd_r, vwnd_n)
  
    delete(uwnd_r)
    delete(vwnd_r)
    print("开始水平插值")  
    uwnd_necp2=linint2_Wrap(uwnd_n&lon, uwnd_n&lat, uwnd_n, True, xo, yo, 0)
    vwnd_necp2=linint2_Wrap(vwnd_n&lon, vwnd_n&lat, vwnd_n, True, xo, yo, 0)
    uwnd_necp2!1="lat"
    uwnd_necp2!2="lon"
    vwnd_necp2!1="lat"
    vwnd_necp2!2="lon"
    
    delete(uwnd_n)
    delete(vwnd_n)
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
    delete(fils2)
;;;
  
;;; ;;;读入MERRA的数据
   print("**************************MERRA****************************8")
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
 
 print(date(timeind))
    uwnd_r=lonPivot(f1[:]->U(timeind,{1000},::-1,:),0)

    vwnd_r=lonPivot(f1[:]->V(timeind,{1000},::-1,:),0)



    print("开始水平插值")  
    uwnd_merra=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, True, xo, yo, 0)
    vwnd_merra=linint2_Wrap(vwnd_r&lon, vwnd_r&lat, vwnd_r, True, xo, yo, 0)
    uwnd_merra!1="lat"
    uwnd_merra!2="lon"
    vwnd_merra!1="lat"
    vwnd_merra!2="lon"
    delete(uwnd_r)
    delete(vwnd_r)

    delete(f1)
    delete(timeind)
    delete(date)
   


;***********************************************************************
;data_mean
;************************************************************************
    uwnd=uwnd_necp1
    vwnd=vwnd_necp1

    uwnd=(/(uwnd_necp1+uwnd_20cr+uwnd_era+uwnd_jra55+uwnd_merra+uwnd_necp2)/6.0/)
    vwnd=(/(vwnd_necp1+vwnd_20cr+vwnd_era+vwnd_jra55+uwnd_merra+uwnd_necp2)/6.0/)
   delete(uwnd_necp1)
   delete(uwnd_20cr)
   delete(uwnd_era)
   delete(uwnd_jra55)
   delete(uwnd_merra)
   delete(uwnd_necp2)

   delete(vwnd_necp1)
   delete(vwnd_20cr)
   delete(vwnd_era)
   delete(vwnd_jra55)
   delete(vwnd_merra)
   delete(vwnd_necp2)

  print("write nc")
  ntime=dimsizes(time)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)

;***********************************************************************
;write nc
;***********************************************************************

    diro=diri_out
    filo1="mre_u_v_chazhi_1000.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "1979-2012,era40,ncep1,jra-55,20cr,merra,ncep1的风场平均结果，精度为2.5乘2.5"    
    fAtt@creation_date = "2017.12.21"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    filevardef(fout, "lat", typeof(lat),(/"lat"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(time),(/"time"/))
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "uwnd", typeof(uwnd),(/"time","lat","lon"/))
    filevardef(fout, "vwnd", typeof(vwnd),(/"time","lat","lon"/))
   

    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "lat", lat)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "uwnd", uwnd)
    filevarattdef(fout, "vwnd", vwnd)


    fout->lon=(/lon/)
    fout->lat=(/lat/)
    fout->time=(/time/)
    fout->uwnd=(/uwnd(:,::-1,:)/)
    fout->vwnd=(/vwnd(:,::-1,:)/)


end 
   
