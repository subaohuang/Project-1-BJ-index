;;*************************************************************
;;;;MRE_DATA  combine: slp
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
    fils1="/mnt/f/share/ncep/ncep1/slp.mon.mean.nc"
    f1=addfile(fils1, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)
    time=f1->time(timeind)
    slp_r=f1->slp(timeind,::-1,:)   ;;;插值前先逆着取
  
    print("开始水平插值")  
    slp_necp1=linint2_Wrap(slp_r&lon, slp_r&lat, slp_r, True, xo, yo, 0)
  
    slp_necp1!1="lat"
    slp_necp1!2="lon"

    lat=slp_necp1&lat(::-1)
    lon=slp_necp1&lon
    lat@units="degree north"
    lon@units="degree east"

    delete(slp_r)
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)

;;;

print("**************20cr******************************")
;;;20CR
    fils1="/mnt/f/share/ncep/20CR/prmsl.mon.mean.nc"
    
    f1=addfile(fils1, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)

    slp_r=f1->prmsl(timeind,::-1,:)
    slp_r=(/slp_r/100.0/)
    delete(f1)
    delete(fils1)

    ;;;interpolation

    print("开始水平插值")  
    slp_20cr=linint2_Wrap(slp_r&lon, slp_r&lat, slp_r, True, xo, yo, 0)
    delete(slp_r)

    slp_20cr!1="lat"
    slp_20cr!2="lon"  
    printVarSummary(slp_20cr)
    ;;;
  
    delete(date)
    delete(timeind)


;;;
    print("*******************jra55***********************************")
    ;;;JRA-55
    fils1="/mnt/f/share/JRA-55/prmsl.mon.mean.nc"

    f1=addfile(fils1, "r")

    slp_r=f1->prmsl(252:659,::-1,:)  ;;;195801-201212
    slp_r=(/slp_r/100.0/)
    delete(f1)

    print("开始水平插值")  
    slp_jra55=linint2_Wrap(slp_r&lon, slp_r&lat, slp_r, True, xo, yo, 0)
    slp_jra55!1="lat"
    slp_jra55!2="lon"

    delete(slp_r)
    printVarSummary(slp_jra55)
    delete(fils1)

;;;

    ;;;读入era_in的数据
    fils1= "/mnt/f/share/era/data/data_series/era_sur_global.nc"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212) 
    slp_r=f1->msl(timeind,::-1,:)
    
    
    slp_in=dble2flt(slp_r*slp_r@scale_factor+slp_r@add_offset)
    copy_VarCoords(slp_r, slp_in)
    slp_in=(/slp_in/100.0/)
    delete(slp_r)
    dnames = getVarDimNames(slp_in)
  
    if any(dnames.eq."latitude")
    slp_in!1="lat"
    slp_in!2="lon"
    end if

    print("开始水平插值")  
    slp_era=linint2_Wrap(slp_in&lon, slp_in&lat, slp_in, True, xo, yo, 0)
    slp_era!1="lat"
    slp_era!2="lon"
    delete(slp_in)

    delete(fils1)
    delete(date)
    delete(timeind)
;;;

;;;
    print("*********************ncep2*****************************")
    ;;NCEP1 U,V 
    fils1="/mnt/f/share/ncep/ncep2/mslp.mon.mean.nc"
  
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201212)

    time=f1->time(timeind)
    slp_r=f1->mslp(timeind,::-1,:)   ;;;插值前先逆着取
  
    slp_n=slp_r*slp_r@scale_factor+slp_r@add_offset
    slp_n=(/slp_n/100.0/)
    copy_VarCoords(slp_r, slp_n)
    delete(slp_r)
  
    print("开始水平插值")  
    slp_necp2=linint2_Wrap(slp_n&lon, slp_n&lat, slp_n, True, xo, yo, 0)
 
    slp_necp2!1="lat"
    slp_necp2!2="lon"
  
    delete(slp_n)
    delete(f1)
    delete(timeind)
    delete(date)
    delete(fils1)
 
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
    slp_r=lonPivot(f1[:]->SLP(timeind,::-1,:),0)
    slp_r=(/slp_r/100.0/)
    print("开始水平插值")  
    slp_merra=linint2_Wrap(slp_r&lon, slp_r&lat, slp_r, True, xo, yo, 0)
 
    slp_merra!1="lat"
    slp_merra!2="lon"

    delete(slp_r)
    delete(f1)
    delete(timeind)
    delete(date)
   


;***********************************************************************
;data_mean
;************************************************************************
    slp=slp_necp1

    slp=(/(slp_necp1+slp_20cr+slp_era+slp_jra55+slp_merra+slp_necp2)/6.0/)

   delete(slp_necp1)
   delete(slp_20cr)
   delete(slp_era)
   delete(slp_jra55)
   delete(slp_merra)
   delete(slp_necp2)

  print("write nc")
  ntime=dimsizes(time)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)

;***********************************************************************
;write nc
;***********************************************************************

    diro=diri_out
    filo1="mre_slp_chazhi.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "SLP的1979-2012,era40,ncep1,jra-55,20cr,merra,ncep1的风场平均结果，精度为2.5乘2.5"    
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
    filevardef(fout, "slp", typeof(slp),(/"time","lat","lon"/))
   

    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "lat", lat)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "slp", slp)

    fout->lon=(/lon/)
    fout->lat=(/lat/)
    fout->time=(/time/)
    fout->slp=(/slp(:,::-1,:)/)



end 
   
