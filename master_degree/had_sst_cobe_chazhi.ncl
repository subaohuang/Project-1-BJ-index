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

;;读入数据
    print("*****************HadISST*********************")
    yearstart=197901
    yearend=201212
    
    year1=1979
    year2=2012
    

    fils1=systemfunc("ls /mnt/f/share/ersst/HadISST_sst.nc")
    f1=addfile(fils1, "r")
    

    vNames= getfilevarnames(f1) ; get variable names of all groups on file
    print(vNames)
    
    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    sst_r_had=lonFlip(f1->sst(timeind,::-1,:))
    dNames=getVarDimNames(sst_r_had)
  
    if any(dNames.eq."latitude") then 
    sst_r_had!0="time"
    sst_r_had!1="lat"
    sst_r_had!2="lon"
    end if 
    

    printVarSummary(sst_r_had)
    
    delete(sst_r_had&time)

 
    print("*****************开始水平插值************************")  
    sst_had=linint2_Wrap(sst_r_had&lon, sst_r_had&lat, sst_r_had, \
                      True, xo, yo, 0)


     sst_had!1="lat"
     sst_had!2="lon"
    
    sst_had&lat@units="degrees_north"
    sst_had&lon@units="degrees_east"

    delete(sst_r_had)
    delete(date)
    delete(f1)
    delete(fils1)
    printVarSummary(sst_had)    

  ;;;*******************************************************************
    print("*****************ERSST*********************")
    ;;;读入ERSST海温数据
    m=""+ispan(1979,2012,1)+""
    diri="/mnt/f/share/ersst/"

    fils=new((/34,12/),string)
    do i=0,33
    fils1=systemfunc("ls "+diri+"ersst."+m(i)+"*"+".nc")
    fils(i,:)= fils1
    end do
    fils_new=onedtond(fils, 34*12)
    fileh=new(34*12,string)
    fileh(0:34*12-1)=fils_new

    f1=addfiles(fileh, "join")
    printVarSummary(f1)
    date=cd_calendar(f1[:]->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    sst_r_ersst=f1[:]->sst(timeind,0,:,:)
    ; sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r_ersst)

    delete(sst_r_ersst&time)


    print("*****************开始水平插值************************")  
    sst_ersst=linint2_Wrap(sst_r_ersst&lon, sst_r_ersst&lat, sst_r_ersst, \
                      True, xo, yo, 0)


     sst_ersst!1="lat"
     sst_ersst!2="lon"
    
    sst_ersst&lat@units="degrees_north"
    sst_ersst&lon@units="degrees_east"

    delete(sst_r_ersst)
    delete(date)
    delete(f1)
    delete(fils1)
    printVarSummary(sst_ersst)    

 print("*****************COBE*********************")

    fils1=systemfunc("ls /mnt/f/share/ersst/COBE_SST.nc")
    f1=addfile(fils1, "r")

    date=cd_calendar(f1->time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    sst_r_cobe=f1->sst(timeind,::-1,:)
    ; sst_r=(/sst_r-273.15/)
    printVarSummary(sst_r_cobe)
    
    delete(sst_r_cobe&time)
   
    print("*****************开始水平插值************************")  
    sst_cobe=linint2_Wrap(sst_r_cobe&lon, sst_r_cobe&lat, sst_r_cobe, \
                      True, xo, yo, 0)


     sst_cobe!1="lat"
     sst_cobe!2="lon"
    
    sst_cobe&lat@units="degrees_north"
    sst_cobe&lon@units="degrees_east"
    
    time=f1->time(timeind)
    lat=sst_cobe&lat(::-1)
    lon=sst_cobe&lon
    print(lon)
 
    delete(sst_r_cobe)
    delete(date)
    delete(f1)
    delete(fils1)
    printVarSummary(sst_cobe)    

    sst_mean=(/sst_had+sst_cobe+sst_ersst/)/3.0
    copy_VarCoords(sst_cobe, sst_mean)
    
  print("write nc")
  ntime=dimsizes(time)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)

;***********************************************************************
;write nc
;***********************************************************************

    diro=diri_out
    filo1="haderscobe_sst_chazhi.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "Sst的海温场平均的结果，精度为2.5乘2.5"    
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
    filevardef(fout, "sst", typeof(sst_mean),(/"time","lat","lon"/))
   

    filevarattdef(fout, "lon", lon)
    filevarattdef(fout, "lat", lat)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "sst", sst_mean)

    fout->lon=(/lon/)
    fout->lat=(/lat/)
    fout->time=(/time/)
    fout->sst=(/sst_mean(:,::-1,:)/)



end 
   
