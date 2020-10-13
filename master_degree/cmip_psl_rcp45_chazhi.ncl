  begin
;;读入数据
   
   
   model=(/"INMCM4","IPSL-CM5A-MR",\
   "MIROC5","MPI-ESM-MR"/)

   do i=0,dimsizes(model)-1 
    diri_output="/mnt/i/CMIP5/RCP45/slp_chazhi/"
 
    yearstart=200501
    yearend=201212
     
    
    if model(i).eq."FGOALS-s2".or.model(i).eq."HadGEM2-ES".or.model(i).eq."HadGEM2-CC"   then
    yearstart=197901
    yearend=200412
    end if

    fils1=systemfunc("ls /mnt/i/CMIP5/RCP45/psl/"+model(i)+"/r1i1p1/psl_Amon_"+model(i)+"_rcp45_r1i1p1_"+"*"+".nc")
    print(fils1)
    f1=addfiles(fils1, "r")
    delete(fils1)
    time=f1[:]->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    delete(date)
  
      psl_r=lonFlip(f1[:]->psl(timeind,:,:))
      psl_r!0="time"
      psl_r!1="lat"
      psl_r!2="lon"


    xo=fspan(-180,177.5,144)
    yo=fspan(-87.5, 87.5, 71)
    
     print("*****************开始水平插值************************")  
    psl=linint2_Wrap(psl_r&lon, psl_r&lat, psl_r, \
                      True, xo, yo, 0)


     psl!1="lat"
     psl!2="lon"
    
    psl&lat@units="degrees_north"
    psl&lon@units="degrees_east"

    delete(psl_r)


    printVarSummary(psl)
    
    ntime=dimsizes(psl&time)
    nlat=dimsizes(psl&lat)
    nlon=dimsizes(psl&lon)
  
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="psl_"+model(i)+"_rcp45_r1i1p1_200601-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "psl插值后的结果，取1000hpa"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(psl&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(psl&lon),(/"lon"/))
    filevardef(fout, "time", typeof(psl&time),(/"time"/))

    filevardef(fout, "psl", typeof(psl),(/"time","lat","lon"/))



    filevarattdef(fout, "lat", psl&lat)
    filevarattdef(fout, "lon", psl&lon)
    filevarattdef(fout, "time", time)
    
    filevarattdef(fout, "psl", psl)

    
    setfileoption(fout,"DefineMode",False)

    fout->lat=(/psl&lat/)
    fout->lon=(/psl&lon/)
    fout->time=(/psl&time/)
    fout->psl=(/psl/)

    delete(psl)

    delete(ntime)
    delete(nlon)
    delete(nlat)
    delete(time)
 end do

end 



