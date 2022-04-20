  begin
;;读入数据

  
   model=(/"CSIRO-Mk3-6-0","GFDL-CM3","INMCM4","IPSL-CM5A-MR",\
   "MIROC5","MPI-ESM-MR"/)


 ; model="Ipr-CM5A-MR"
    do i=0,dimsizes(model)-1
    diri_output="/mnt/i/CMIP5/rcp45/pr_chazhi/"
 
    yearstart=200501
    yearend=201212
     

    if model(i).eq."FGOALS-s2".or.model(i).eq."HadGEM2-ES".or.model(i).eq."HadGEM2-CC"   then
    yearstart=197901
    yearend=200412
    end if


    fils1=systemfunc("ls /mnt/i/CMIP5/rcp45/pr/"+model(i)+"/r1i1p1/pr_Amon_"+model(i)+"_rcp45_r1i1p1_"+"*"+".nc")
    print(fils1)
    f1=addfiles(fils1, "r")
    delete(fils1)

    time=f1[:]->time
    date=cd_calendar(time, 1)
    print("******************日期*********************")
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    delete(date)
   
    pr_r=f1[:]->pr(timeind,:,:)
    pr_r!0="time"
    pr_r!1="lat"
    pr_r!2="lon"
    
    xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)
    
     print("*****************开始水平插值************************")  
    pr=linint2_Wrap(pr_r&lon, pr_r&lat, pr_r, \
                      True, xo, yo, 0)


      pr!1="lat"
      pr!2="lon"

      pr&lat@units="degrees_north"
      pr&lon@units="degrees_east"

    delete(pr_r)


    printVarSummary(pr)
    
    ntime=dimsizes(pr&time)
    nlat=dimsizes(pr&lat)
    nlon=dimsizes(pr&lon)
  
    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="pr_"+model(i)+"_rcp45_r1i1p1_200601-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "pr插值后的结果，取1000hpa"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(pr&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(pr&lon),(/"lon"/))
    filevardef(fout, "time", typeof(pr&time),(/"time"/))
    filevardef(fout, "pr", typeof(pr),(/"time","lat","lon"/))



    filevarattdef(fout, "lat", pr&lat)
    filevarattdef(fout, "lon", pr&lon)
    filevarattdef(fout, "time", time)
    filevarattdef(fout, "pr", pr)

    
    setfileoption(fout,"DefineMode",False)

    fout->lat=(/pr&lat/)
    fout->lon=(/pr&lon/)
    fout->time=(/pr&time/)
    fout->pr=(/pr/)

    delete(pr)

    delete(ntime)
    delete(nlon)
    delete(nlat)
    delete(time)
 end do

end 



