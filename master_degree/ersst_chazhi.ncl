
 begin

    f  =addfile("/mnt/f/share/ersst/sst.mnmean.v4.nc", "r")
    yearstart=197901
    yearend=201212  
    date=cd_calendar(f->time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    sst_r=f->sst(timeind,::-1,:)
      ;;;插值 sst 

    xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)

     print("*****************开始水平插值************************")  
     sst=linint2_Wrap(sst_r&lon, sst_r&lat, sst_r, \
                      True, xo, yo, 0)

     sst!1="lat"
     sst!2="lon"
     
     sst&lat@units="degrees_north"
     sst&lon@units="degrees_east"

     delete(sst_r)
     printVarSummary(sst)

    ntime=dimsizes(sst&time)
    nlon=dimsizes(sst&lon)
    nlat=dimsizes(sst&lat)

 

;***************************************************
;写入nc文件
;*************************************************
    diro="/mnt/f/share/ersst/"
    filo1="ersst_v4b.chazhi.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "ersst 插值到2.5格点上 插值的结果"    
    fAtt@creation_date = "2018.12.7"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
   
    filevardef(fout, "time", typeof(sst&time),(/"time"/))
    filevardef(fout, "lat", typeof(sst&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(sst&lon),(/"lon"/))
    filevardef(fout, "sst", typeof(sst),(/"time","lat","lon"/))
 
     filevarattdef(fout, "time", sst&time)
     filevarattdef(fout, "lat", sst&lat)
     filevarattdef(fout, "lon", sst&lon)
    setfileoption(fout,"DefineMode",False)

    fout->time=(/sst&time/)
    fout->lat=(/sst&lat/)
    fout->lon=(/sst&lon/)
    fout->sst=(/sst/)
   delete(sst)
   delete(filo1)
   delete(ntime)
   delete(nlat)
   delete(nlon) 

end 