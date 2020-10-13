  begin
;;读入数据
   
    ; model="ACCESS1-0"
    ; model="CanESM2"
    model="FGOALS-s2"
    diri_output="/mnt/f/share/CMIP5/AMIP/ts_new/"
   ; system("if ! test -d " + diri_output +" ; then mkdir " + diri_output + " ; fi")   ;;;make dir
   ; diri_output="/mnt/f/share/CMIP5/msf_history/"+model
   yearstart=197901
   yearend=200812
  
    fils1=systemfunc("ls /mnt/f/share/CMIP5/AMIP/ts/"+model+"/r1i1p1/ts_Amon_"+model+"_amip_r1i1p1"+"_197901-200812"+".nc")

    print(fils1)

    f1=addfiles(fils1, "join")
    time=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date=cd_calendar(time, 1)
    timeind=ind(date.le.yearend.and.date.ge.yearstart)


    uwnd_r=f1[:]->ts(timeind,:,:)  ;;;195801-201212
    delete(f1)


    xo=fspan(0,357.5,144)
    yo=fspan(-87.5, 87.5, 71)

     print("*****************开始水平插值************************")  
    uwnd=linint2_Wrap(uwnd_r&lon, uwnd_r&lat, uwnd_r, \
                      True, xo, yo, 0)

     uwnd!1="lat"
     uwnd!2="lon"
     
   uwnd&lat@units="degrees_north"
   uwnd&lon@units="degrees_east"

   delete(uwnd_r)
   printVarSummary(uwnd)

    ntime=dimsizes(uwnd&time)
    nlat=dimsizes(uwnd&lat)
    nlon=dimsizes(uwnd&lon)

    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="ts_"+model+"_amip_r1i1p1_197901-200812.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "插值到格点上的结果"    
    fAtt@creation_date = "2018.11.8"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(uwnd&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(uwnd&lon),(/"lon"/))
    filevardef(fout, "time", typeof(uwnd&time),(/"time"/))

    filevardef(fout, "ts", typeof(uwnd),(/"time","lat","lon"/))
  

    filevarattdef(fout, "lat", uwnd&lat)
    filevarattdef(fout, "lon", uwnd&lon)
    filevarattdef(fout, "time", uwnd&time)

    filevarattdef(fout, "ts", uwnd)

    setfileoption(fout,"DefineMode",False)

    fout->lat=(/uwnd&lat/)
    fout->lon=(/uwnd&lon/)
    fout->time=(/uwnd&time/)
    fout->ts=(/uwnd/)


end 


