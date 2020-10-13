    begin
    model="FGOALS-s2"
    diri_output="/mnt/j/CMIP5/history/tos_new"
    ; system("if ! test -d " + diri_output +" ; then mkdir " + diri_output + " ; fi")   ;;;make dir
    ; diri_output="/mnt/f/share/CMIP5/msf_history/"+model
    yearstart=197901
    yearend=200412

    year1=1979
    year2=2004
    fils1=systemfunc("ls /mnt/j/CMIP5/history/tos/"+model+"/r1i1p1/tos_Omon_"+model+"_historical_r1i1p1_"+"*"+".nc")
    print(fils1)
    f1=addfile(fils1, "r")
    printVarSummary(f1)

    lon=f1->lon
    lat=f1->lat
    time=f1->time

    date=cd_calendar(time, 1)
    print("******************日期*********************")


    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(timeind))
    

    sst_r=f1->tos(timeind,:,:)

    print("*****************开始水平插值************************")  
    xo=fspan(0.5,359.5,360)
    yo=fspan(-89.5, 89.5, 180)
   
    sst=linint2_Wrap(sst_r&lon, sst_r&lat, sst_r, \
                      True, xo, yo, 0)
    sst!0="time"
    sst!1="lat" 
    sst!2="lon"
    sst&lat@units="degrees_north"
    sst&lon@units="degrees_east"
    delete(sst_r)
    printVarSummary(sst)
   
    ntime=dimsizes(sst&time)
    nlat=dimsizes(sst&lat)
    nlon=dimsizes(sst&lon)
   
    print("*****************写nc文件************************")  
   ;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="tos_Omon_"+model+"_historical_r1i1p1_197901-200512.nc"
    system("/bin/rm -f "+filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "FGOALS-G2 插值 1979-200512"    
    fAtt@creation_date = "2018.4.11"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"TIME","LAT","LONN179_180"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "LAT",typeof(sst&lat), (/"LAT"/))
    filevardef(fout, "LONN179_180", typeof(sst&lon),(/"LONN179_180"/))
    filevardef(fout, "TIME", typeof(sst&time),(/"TIME"/))

   
    filevardef(fout, "TOS", typeof(sst),(/"TIME","LAT","LONN179_180"/))


    filevarattdef(fout, "LAT", sst&lat)
    filevarattdef(fout, "LONN179_180", sst&lon)
    filevarattdef(fout, "TIME", sst&time)
    filevarattdef(fout, "TOS", sst)
 
   
    setfileoption(fout,"DefineMode",False)

    fout->LAT=(/sst&lat/)
    fout->LONN179_180=(/sst&lon/)
    fout->TIME=(/sst&time/)
    fout->TOS=(/sst/)
 end

    