begin
    ;画时间序列图的时候，把BCC-CSM1-1给去掉
    models=(/"20cr","era","jra55","MERRA","ncep1","ncep2","mre"/)

    yearstart=197901
    yearend=201212  
    year1=1979
    year2=2012

;******************************************************************************************************
;循环读入数据
;******************************************************************************************************
do i=0,dimsizes(models)-1
 ; model="CanESM2"
 model=models(i)
 print("************************"+model+"*****************************")
 diri="/mnt/f/share/rea_msf/"
 diri_output="/mnt/f/share/rea_msf/"
 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_global.nc")

;;;********************************************* 
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    time_r=f1->time
    timeind=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)
    time=time_r(timeind)
    delete(time_r)
    msf_time= f1->msf_time(:,:,timeind)
    delete(f1)
    printVarSummary(msf_time)

    ;;;垂直方向插值
    print("开始垂直插值")
    z0=(/1000.0,925.0,850.0,700.0,600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0/)
    msf_time_new1=int2p_n_Wrap(msf_time&level, msf_time, z0, 1, 0)
    delete(msf_time)
    ;;;水平插值
    xo=fspan(0,357.5,144)
    msf_time_new=linint1_n_Wrap(msf_time_new1&lon,msf_time_new1,False,xo,0,1)
    delete(msf_time_new1)

    msf_time_new!0="level"
    msf_time_new!1="lon"
    msf_time_new&level@units="hpa"
    msf_time_new&level=z0
    msf_time_new&lon=xo
    msf_time_new&lon@units="degree_east"
    printVarSummary(msf_time_new)
    fi_time_new=msf_time_new
    delete(msf_time_new)

    nlev=dimsizes(fi_time_new&level)
    ntime=dimsizes(fi_time_new&time)
    nlon=dimsizes(fi_time_new&lon)

    print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    diro=diri_output
    filo1="msf_time_"+model+"_global_new.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf,w,ud经向变化的量,有时间变量1950-200512，取5到-5平均的结果"    
    fAtt@creation_date = "2018.1.2"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(fi_time_new&level), (/"level"/))
    filevardef(fout, "lon", typeof(fi_time_new&lon),(/"lon"/))
    filevardef(fout, "time", typeof(fi_time_new&time),(/"time"/))

    filevardef(fout, "msf_time", typeof(fi_time_new),(/"level","lon","time"/))

    filevarattdef(fout, "lon", fi_time_new&lon)
    filevarattdef(fout, "time", time)

    filevarattdef(fout, "msf_time", fi_time_new)
    setfileoption(fout,"DefineMode",False)

    fout->level=(/fi_time_new&level/)
    fout->lon=(/fi_time_new&lon/)
    fout->time=(/fi_time_new&time/)
    fout->msf_time=(/fi_time_new/)

    delete(time)
    delete(fi_time_new)
    delete(fout)
    delete(model)
   
end do 
end
    


    