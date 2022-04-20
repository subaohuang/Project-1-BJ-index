  ;;;************************************************************
  ;; 计算旋转风
  ;; 作为funcion调用
  ;; using function 
  ;; input values must be in ascending latitude order
  ;; input array must be on a global grid
  ;; input data must  have the dimension name of lat ,lon ,time ,lev,can use cdo setname 
  ;; to change ,the data variable must be named U,V
  ;; linshuheng 2020.11.03
  ;; 
  ;;;********************************************************、
function get_msf(path,filename,timestart,timeend,outpath,outputname)
local a, g, pi ,dfi
  begin
   
  ;;------------------读入数据---------------------------------
        diri = path   ;;; 变量
        fils1 = systemfunc("ls "+diri+filename)
     
    
        f1 = addfile(fils1, "r")
        time = f1->time
        ; time = (/time - 1/)

        date = cd_calendar(time, 1)    ;;;output calendar date 
        timeind = ind(date.ge.timestart.and.date.le.timeend)
         
         lat1 = -90
         lat2 = 90

        u = f1->U(timeind,{1000:100},{lat1:lat2},:)
        v = f1->V(timeind,{1000:100},{lat1:lat2},:)

        time = f1->time(timeind)
        delete(f1)

      

        print("计算涡度")

        vp = uv2vrF_Wrap(u, v)
        delete(u)
        delete(v)

        print("计算旋转风")

        uvr = vr2uvF_Wrap(vp)

        ur = uvr(0,:,:,:,:)
        vr = uvr(1,:,:,:,:)

        delete(uvr)
        printVarSummary(ur)
        printVarSummary(vr)
;;;;;;
      


;;;写nc文件********************************************************************

    dim=dimsizes(ur)
    ntime=dim(0)
    nlev=dim(1)
    nlat=dim(2)
    nlon=dim(3)

;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

    diro = outpath
    filo1 = outputname
    system("/bin/rm -f "+diro + filo1+"")
    ; setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "模式计算的旋转风ur,vr"    
    ; fAtt@creation_date = systemfunc("date")
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lev","lat","lon"/)
    dimsize=(/ntime,nlev, nlat,nlon/)
    dimunlim=(/False,False,False, False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数

    filevardef(fout, "time", typeof(ur&time),(/"time"/))
    filevardef(fout, "lev",typeof(ur&lev), (/"lev"/))
    filevardef(fout, "lat", typeof(ur&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(ur&lon),(/"lon"/))


    filevardef(fout, "ur", typeof(ur),(/"time","lev","lat","lon"/))
    filevardef(fout, "vr", typeof(vr),(/"time","lev","lat","lon"/))  
 
  
    filevarattdef(fout, "time", ur&time)   ;;;记得给时间附上时间单位
    filevarattdef(fout, "lev", ur&lev)
    filevarattdef(fout, "lat", ur&lat)
    filevarattdef(fout, "lon", ur&lon)

    filevarattdef(fout, "ur", ur)
    filevarattdef(fout, "vr", vr)
     
    setfileoption(fout,"DefineMode",False)
    
    fout->time=(/ur&time/)
    fout->lat=(/ur&lat/)
    fout->lon=(/ur&lon/)
    fout->lev=(/ur&lev/)

    fout->ur=(/ur/)
    fout->vr=(/vr/)
    
    echo = "rotation wind creat successfully"
return(echo)
end 

begin

;function(path,filename,timestart,timeend,outpath,outputname)
  
       ; inpath="/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmatlantic/"
       ; filename="F_1850-PDAY-rmatlan.cam.h0.197801-201212_chazhi.nc"
       ; outputpath="~/lsh/data/wc-reslut/"
       ; outputname="msf_F_1850-PDAY-rmatlan_197912-201212.nc"
  m = get_msf(inpath,filename,101,4012,outputpath,outputname)

end 

