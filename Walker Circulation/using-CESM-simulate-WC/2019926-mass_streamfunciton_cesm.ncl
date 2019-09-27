  ;;;************************************************************
  ;; 计算质量流函数
  ;; 作为funcion调用计算MSF
  ;; using function 
  ;; input values must be in ascending latitude order
  ;; input array must be on a global grid
  ;; input data must  have the dimension name of lat ,lon ,time ,lev,can use cdo setname 
  ;; to change ,the data variable must be named U,V,OMEGA
  ;; dv2uvF_Wrap
  ;; linshuheng 2019.9.25
  ;; 
  ;;;********************************************************、
function get_msf(path,filename,timestart,timeend,outpath,outputname)
local a, g, pi ,dfi
  begin
   ;;;流函数公式需要的参数
    a = 6371000
    g = 9.8
    pi = get_pi("double")
    dfi = (pi*10)/180
  ;;------------------读入数据---------------------------------
        diri = path   ;;; 变量 
        fils1 = systemfunc("ls "+diri+filename)
     
    
        f1 = addfile(fils1, "r")
        time = f1->time
        lat = f1->lat
        lon = f1->lon
        lev= f1->lev   ;level反着取
        levelw=f1->lev

        lev@units = "hpa"
        lev@long_name = "pressure level"
        
        

        date = cd_calendar(time, 1)    ;;;output calendar date 

        timeind = ind(date.ge.timestart.and.date.le.timeend)
        time_wr = time(timeind)
        
        ;;;select region
        latind1 = ind_nearest_coord(5,lat, 0)
        latind2 = ind_nearest_coord(-5,lat, 0)

        ; lonind1 = ind(lon.eq.0)
        ; lonind2 = ind(lon.eq.357.5)

        lonind1 = ind_nearest_coord(0, lon,0)
        lonind2 = ind_nearest_coord(360,lon,0)

        levindw1 = ind(levelw.eq.100)
        levindw2 = ind(levelw.eq.1000)

        uwnd = f1->U(timeind,:,:,:)   ;;;variable must have ascending lat and lon order coordinates
        vwnd = f1->V(timeind,:,:,:)   ;;;variable must have ascending lat and lon order coordinates
         
        w = f1->OMEGA(timeind,:,latind2:latind1,:)   ;;;ensure the lat order of data from south to north 

        
        w_avg = dim_avg_n_Wrap(w(time|:,lat|:,lev|:,lon|:) , (/0,1/)) 
        w_time = dim_avg_n_Wrap(w(time|:,lat|:,lev|:,lon|:), 1)
        w_time_new = w_time*(-50)
        copy_VarMeta(w_time, w_time_new)

        ;;; 
   ;; 

   ;;; -------------Calcaulate divergence wind -------------------------------
    
        vrdvx = uv2vrdvF(uwnd,vwnd)
        uvd = dv2uvF(vrdvx(1,:,:,:,:))
        delete(vrdvx)
        ud = uvd(0,:,:,:,:)
        delete(uvd)
        ud!0 = "time"
        ud!1 = "level"
        ud!2 = "lat"
        ud!3 = "lon"
        printVarSummary(ud)
        ud_new = ud(time|:,lat|latind2:latind1,lon|:,level|:)  
        delete(ud)

        ; ud_new&lon = w&lon
        ; ud_new&lat = fspan(-5,5,dimsizes(w&lat))
        ; ud_new&time = ispan(1,dimsizes(w&time),1)
        ; ud_new&level = level

        copy_VarCoords(w(time|:,lat|:,lon|:,lev|:), ud_new)
    
        ud_avg = dim_avg_n_Wrap(ud_new,(/0,1/))    ;;; timemean and meridional mean 
        ud_avg_new = ud_avg(lev|0:dimsizes(levelw)-1,lon|:)  ;;


        ud_time = dim_avg_n_Wrap(ud_new(time|:,lat|:,lev|0:dimsizes(levelw)-1,lon|:), 1)
        ud_time_new = ud_time(lev|:,lon|:,time|:)

        w_avg_new = w_avg*(-50)
        copy_VarCoords(w_avg, w_avg_new)
    ;;;; 

    ;;;---------Calcaulate mass stream function-----------------------------
   
        pbot = w_avg_new&lev
        pres_new = new((/dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)
        pres_new = 1013.0
        vint = new((/dimsizes(w_avg_new&lev),dimsizes(ud_new&time),dimsizes(ud_new&lat),dimsizes(ud_new&lon)/),float)


        do i = 0,(dimsizes(w_avg_new&lev)-1)

          vint(i,:,:,:) = vibeta(lev, ud_new, 1, pres_new,pbot(i),0)    ;;;intergrate from specific level to 0 ,the order of ud_new's level must 
                                                                          ;;; Bottom to top
        end do
        
        fi = ((a*dfi)/g)*vint*100
        delete(vint)

        fi_mean = dim_avg_n_Wrap(fi, (/1,2/))     
    ;;;;
    

        fi_mean!0 = "lev"
        fi_mean!1 = "lon"
        fi_mean&lev = pbot
        fi_mean&lev@units = "hpa"
        fi_mean&lon = ud_new&lon

        fi_time = dim_avg_n_Wrap(fi, 2)
        delete(fi)
        fi_time!0 = "lev"
        fi_time!1 = "time"
        fi_time!2 = "lon"
        fi_time&lev = pbot
        fi_time&lev@units = "hpa"
        fi_time&lon = ud_new&lon
        fi_time&time = ud_new&time
        fi_time_new = fi_time(time|:,lev|:,lon|:)

        ntime=dimsizes(fi_time_new&time)
        nlev=dimsizes(fi_time_new&lev)
        nlon=dimsizes(fi_time_new&lon)
 
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

        diro = outpath
        filo1 = outputname
        system("/bin/rm -f "+diro + filo1+"")
        fout=addfile(diro + filo1,"c")
        setfileoption(fout,"DefineMode",True)
        ;======================================文件名
        fAtt               = True            ; assign file attributes
        fAtt@title         = "MSF of CESM "    
        fAtt@creation_date = systemfunc("date")
        fileattdef( fout, fAtt)            ; copy file attributes    
        ;----------------------------坐标变量的定义
        dimnames=(/"time","lev","lon"/)
        dimsize=(/ntime,nlev,nlon/)
        dimunlim=(/False,False,False/)
        filedimdef(fout, dimnames, dimsize, dimunlim)
        ;-----------------------------定义变量名，类型和维数
        filevardef(fout, "lev",typeof(lev), (/"lev"/))
        filevardef(fout, "lon", typeof(lon),(/"lon"/))
        filevardef(fout, "time", typeof(time),(/"time"/))

        filevardef(fout, "msf_time", typeof(fi_time_new),(/"time","lev","lon"/))
        filevardef(fout, "ud_time", typeof(ud_time_new),(/"time","lev","lon"/))
        filevardef(fout, "w_time", typeof(w_time_new),(/"time","lev","lon"/))

        filevarattdef(fout, "lev", lev)
        filevarattdef(fout, "lon", lon)
        filevarattdef(fout, "time", time_wr)

      
        filevarattdef(fout, "msf_time", fi_time_new)
        filevarattdef(fout, "ud_time", ud_time_new)
        filevarattdef(fout, "w_time", w_time_new)


        setfileoption(fout,"DefineMode",False)

        fout->lev=(/fi_mean&lev/)
        fout->lon=(/fi_mean&lon/)
        fout->time=(/time_wr/)

        fout->msf_time=(/fi_time_new/)
        fout->ud_time=(/ud_time_new(time|:,lev|:,lon|:)/)
        fout->w_time=(/w_time_new/)
    echo = "msf creat successfully"
return(echo)
end 

begin

;function(path,filename,timestart,timeend,outpath,outputname)
  path = "~/lsh/CESM-data/"
  filename = "f19_f19_FAMIP.cam.h0.197901-200611_chazhi.nc"
  m = get_msf(path,filename,197901,201212,"~/lsh/data/wc-reslut/","msf_f19-f19-FAMIP_1979-200611.nc")

end 

