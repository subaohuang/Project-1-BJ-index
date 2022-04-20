;;;************************************************************
;;
;; linshuheng 2021.03.04
;; 计算质量流函数
;;;*********************************************************
begin
;;;region and time period 
   timestart  = 197901
   timeend  =202012
;;;;
 
  a=6371000
  g=9.8
  pi=get_pi("double")
  dfi=(pi*10)/180

;;;读入数据
  diri="/home/yangsong3/data-observation/ERA5-monthly/pressure/"
  fils1 = systemfunc("ls "+diri+"u_mon_1x1_1979_2020.nc")
  fils2 = systemfunc("ls "+diri+"v_mon_1x1_1979_2020.nc")
  fils3 = systemfunc("ls "+diri+"w_mon_1x1_1979_2020.nc")
  

  f1 = addfile(fils1, "r")
  f2 = addfile(fils2, "r")
  f3 = addfile(fils3, "r")
 

  time = f1->time
  lat = f1->latitude
  lon = f1->longitude
  lev = f1->level
  levwant= (/1000.0, 925.0, 850.0, 700.0, 600.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0,50.0,10.0/)   ;level反着取
  
  levind1 = get1Dindex(lev,levwant)
  
  level  = f1->level(levind1)


  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)
  time_wr = time(timeind)


  ;;;select region
  latind1 = ind_nearest_coord(5,lat, 0)
  latind2 = ind_nearest_coord(-5,lat, 0)

  ; lonind1 = ind(lon.eq.0)
  ; lonind2 = ind(lon.eq.357.5)

  lonind1 = ind_nearest_coord(0, lon,0)
  lonind2 = ind_nearest_coord(360,lon,0)

  levw = f3->level
  levind2 = get1Dindex(levw,levwant)

  w = short2flt(f3->w(timeind,levind2,latind2:latind1,:))   ;;;ensure the lat order of data from south to north 
  w_avg = dim_avg_n_Wrap(w(time|:,latitude|:,level|:,longitude|:) , (/0,1/)) 
  w_time = dim_avg_n_Wrap(w(time|:,latitude|:,level|:,longitude|:), 1)
  w_time_new = w_time*(-50)
  copy_VarMeta(w_time, w_time_new)
  


;;; 
;; 

   ;;; -------------Calcaulate divergence wind -------------------------------

      uwnd = short2flt(f1->u(timeind,levind1,::-1,:))   ;;;variable must have ascending lat and lon order coordinates
      vwnd = short2flt(f2->v(timeind,levind1,::-1,:))   ;;;variable must have ascending lat and lon order coordinates
      printVarSummary(uwnd)
    
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

      copy_VarCoords(w(time|:,latitude|:,longitude|:,level|:), ud_new)
      printVarSummary(ud_new)
      print(ud_new&level)
  
   
      ud_time = dim_avg_n_Wrap(ud_new(time|:,latitude|:,level|:,longitude|:), 1)
      ud_time_new = ud_time(level|:,longitude|:,time|:)

      w_avg_new = w_avg*(-50)
      copy_VarCoords(w_avg, w_avg_new)
    ;;;; 

    ;;;---------Calcaulate mass stream function-----------------------------
   
        pbot = w_avg_new&level
        pres_new = new((/dimsizes(ud_new&time),dimsizes(ud_new&latitude),dimsizes(ud_new&longitude)/),float)
        pres_new = 1013.0
        vint = new((/dimsizes(w_avg_new&level),dimsizes(ud_new&time),dimsizes(ud_new&latitude),dimsizes(ud_new&longitude)/),float)


        do i = 0,(dimsizes(w_avg_new&level)-1)

          vint(i,:,:,:) = vibeta(levwant, ud_new, 1, pres_new,pbot(i),0)    ;;;intergrate from specific level to 0 ,the order of ud_new's level must 
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
        fi_mean&lon = ud_new&longitude

        fi_time = dim_avg_n_Wrap(fi, 2)
        delete(fi)
        fi_time!0 = "lev"
        fi_time!1 = "time"
        fi_time!2 = "lon"
        fi_time&lev = pbot
        fi_time&lev@units = "hpa"
        fi_time&lon = ud_new&longitude
        fi_time&time = ud_new&time
        fi_time_new = fi_time(time|:,lev|:,lon|:)

        ntime=dimsizes(fi_time_new&time)
        nlev=dimsizes(fi_time_new&lev)
        nlon=dimsizes(fi_time_new&lon)

;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

        diro = "/home/yangsong3/data-observation/linshh/data/rea_msf/"
        filo1 = "msf_time_era5_global_2020.nc"
        system("/bin/rm -f "+diro + filo1+"")
        ; setfileoption("nc","Format","LargeFile")
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
        fout->ud_time=(/ud_time_new(time|:,level|:,longitude|:)/)
        fout->w_time=(/w_time_new/)
print("msf creat successfully")

end 




 