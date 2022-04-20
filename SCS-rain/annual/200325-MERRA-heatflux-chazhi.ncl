;;********************************************************
;;差值MERRA的通量数据
;;
;;                                         linshuheng 2020.02.18
;**********************************************************
function chazhi(data)
 begin
   xi = data&lon
   yi = data&lat
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
   data_chazhi!0 = "time"
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
   data_chazhi&time  = data&TIME
   data_chazhi&lat = yo
   data_chazhi&lon = xo
   data_chazhi&lat@units = "degrees north"
   data_chazhi&lon@units = "degrees east"
   copy_VarAtts(data, data_chazhi)
 
   return data_chazhi
end 

;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;



;;;read NCEP-NCAR heat flux data
  ;;;读入MERRA的数据
    m=""+ispan(1979,2013,1)+""

    fils=new((/35,12/),string)
    do i=0,34

    if i.le.13 then
    fils3=systemfunc("ls /home/yangsong3/data-observation/MERRA/mon/heatflux/MERRA100.prod.assim.tavgM_2d_ocn_Nx."+m(i)+"*"+".SUB.nc")
    end if
  
    if i.ge.14.and.i.le.21 then 
    fils3=systemfunc("ls /home/yangsong3/data-observation/MERRA/mon/heatflux/MERRA200.prod.assim.tavgM_2d_ocn_Nx."+m(i)+"*"+".SUB.nc")
    end if

    if i.ge.22 then 
    fils3=systemfunc("ls /home/yangsong3/data-observation/MERRA/mon/heatflux/MERRA300.prod.assim.tavgM_2d_ocn_Nx."+m(i)+"*"+".SUB.nc")
    end if
    fils(i,:)= fils3
    end do
  
  
    f3 = addfiles(fils, "r")

    date=cd_calendar(f3[:]->TIME, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->TIME(timeind)
    lat = f3[1]->lat
    lon = f3[1]->lon
    
  
    ;;;region 
      lat1 = -90
      lat2 = 90
      lon1 = 0
      lon2 = 360
    ;;;

    lhtfl = f3[:]->EFLUXWTR(:,{lat1:lat2},{lon1:lon2})
    shtfl = f3[:]->HFLUXWTR(:,{lat1:lat2},{lon1:lon2})
    nswrs = f3[:]->SWGNTWTR(:,{lat1:lat2},{lon1:lon2})
    nlwrs = f3[:]->LWGNTWTR(:,{lat1:lat2},{lon1:lon2})

  


    ;;;chazhi 
     
     lhtfl_chazhi = chazhi(lhtfl)
     shtfl_chazhi = chazhi(shtfl)

     nswrs_chazhi = chazhi(nswrs)
     nlwrs_chazhi = chazhi(nlwrs)
    

    printVarSummary(lhtfl_chazhi)

    nlat = dimsizes(lhtfl_chazhi&lat)
    nlon = dimsizes(lhtfl_chazhi&lon)
    ;;;

;;;写入nc文件

    outpath = "/home/ys17-19/lsh/Project/SCS-rain/MERRA-heatflux-1degree.nc"
       ;;-----------------------写入nc文件---------------------------
    system("rm " + outpath)   ; remove any pre-existing file
    fout = addfile(outpath ,"c")  ; open output netCDF file
    setfileoption(fout,"DefineMode",True)

   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================

    fAtt               = True            ; assign file attributes
    fAtt@title         = "1979-2013，精度为1乘1"    
    fAtt@creation_date = "2020.3.04"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    dimsize=(/-1,nlat,nlon/)
    dimunlim=(/True,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    filevardef(fout, "lat", typeof(lhtfl_chazhi&lat),(/"lat"/))
    filevardef(fout, "lon", typeof(lhtfl_chazhi&lon),(/"lon"/))
    filevardef(fout, "time", typeof(lhtfl_chazhi&time),(/"time"/))
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lhtfl",typeof(lhtfl_chazhi),(/"time","lat","lon"/))
    filevardef(fout, "shtfl",typeof(lhtfl_chazhi),(/"time","lat","lon"/))
    filevardef(fout, "nswrs",typeof(lhtfl_chazhi),(/"time","lat","lon"/))
    filevardef(fout, "nlwrs",typeof(lhtfl_chazhi),(/"time","lat","lon"/))


      ;;属性
      filevarattdef(fout, "lon", lhtfl_chazhi&lon)
      filevarattdef(fout, "lat", lhtfl_chazhi&lat)
      filevarattdef(fout, "time", time)

      filevarattdef(fout, "lhtfl", lhtfl)
      filevarattdef(fout, "shtfl", shtfl)
      filevarattdef(fout, "nswrs", nswrs)
      filevarattdef(fout, "nlwrs", nlwrs)

      ;;

      fout->lon=(/lhtfl_chazhi&lon/)
      fout->lat=(/lhtfl_chazhi&lat/)
      fout->time=(/time/)

      fout->lhtfl = (/lhtfl_chazhi/)
      fout->shtfl = (/shtfl_chazhi/)
      fout->nlwrs = (/nlwrs_chazhi/)
      fout->nswrs = (/nswrs_chazhi/)

       print("=====finish NC write==========")
    ;;;

end








    




