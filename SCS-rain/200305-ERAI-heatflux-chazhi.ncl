;;********************************************************
;;计算SST倾向以及SST值与净长波辐射，短波辐射，感热，潜热以及平流和对流
;;的相关系数场，热通量数据是NCEP-NCAR数据，海洋数据是SODA2.2.4,SST是哈哈得来SST数据
;;由于数据之间的精度不同，需要进行插值预处理,SODA精度是0.5×0.5,HADISST是1×1
;;NCEP-NCAR的精度为1.87×1.9
;;
;;                                         linshuheng 2020.02.18
;**********************************************************
function chazhi(data)
 begin
   xi = data&longitude
   yi = data&latitude
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
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


    fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/heatflux/heatflux_monthly_"+"*"+".nc")
    f3 = addfiles(fils3, "r")

    date=cd_calendar(f3[:]->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude
    lon = f3[1]->longitude
   print(f3)
   exit()
    ;;;region 
      lat1 = -90
      lat2 = 90
      lon1 = 0
      lon2 = 360
    ;;;

      lhtfl = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
      shtfl = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
      nswrs = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
      nlwrs = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)

    do i = 0,dimsizes(timeind)/12-1
      lhtfl(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->slhf(:,{lat1:lat2},{lon1:lon2}))
      shtfl(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->sshf(:,{lat1:lat2},{lon1:lon2}))
      nswrs(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->ssr(:,{lat1:lat2},{lon1:lon2}))
      nlwrs(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->str(:,{lat1:lat2},{lon1:lon2}))
    end do

  


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

    outpath = "/home/ys17-19/lsh/Project/SCS-rain/ERAI-heatflux-1degree.nc"
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








    




