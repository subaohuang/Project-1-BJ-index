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
   xi = data&lon
   yi = data&lat
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
   data_chazhi&lat = yo
   data_chazhi&lon = xo
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

    path = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/hflx/"
    var_name = (/"lhtfl","shtfl","nswrs","nlwrs"/)
    filein = path + var_name +".sfc.mon.mean.nc"
    
    ;;;region 
      lat1 = -90
      lat2 = 90
      lon1 = 0
      lon2 = 360
    ;;;

   
    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    lhtfl = f1->lhtfl(timeind,{lat1:lat2},{lon1:lon2})

    f2 = addfile(filein(1),"r")
    shtfl = f2->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    
    f3 = addfile(filein(2),"r")
    nswrs = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})  
    
    f4 = addfile(filein(3),"r")
    nlwrs = f4->nlwrs(timeind,{lat1:lat2},{lon1:lon2})  
    

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

    outpath = "/home/ys17-19/lsh/Project/SCS-rain/NCEP-NCAR-heatflux-1degree.nc"
       ;;-----------------------写入nc文件---------------------------
    system("rm " + outpath)   ; remove any pre-existing file
    fout = addfile(outpath ,"c")  ; open output netCDF file
    setfileoption(fout,"DefineMode",True)

   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================

    fAtt               = True            ; assign file attributes
    fAtt@title         = "1979-2013，精度为1乘1"    
    fAtt@creation_date = "2020.2.18"
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
      filevarattdef(fout, "time", lhtfl_chazhi&time)

      filevarattdef(fout, "lhtfl", lhtfl_chazhi)
      filevarattdef(fout, "shtfl", shtfl_chazhi)
      filevarattdef(fout, "nswrs", nswrs_chazhi)
      filevarattdef(fout, "nlwrs", nlwrs_chazhi)

      ;;

      fout->lon=(/lhtfl_chazhi&lon/)
      fout->lat=(/lhtfl_chazhi&lat/)
      fout->time=(/lhtfl_chazhi&time/)

      fout->lhtfl = (/lhtfl_chazhi/)
      fout->shtfl = (/shtfl_chazhi/)
      fout->nlwrs = (/nlwrs_chazhi/)
      fout->nswrs = (/nswrs_chazhi/)

       print("=====finish NC write==========")
    ;;;

end








    




