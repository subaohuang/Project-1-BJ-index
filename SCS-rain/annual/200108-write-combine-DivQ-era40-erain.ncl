;;********************************************************
;;将ERA40计算的DIVQ和ERAI计算的Q结合
;;数据的精度是1.5×1.5°
;;
;;linshuheng 2020.01.08
;**********************************************************
begin
  ;;;read era40 data
  diri = "/home/ys17-19/lsh/Project/SCS-rain/"
  
  f1 = addfile(diri+"era-40-DivQ.nc","r")

  date = cd_calendar(f1->time, 1)


  ;;;time span
  startmon =195701
  endmon = 197812
  timeind=ind(date.ge.startmon.and.date.le.endmon)

  divQ_era40 = f1->DivQ(timeind,:,:)

  ;;;read eraI data
  f2 = addfile(diri+"era-interium-DivQ.nc","r")

  date := cd_calendar(f2->time, 1)
  
  ;;;time span
  startmon =197901
  endmon = 201312

  timeind := ind(date.ge.startmon.and.date.le.endmon)
  divQ_eraI = f2->DivQ(timeind,:,:)
  ;;;

;;;差值，统一插值到1.9乘2.5网格,匹配模式结果
   xo=fspan(0.0, 357.5, 144)
   yo=fspan(-90,90,96)
  ;;
   
   divQ_era40_chazhi = linint2_Wrap(divQ_era40&lon,divQ_era40&lat, divQ_era40, \
                          True, xo, yo, 0)

   divQ_eraI_chazhi = linint2_Wrap(divQ_eraI&lon,divQ_eraI&lat, divQ_eraI, \
                          True, xo, yo, 0)
;;;
dim_era40 = dimsizes(divQ_era40_chazhi)
dim_eraI  = dimsizes(divQ_eraI_chazhi)

new_data = new((/dim_era40(0)+dim_eraI(0), dim_era40(1),dim_era40(2)/), typeof(divQ_era40))

new_data(0:dim_era40(0)-1,:,:) = divQ_era40_chazhi
new_data(dim_era40(0):,:,:) = divQ_eraI_chazhi

new_data!0 = "time"
new_data!1 = "lat"
new_data!2 = "lon"

new_data&time(0:dim_era40(0)-1) = divQ_era40&time 
new_data&time(dim_era40(0):) = divQ_eraI&time 

new_data&lat@units = divQ_era40&lat@units
new_data&lon@units = divQ_era40&lon@units
new_data&time@units = divQ_era40&time@units


nlat  = dimsizes(new_data&lat)
nlon = dimsizes(new_data&lon)
;;;写成nc 文件输出

outpath = "/home/ys17-19/lsh/Project/SCS-rain/era40-Interium-DivQ.nc"
   ;;-----------------------写入nc文件---------------------------
system("rm " + outpath)   ; remove any pre-existing file
fout = addfile(outpath ,"c")  ; open output netCDF file
setfileoption(fout,"DefineMode",True)

;===================================================================
; create global attributes of the file (optional)
;===================================================================

fAtt               = True            ; assign file attributes
fAtt@title         = "1957-2013年,era40和eraI精度为1.9乘2.5"    
fAtt@creation_date = systemfunc("date")
fileattdef( fout, fAtt)            ; copy file attributes    
;----------------------------坐标变量的定义
dimnames=(/"time","lat","lon"/)
dimsize=(/-1,nlat,nlon/)
dimunlim=(/True,False,False/)
filedimdef(fout, dimnames, dimsize, dimunlim)
filevardef(fout, "lat", typeof(new_data&lat),(/"lat"/))
filevardef(fout, "lon", typeof(new_data&lon),(/"lon"/))
filevardef(fout, "time", typeof(new_data&time),(/"time"/))
;-----------------------------定义变量名，类型和维数
filevardef(fout, "DivQ",typeof(new_data),(/"time","lat","lon"/))

  ;;属性
  filevarattdef(fout, "lon", new_data&lon)
  filevarattdef(fout, "lat", new_data&lat)
  filevarattdef(fout, "time", new_data&time)
  filevarattdef(fout, "DivQ", new_data)

  ;;

  fout->lon=(/new_data&lon/)
  fout->lat=(/new_data&lat/)
  fout->time=(/new_data&time/)
  fout->DivQ=(/new_data/)
   print("=====finish NC write==========")
;;;

end 



 

   



