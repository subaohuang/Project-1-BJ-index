begin
; f1=addfile("/mnt/i/CMIP5/history/tos_new/BNU-ESM/r1i1p1/tos_Omon_BNU-ESM_historical_r1i1p1_185001-200512.nc", "r")
f1=addfile("/mnt/i/CMIP5/RCP45/tos_new/BNU-ESM/r1i1p1/tos_Omon_BNU-ESM_rcp45_r1i1p1_200601-210012.nc", "r")
lat2d=f1->lat
lon2d=f1->lon
time=f1->time
date=cd_calendar(time, 1)
; timeind=ind(date.ge.197901.and.date.le.200512)
timeind=ind(date.ge.200601.and.date.le.201212)
print(date(timeind))
tos=f1->tos(timeind,:,:)
printVarSummary(lat2d)
printVarSummary(lon2d)
print(lat2d(:,1))
print(lon2d(0,:))
; exit()
delete(tos&j)
delete(tos&i)
tos!1="lat"
tos!2="lon"
tos&lat=lat2d(:,1)
tos&lon=lon2d(0,:)
printVarSummary(tos)
tos_r=lonFlip(tos)
printVarSummary(tos_r)
;;;;插值要插值到-179.5-179.5 精度1度
;;;LAT要插值到 -89.5 - 89.5 1度
xo=fspan(-179.5,179.5,360)
yo=fspan(-89.5, 89.5, 180)

print("*****************开始水平插值************************")  
TOS=linint2_Wrap(tos_r&lon,tos_r&lat, tos_r, \
              True, xo, yo, 0)
TOS!0="TIME"
TOS!1="LAT"
TOS!2="LONN179_180"

TOS&LAT@units="degrees_north"
TOS&LONN179_180@units="degrees_east"
printVarSummary(TOS)

dim=dimsizes(TOS)
ntime=dim(0)
nlat=dim(1)
nlon=dim(2)

print("***************开始写nc文件****************************")
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
    ; diro="/mnt/i/CMIP5/history/tos_new/BNU-ESM/r1i1p1/"
    diro="/mnt/i/CMIP5/RCP45/tos_new/BNU-ESM/r1i1p1/"
    filo1="tos_Omon_BNU-ESM_rcp45_r1i1p1_200601-201212.nc"
    system("/bin/rm -f "+diro + filo1+"")
    setfileoption("nc","Format","LargeFile")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "插值后的BNU-ESM"    
    fAtt@creation_date = "2019.2.21"
    fileattdef( fout, fAtt)            ; copys file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"TIME","LAT","LONN179_180"/)
    dimsize=(/ntime,nlat,nlon/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "LAT",typeof(TOS&LAT), (/"LAT"/))
    filevardef(fout, "LONN179_180", typeof(TOS&LONN179_180),(/"LONN179_180"/))
    filevardef(fout, "TIME", typeof(TOS&TIME),(/"TIME"/))
    filevardef(fout, "TOS", typeof(TOS),(/"TIME","LAT","LONN179_180"/))



    filevarattdef(fout, "LAT", TOS&LAT)
    filevarattdef(fout, "LONN179_180",TOS&LONN179_180)
    filevarattdef(fout, "TIME", time)
    filevarattdef(fout, "TOS", TOS)

    
    setfileoption(fout,"DefineMode",False)

    fout->LAT=(/TOS&LAT/)
    fout->LONN179_180=(/TOS&LONN179_180/)
    fout->TIME=(/TOS&TIME/)
    fout->TOS=(/TOS/)

end 