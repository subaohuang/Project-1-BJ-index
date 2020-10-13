begin
diri="/mnt/f/share/MRE/"
fils=systemfunc("ls "+diri+"mre_u_v.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201212.and.date.ge.195801)
u=f1->uwnd(tind,:,:,:)
v=f1->vwnd(tind,:,:,:)

result = uv2sfvpF(u(:,:,::-1,:),v(:,:,::-1,:))
printVarSummary(result)
pv=result(1,:,:,::-1,:)
delete(result)
copy_VarCoords(u, pv)
delete(u)
delete(v)
printVarSummary(pv)
; ntime=dimsizes(data)
;;输出至Netcdf文件
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
   diro=diri
    filo1="velocity_potensial_mre.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "速度势"    
    fAtt@creation_date = "2017.9.26"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","level","lat","lon"/)
    ; dimsize=(/ntime,nlat,nlon/)
    dimsize=dimsizes(pv)
    dimunlim=(/False,False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(pv&level), (/"level"/))
    filevardef(fout, "lat",typeof(pv&lat), (/"lat"/))
    filevardef(fout, "lon", typeof(pv&lon),(/"lon"/))
    filevardef(fout, "time", typeof(pv&time),(/"time"/))
    filevardef(fout, "pv", typeof(pv),(/"time","level","lat","lon"/))

    filevarattdef(fout, "lat", pv&lat)
    filevarattdef(fout, "lon", pv&lon)
    filevarattdef(fout, "level", pv&level)
    filevarattdef(fout, "time", pv&time)
    
    setfileoption(fout,"DefineMode",False)
    
    fout->level=(/pv&level/)
    fout->lat=(/pv&lat/)
    fout->lon=(/pv&lon/)
    fout->time=(/pv&time/)

    fout->pv=(/pv/)
   
end