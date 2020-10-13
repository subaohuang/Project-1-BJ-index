begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_global_u_v_1.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)
u1=f1->u(tind,{850},:,:)
v1=f1->v(tind,{850},:,:)

u=u1*u1@scale_factor+u1@add_offset
v=v1*v1@scale_factor+v1@add_offset
copy_VarCoords(u1, u)
copy_VarCoords(v1, v)
delete(u1)
delete(v1)

result = uv2sfvpF(u(:,::-1,:),v(:,::-1,:))


printVarSummary(result)
pv=result(1,:,::-1,:)
delete(result)
copy_VarCoords(u(:,::-1,:), pv)
delete(u)
delete(v)
printVarSummary(pv)
; ntime=dimsizes(data)
;;输出至Netcdf文件
;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。
   diro=diri
    filo1="velocity_potensial.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "速度势"    
    fAtt@creation_date = "2017.9.26"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"time","lat","lon"/)
    ; dimsize=(/ntime,nlat,nlon/)
    dimsize=dimsizes(pv)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "lat",typeof(pv&latitude), (/"lat"/))
    filevardef(fout, "lon", typeof(pv&longitude),(/"lon"/))
    filevardef(fout, "time", typeof(pv&time),(/"time"/))
    filevardef(fout, "pv", typeof(pv),(/"time","lat","lon"/))

    filevarattdef(fout, "lat", pv&latitude)
    filevarattdef(fout, "lon", pv&longitude)


    
    setfileoption(fout,"DefineMode",False)

    fout->lat=(/pv&latitude/)
    fout->lon=(/pv&longitude/)
    fout->time=(/pv&time/)
   
    fout->pv=(/pv/)
   
end