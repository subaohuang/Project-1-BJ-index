
begin
 diri="/mnt/f/share/era/data/data_series/"
 fils1=systemfunc("ls "+diri+"era_pre_global.nc")
 f1=addfile(fils1, "r")
 time=f1->time
 lon=f1->longitude
 lat=f1->latitude
 level=f1->level
 date=cd_calendar(time, 1)
 timeind=ind(date.ge.197901.and.date.le.201212)
 levind=ind(level.ge.100.and.level.le.1000)

 w_r=f1->w(timeind,levind,{5:-5},:)
 

if isatt(w_r, "add_offset") then 
w=w_r@scale_factor*w_r + w_r@add_offset
print("**************"+"jius"+"***")
else 
w=w_r
end if
copy_VarCoords(w_r, w)
copy_VarAtts(w_r, w)
w!0="time"
w!1="level"
w!1="lat"
w!2="lon"
printVarSummary(w)
delete(w_r)
;;;读入T
t_r=f1->t(timeind,levind,{-5:5},:)

if t_r@units.eq."degC" then 
   t_r=(/t_r+273.16/)
end if

if isatt(t_r, "add_offset") then 
t=t_r@scale_factor*t_r + t_r@add_offset
print("**************"+"jius"+"***")
else 
t=t_r
end if
copy_VarCoords(t_r, t)
copy_VarAtts(t_r, t)
t!0="time"
t!1="level"
t!2="lat"
t!3="lon"
printVarSummary(t)
delete(t_r)

 p=level(levind)*100.0

 pt=pot_temp(p, t, 1, False)
 printVarSummary(pt)
 delete(t)
 printMinMax(pt, 1)                                        
 dtdp=center_finite_diff_n(pt, p,False , 0, 1)                                                         ;计算稳定度
 copy_VarCoords(pt, dtdp)
 printMinMax(w, 1)
 printVarSummary(dtdp)
 Q=new((/dimsizes(pt&time),dimsizes(pt&level),dimsizes(pt&lat),dimsizes(pt&lon)/),float)
 do i=0,dimsizes(p)-1
 Q(:,i,:,:)=tofloat(((p(i)/100000.0)^0.286)*w(:,i,:,:)*dtdp(:,i,:,:)*3600.0*24.0)
 end do
 delete(w)
 copy_VarCoords(pt, Q)
 Q_time=dim_avg_n_Wrap(Q(time|:,lat|:,level|:,lon|:), 1)
 Q_time_new=Q_time(level|:,lon|:,time|:)
nlev=dimsizes(Q_time_new&level)
nlon=dimsizes(Q_time_new&lon)
ntime=dimsizes(Q_time_new&time)

;****************************************************************************
;随时间变化的量
    diro=diri
    filo1="Q_time_era.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "Q的值，取5到-5平均了"    
    fAtt@creation_date = "2017.8.25"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(level), (/"level"/))
    filevardef(fout, "lon", typeof(lon),(/"lon"/))
    filevardef(fout, "time", typeof(time),(/"time"/))
    filevardef(fout, "Q_time", typeof(Q_time_new),(/"level","lon","time"/))    
 
    
  
    filevarattdef(fout, "level", Q_time_new&level)
    filevarattdef(fout, "time", Q_time_new&time)
    filevarattdef(fout, "lon", Q_time_new&lon)

    setfileoption(fout,"DefineMode",False)

    fout->level=(/Q_time_new&level/)
    fout->lon=(/Q_time_new&lon/)
    fout->time=(/Q_time_new&time/)
    fout->Q_time=(/Q_time_new/)

 end