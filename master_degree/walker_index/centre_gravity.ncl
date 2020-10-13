begin
 diri="/mnt/f/share/era/"
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 ; fils2=systemfunc("ls "+diri+"msf_era_level_location1979-2012.nc")
 ; f2=addfile(fils2, "r")
 ; locationxi=f2->locationxi
 ; locationdong=f2->locationdong
 ; printVarSummary(locationxi)
 level=f1->level(:)

 lon=f1->lon(:)
 lonind1=ind(lon.eq.140)
 lonind2=ind(lon.eq.280)
  msf_time=f1->msf_time
 ; ud=f1->ud(:,lonind1:lonind2)
 ; printVarSummary(ud)
 msf_time_scale=doubletofloat(msf_time(:,lonind1:lonind2,:)/10^9)
 copy_VarMeta(msf_time(:,lonind1:lonind2,:),msf_time_scale)
 delete(msf_time)
 printVarSummary(msf_time_scale)
 pout=fspan(100,1000,901)

 msf_int=int2p_n_Wrap(msf_time_scale&level, msf_time_scale,pout(::-1), 1, 0)
 delete(msf_time_scale)
 delete(pout)
 xout=fspan(140, 280, 141)
 msf_int_f=linint1_n_Wrap(msf_int&lon, msf_int, False,xout,0,1)
 delete(msf_int)
 delete(xout)
 msf_int_f=where(msf_int_f.ge.1.0, msf_int_f, 0.0)

 
  m=new((/450,70,dimsizes(msf_int_f&time)/),float)
 do i=1,dimsizes(msf_int_f&LEVEL)-2,2
    do j=1,dimsizes(msf_int_f&LON)-2,2
        do t=0,dimsizes(msf_int_f&time)-1
           m((i-1)/2,(j-1)/2,t)=abs((msf_int_f(i-1,j-1,t)+msf_int_f(i-1,j,t)+msf_int_f(i-1,j+1,t)-(msf_int_f(i+1,j-1,t)+msf_int_f(i+1,j,t)+msf_int_f(i+1,j+1,t)))/3.0)
     end do
 end do
end do
m!0="level"
m!1="lon"
m!2="time"
hlevel=ispan(101,999,2)
m&level=hlevel(::-1)
m&lon=fspan(141, 279,70)
m&time=msf_int_f&time
M=dim_sum_n_Wrap(m, (/0,1/))
delete(msf_int_f)
x=new((/dimsizes(m&time)/),float)
y=new((/dimsizes(m&time)/),float)
k=0.0
printVarSummary(k)
do t=0,dimsizes(m&time)-1
 do i=0,dimsizes(m&level)-1
  do j=0,dimsizes(m&lon)-1
     k=m(i,j,t)*m&level(i)+k
   end do
 end do   
 y(t)=k/M(t)
 k=0.0
 end do

do t=0,dimsizes(m&time)-1
 do i=0,dimsizes(m&level)-1
  do j=0,dimsizes(m&lon)-1
     k=m(i,j,t)*m&lon(j)+k
   end do
 end do   
 x(t)=k/M(t)
 k=0.0
 end do
x!0="time"
x&time=m&time

y!0="time"
y&time=m&time
ntime=dimsizes(x&time)
nlevel=dimsizes(m&level)
nlon=dimsizes(m&lon)
;---------------------------------------------------------------------
;写文件
;--------------------------------------------------------------------
    filo="msf_era_gravity_center_location1979-2017.nc"
    system("/bin/rm -f "+diri + filo+"")
    fout=addfile(diri + filo,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "重心位置"    
    fAtt@creation_date = "2017.8.18"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlevel,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
    ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level", "float",(/"level"/))
    filevardef(fout, "lon", "float",(/"lon"/))
    filevardef(fout, "time", "float",(/"time"/))
    filevardef(fout, "x", typeof(x),getvardims(x))
    filevardef(fout, "y", typeof(y),getvardims(y))
    filevardef(fout, "m", typeof(m),getvardims(m))
   

    filevarattdef(fout, "x", x)
    filevarattdef(fout, "y", y)
    filevarattdef(fout, "m", m)
   

    setfileoption(fout,"DefineMode",False)

    fout->x=(/x/)
    fout->y=(/y/)
    fout->m=(/m/)
    fout->level=(/m&level/)
    fout->lon=(/m&lon/)
    fout->time=(/m&time/)
end 



















