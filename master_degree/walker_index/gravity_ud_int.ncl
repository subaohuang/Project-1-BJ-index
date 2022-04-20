function gain_data(data1)
 begin
 
 rc1=regCoef(data1&time,data1)


 data_fuc=new((/3,dimsizes(data1&time)/),float)


 data_fuc(0,:) = data1(:)
 data_fuc(1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(2,:) =avg(data1)

return(data_fuc)
end 

begin
diri="/mnt/f/share/era/"
fils1=systemfunc("ls "+diri+"msf_time_era.nc")
f1=addfile(fils1, "r")
lon=f1->lon(:)
ud_time=f1->ud_time
printVarSummary(ud_time)
 print("level="+ud_time&level)
plev=ud_time&level
plev@units="hpa"
 psfc = 1008.0
 ptop = 100 
  dp   = dpres_plevel(plev, psfc, ptop, 0)
print("dp=" +dp)

m=new((/dimsizes(ud_time&level),dimsizes(ud_time&lon),dimsizes(ud_time&time)/),float)
do l=0,dimsizes(ud_time&lon)-1
    do t=0,dimsizes(ud_time&time)-1
 m(:,l,t)=abs(ud_time(:,l,t)*dp)
end do
end do
m!0="level"
m!1="lon"
m!2="time"
m&level=ud_time&level
m&lon=ud_time&lon
m&time=ud_time&time

M=dim_sum_n_Wrap(m, (/0,1/))
print(M)
print("m="+m(:,1,1))
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

 print("x="+x)
 print("y="+y)
x!0="time"
x&time=m&time

y!0="time"
y&time=m&time

data=gain_data(y)
data2=gain_data(x)
wks1=gsn_open_wks("png", ""+diri+"ceshi_gravity_y")
wks2=gsn_open_wks("png", ""+diri+"ceshi_gravity_x")
res=True
plot1=gsn_csm_xy(wks1, m&time, data, res)
plot2=gsn_csm_xy(wks2, m&time, data2, res)

end 
