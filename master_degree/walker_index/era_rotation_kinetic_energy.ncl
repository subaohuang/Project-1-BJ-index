function gain_data(data1)
 begin
 
 rc1=regCoef(data1&time,data1)


 data_fuc=new((/3,dimsizes(data1&time)/),double)


 data_fuc(0,:) = data1(:)
 data_fuc(1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(2,:) =avg(data1)

return(data_fuc)
end 
function gain_data_an(data1)
 begin
 
 intesion_an=runave_n_Wrap(dim_rmvmean(month_to_annual(data1(0:455), 1)),3,0,0)
 intesion_an!0="time"
 intesion_an&time=ispan(1,38,1)


 rc1=regCoef(intesion_an&time,intesion_an)

 data_fuc=new((/3,dimsizes(intesion_an&time)/),double)


 data_fuc(0,:) = intesion_an(:)
 data_fuc(1,:) = rc1*intesion_an&time+rc1@yintercept
 data_fuc(2,:) =avg(intesion_an)

return(data_fuc)
end 

begin
diri="/mnt/f/share/era/"
 fils1=systemfunc("ls "+diri+"msf_era_gravity_center_location1979-2017.nc")
 f1=addfile(fils1, "r")
 m=f1->m
 printVarSummary(m)
 fils2=systemfunc("ls "+diri+"era_geopoten.nc")
 f2=addfile(fils2, "r")

fils3=systemfunc("ls "+diri+"msf_time_era.nc")
f3=addfile(fils3, "r")
ud_time=f3->ud_time

  pout=ispan(101,999,2)
  ud_int = int2p_n_Wrap(ud_time&level,ud_time,pout(::-1),2,0)
  printVarSummary(ud_int)
  xo = ispan(141,279,2)
  ud_int_new=linint1_n_Wrap(ud_int&lon, ud_int, False, xo, 0, 1)
  printVarSummary(ud_int_new) 
  

ro_ki=new((/dimsizes(ud_int_new&time)/),float)
 k=0.0
  do t=0,dimsizes(m&time)-1
     do i=0,dimsizes(m&level)-1
      do j=0,dimsizes(m&lon)-1
      k=0.5*(ud_int(i,j,t))^2*m(i,j,t)+k
  end do 
 end do 
 ro_ki(t)=k
 k=0.0
end do
ro_ki!0="time"
ro_ki&time = ud_int_new&time
print(ro_ki)
data=gain_data(ro_ki)
data2=gain_data_an(ro_ki)

res=True
wks=gsn_open_wks("png", ""+diri+"ro_ki")
wks1=gsn_open_wks("png", ""+diri+"ro_ki_an")
plot=gsn_csm_xy(wks, ro_ki&time, data, res)
plot1=gsn_csm_xy(wks1, ispan(1,38,1), data2, res)
 end