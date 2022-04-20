;***************11
undef("overlay_x_longitude")
function overlay_x_longitude(wks,plot,latlon,color)
begin

f=addfile("/data6/zhanggsh/data/topography/orography_etopo_to_surface_pressure.nc","r")
pres=f->ps
lat=f->lat
lon=f->lon

 ydim = dimsizes(pres)
;x-axis is longitude
 yp  = new((/ydim(0),2*ydim(1)/),"double")
 xp  = new((/2*ydim(1)/),"double")
 do k = 0,ydim(1)-1
 yp(:,k) = pres(:,k)
 yp(:,2*ydim(1)-1-k) = 1000
 xp(k) =  lon(k)
 xp(2*ydim(1)-1-k) = lon(k)
 end do

  gsres  = True
  gsres@gsFillColor= color
  gsres@gsLineColor =color

 ;yp_ave=dim_avg_n_Wrap(yp({28:38},:), 0)
 
;dum =gsn_add_polygon(wks,plot,xp,yp_ave,gsres) 
dum =gsn_add_polygon(wks,plot,xp,yp({latlon},:),gsres)
return(dum)
end




undef("overlay_x_latitude")
function overlay_x_latitude(wks,plot,latlon,color)
begin

f=addfile("/data6/zhanggsh/data/topography/orography_etopo_to_surface_pressure.nc","r")
pres=f->ps
lat=f->lat
lon=f->lon
;x-axis is latitude
 ydim = dimsizes(pres)
 yp1  = new((/2*ydim(0),ydim(1)/),"double")
 xp1  = new((/2*ydim(0)/),"double")
 do k = 0,ydim(0)-1
 yp1(k,:) = pres(k,:)
 yp1(2*ydim(0)-1-k,:) = 1000
 xp1(k) =  lat(k)
 xp1(2*ydim(0)-1-k) = lat(k)
 end do
 delete(ydim)

  gsres  = True
  gsres@gsFillColor= color
  gsres@gsLineColor =color

 yp1_ave=dim_avg_n_Wrap(yp1(:,{90:105}), 1)
dum =gsn_add_polygon(wks,plot,xp1,yp1_ave,gsres)

;dum =gsn_add_polygon(wks,plot,xp1,yp1(:,{latlon}),gsres)
return(dum)
end

undef("overlay_tp_horizontal")
function overlay_tp_horizontal(wks,lat1,lat2,lon1,lon2,color)
begin
  
    f1 = addfile("/home/ys17-19/lsh/Project/other/guanshun/dixing_0.25X0.25.nc", "r")
    ; f1= addfile("/home/ys17-19/lsh/data/ETOPO1_Bed_g_gmt4.grd.nc", "r")
    z = short2flt(f1->z)
    z = (z/9.8)
    z_clm  = dim_avg_n_Wrap(z, 0)
    ; z_clm = z



    tbres=True
    tbres@gsnFrame=False
    tbres@gsnDraw=False
    tbres@tmYROn=False          ;关闭右Y轴刻度
    tbres@tmXTOn=False          ;关闭右X轴刻度
    tbres@cnLineLabelsOn =False
    tbres@cnFillOn =True
    tbres@cnLinesOn=False
    tbres@gsnLeftString = ""
    tbres@gsnRightString=" "
    tbres@cnLineColor="orangered4"
    tbres@lbLabelBarOn     =False
    tbres@cnLevelSelectionMode = "ExplicitLevels"
    tbres@cnLevels = (/3000/)
    tbres@cnFillColors=(/"Transparent",color/)
    tbres@tiMainString=""
    tbres@cnInfoLabelOn =False
    ; tbres@cnFillDrawOrder="PostDraw"
    plot1= gsn_csm_contour(wks,z_clm({lat1:lat2},{lon1:lon2}),tbres) 
return plot1
end 