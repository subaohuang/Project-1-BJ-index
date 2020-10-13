begin
 diri="/mnt/f/share/era/data/"
 fils1=systemfunc("ls "+diri+"era_sur_global.nc")
 f1=addfile(fils1, "r")

 fils2=systemfunc("ls "+diri+"era_global_u_v_1.nc")
 f2=addfile(fils2, "r")
 
 sst_r=f1->sst(0:455,:,:)
 u10_r=f2->u(0:455,{850},:,:)
 v10_r=f2->v(0:455,{850},:,:)
 sst=sst_r*sst_r@scale_factor+sst_r@add_offset-273.5
 u10=u10_r*u10_r@scale_factor+u10_r@add_offset
 v10=v10_r*v10_r@scale_factor+v10_r@add_offset
 copy_VarCoords(sst_r, sst)
 copy_VarCoords(u10_r, u10)
 copy_VarCoords(v10_r, v10)
 printVarSummary(sst)
 ssta=rmMonAnnCycTLL(sst)
 u10a=rmMonAnnCycTLL(u10)
 v10a=rmMonAnnCycTLL(v10)
 time=cd_calendar(f1->time(0:455), 1)
 
 delete(ssta&time)
 delete(u10a&time)
 delete(v10a&time)
 ssta&time=time
 u10a&time=time
 v10a&time=time

 
;;;1994年海表温度
 sst_son=dim_avg_n_Wrap(ssta(188:190,{30:-30},:),0)
;;;1997年enino
 sst_1997=dim_avg_n_Wrap(ssta(224:226,{30:-30},:),0)
;;;
 u10_1997=dim_avg_n_Wrap(u10a(224:226,{30:-30},:),0)
 v10_1997=dim_avg_n_Wrap(v10a(224:226,{30:-30},:),0)

 u10_1994=dim_avg_n_Wrap(u10a(188:190,{30:-30},:),0)
 v10_1994=dim_avg_n_Wrap(v10a(188:190,{30:-30},:),0)
;;;画1994到1995年的海温变化情况
  ssta_1994_5=ssta({199405:199502},:,:)
  u10a_1994_5=u10a({199405:199502},:,:) 
  v10a_1994_5=v10a({199405:199502},:,:)

  ssta_1992_3=ssta({199203:199302},:,:)
  u10a_1992_3=u10a({199203:199302},:,:) 
  v10a_1992_3=v10a({199203:199302},:,:)


;;;画1967年的海温情况
    
wks=gsn_open_wks("pdf", ""+diri+"sst_globle1994")
wks1=gsn_open_wks("pdf", ""+diri+"sst_globle1997")
wks2=gsn_open_wks("pdf", ""+diri+"sst_1994到199502季节")
wks3=gsn_open_wks("pdf", ""+diri+"sst_199203到199302季节")

gsn_define_colormap(wks, "BlueWhiteOrangeRed")
gsn_define_colormap(wks1, "BlueWhiteOrangeRed")
gsn_define_colormap(wks2, "BlRe")
gsn_define_colormap(wks3, "BlRe")
  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode ="ManualLevels"
  res@cnMinLevelValF = -1.5
  res@cnMaxLevelValF = 1.5
  res@cnLevelSpacingF = 0.25
  res@lbLabelBarOn=False
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(sst_son&latitude)
  res@mpMaxLatF             = max(sst_son&latitude)
  res@mpMinLonF             = min(sst_son&longitude)
  res@mpMaxLonF             = max(sst_son&longitude)
  res@mpCenterLonF          = avg(sst_son&longitude)
  res@gsnScalarContour        = True

  res@vcRefMagnitudeF          = 3.0           ; define vector ref mag
  res@vcRefLengthF             = 0.02        ; define length of vec ref
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  res@vcMinDistanceF           = 0.02         ; thin vectors
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角

plot=gsn_csm_vector_scalar_map(wks, u10_1994,v10_1994, sst_son, res)
plot1=gsn_csm_vector_scalar_map(wks1, u10_1997,v10_1997,sst_1997, res)

draw(plot1)
draw(plot)
frame(wks)
frame(wks1)

plot2=new(10,graphic)
plot3=new(12,graphic)
do i=0,9
  plot2(i)=gsn_csm_vector_scalar_map(wks2, u10a_1994_5(i,:,:),v10a_1994_5(i,:,:),ssta_1994_5(i,:,:), res)
end do

do i=0,11
  plot3(i)=gsn_csm_vector_scalar_map(wks3, u10a_1992_3(i,:,:),v10a_1992_3(i,:,:),ssta_1992_3(i,:,:), res)
end do

pres=True
pres@gsnPanelLabelBar=True
gsn_panel(wks2, plot2, (/5,2/), pres)

gsn_panel(wks3, plot3, (/6,2/), pres)


end
