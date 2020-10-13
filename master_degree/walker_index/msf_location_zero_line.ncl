function gain_data(data1,data2)
 begin
 
 rc1=regCoef(data1&time,data1)
 rc2=regCoef(data2&time,data2)
print(rc1)
print(rc2)
 data_fuc=new((/2,3,dimsizes(data1&time)/),double)

 data_fuc(0,0,:) = data1(:)
 data_fuc(0,1,:) = rc1*data1&time+rc1@yintercept
 data_fuc(0,2,:) =avg(data1(:))

 data_fuc(1,0,:) = data2(:)
 data_fuc(1,1,:) = rc2*data2&time+rc2@yintercept
 data_fuc(1,2,:) =avg(data2(:))
 
return(data_fuc)
end 

function draw_panel(wks1,wks2,time,data1,data2,title)
 begin
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@tmXBMode="Explicit"
    ; resxy@trXMaxF=456
    resxy@tmXBValues=ispan(100,240,20)
    resxy@tmXBLabels=ispan(100,240,20)+"~S~o~N~E"
    resxy@tmXBMinorOn=False

    
    resxy@tmYLMode="Explicit"
    
    resxy@trYMaxF=dimsizes(data1&time)
    resxy@tmYLValues=ispan(13,dimsizes(data1&time),60)
    resxy@tmYLLabels=ispan(1980,2015,5)

    resxy@vpHeightF=0.7
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.9

    resxy@xyLineThicknesses = (/ 2.0,1.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
    resxy@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/)   

    resxy@xyDashPatterns    =(/0,0,2/)
    

    resxy@tmXBLabelFontAspectF=2.5
    resxy@tmYLLabelFontAspectF=2.5
    resxy@tmXTOn=False
    resxy@tmYROn=False
    resxy@tiYAxisString = title

    ; resxy@tmLabelAutoStride = True                    ; nice stride on labels
    contour=new(2,graphic)
    resxy2=resxy
    delete(resxy2@tmXBValues)
    delete(resxy2@tmXBLabels)

    resxy2@tmXBValues=ispan(180,280,20)
    resxy2@tmXBLabels=""+ispan(180,280,20)+"~S~o~N~E"

    contour(0)=gsn_csm_xy(wks1,data1,time,resxy)
    contour(1)=gsn_csm_xy(wks2,data2,time,resxy2)


    ; contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    ; contour(1)=gsn_csm_xy(wks, time,data2,resxy)

return contour
end 

begin
 diri="/media/sf_share/ncep/ncep1/"
 fils1=systemfunc("ls "+diri+"msf_necp1.nc")
 fils2=systemfunc("ls "+diri+"msf_time_necp1.nc")
 fils4=systemfunc("ls "+diri+"omega.mon.mean.nc")
 f1=addfile(fils1, "r")
 f2=addfile(fils2, "r") 
 f4=addfile(fils4, "r")
  timew=f4->time
  latw=f4->lat
  lonw=f4->lon
  date=cd_calendar(timew, 1)
  tind1=ind(date.eq.197901)
  tind2=ind(date.eq.201701)
  latind1=ind(latw.eq.5)
  latind2=ind(latw.eq.-5)
  lonind1=ind(lonw.eq.100)
  lonind2=ind(lonw.eq.290) 
  w=f4->omega(tind1:tind2,:,latind2:latind1,lonind1:lonind2)

  w_new=w(time|:,lat|:,level|:,lon|:) 
  w_avg=dim_avg_n_Wrap(w_new, 1)
  w_avg_new=w_avg*(-50)
  copy_VarMeta(w_avg, w_avg_new)

 msf_time=f2->msf_time
 ud_time=f2->ud_time
 msf=f1->msf
 level=f1->level(:)
 lon=f1->lon(:)

 msf_time&level=level
 msf_time&lon=lon
 msf_time&level@units="hpa"

 ud_time&level=level
 ud_time&lon=lon
 ud_time&level@units="hpa"

 msf_scale=msf_time/10^9
 copy_VarMeta(msf_time, msf_scale)
 locationxi=new((/dimsizes(msf_time&time)/),float)
 locationdong=new((/dimsizes(msf_time&time)/),float)
 wks1=gsn_open_wks("png", ""+diri+"png/"+"locationxi")
 wks2=gsn_open_wks("png", ""+diri+"png/"+"locationdong")
 wks=gsn_open_wks("png", ""+diri+"png/"+"location")

 do i=0,dimsizes(msf_time&time)-1
 
 

 res = True
   res@gsnDraw=False
   res@gsnFrame=False
   res@gsnLeftString=""
   res@gsnRightString=""
   res@lbLabelStrings=ispan(-6,6,1)
   res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
    res@cnLinePalette="BlueDarkRed18"
    res@cnLineLabelsOn=True
    res@cnMonoFillScale=True
     res@cnInfoLabelOn=False

     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res@cnLabelScalingMode="ScaleFactor"
       ; res@cnLabelScaleValueF=10^9
      
     res@lbOrientation="Vertical"
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
     res@tiYAxisString=""
     res@tmXTOn=False
     res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(100,290,40)
    res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    res@tmXBMinorOn=False
    res@vcRefMagnitudeF = 2.0                ; define vector ref mag
    res@vcRefLengthF    = 0.04        ; define length of vec ref
    res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.02       ; thin out vectors
    res@vcMapDirection  = False

  
  ; res@vcFillArrowHeadXF=0.1
  res@vcFillArrowMinFracWidthF=1.0
  res@vcFillArrowHeadMinFracYF=1.0
  ; res@vcFillArrowHeadXF=1.0
  res@vcFillArrowWidthF=0.05
  res@vcFillArrowHeadYF=0.15
  
  res@vcRefAnnoOrthogonalPosF  = -1.05
  res@vcRefAnnoPerimOn=False


  plot=gsn_csm_pres_hgt_vector(wks,msf_scale(:,:,i),ud_time(:,:,i),w_avg_new(i,:,:),res)



 getvalues plot@contour                                                     
  "pmAnnoManagers" : am_ids1
                                         
  end getvalues                  
                          
  index1 = ind(NhlName(am_ids1).eq."right_axis")  
  if(.not.ismissing(index1)) then                                             
    NhlRemoveAnnotation(plot@contour,am_ids1(index1))                          
  end if 

;-------------------------------
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn=True
res@cnLevels = 0
res@cnLineThicknessF = 4.
res@cnFillOn        = False


plot1 = gsn_csm_pres_hgt(wks, msf_scale(:,:,i), res)

  getvalues plot1@contour                                                     
  "pmAnnoManagers" : am_ids   
                                   
  end getvalues                                                              

  index = ind(NhlName(am_ids).eq."right_axis")  
  if(.not.ismissing(index)) then                                             
    NhlRemoveAnnotation(plot1@contour,am_ids(index))                          
  end if 

overlay(plot,plot1)
draw(plot)
frame(wks)


;------------------------------
isolines=get_isolines(plot@contour, 0)
y=isolines(0,:)
x=isolines(1,:)

y=where(y.gt.300.and.y.lt.600, y, y@_FillValue)
x_xi=where(ismissing(y).ne.True.and.x.lt.220, x,x@_FillValue)
x_dong=where(ismissing(y).ne.True.and.x.lt.280.and.x.gt.180, x,x@_FillValue)

locationxi(i)=dim_avg(x_xi)
locationdong(i)=dim_avg(x_dong)
delete(res)
delete(plot)
delete(isolines)
delete(y)
delete(x)
delete(x_xi)
delete(x_dong)
end do






locationxi!0="time"
locationdong!0="time"
locationxi&time=msf_time&time
locationdong&time=msf_time&time
printVarSummary(locationxi)
data= gain_data(locationxi,locationdong)
    pres1=True
    pres2=True
    pres1@gsnMaximize=True
    contour=draw_panel(wks1,wks2,locationxi&time,data(0,:,:),data(1,:,:),"" )

    draw(contour(0))
    draw(contour(1))
    frame(wks1)
    frame(wks2)



end 
; print(iso@segment_count)
; print(iso(1,:))

   ;   print(fi_mean(0,:))