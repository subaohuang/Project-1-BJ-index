;计算3D Q 非绝热加热项Q的情况
function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin


     data_a=data
    ; data_a2=dim_standardize_n_Wrap(data_a,1,0)        ;标准化
    ; data_a1=doubletofloat(data_a2)

    ; data_a1&time=ispan(1,456,1)
    data_a&time=ispan(1,dimsizes(data&time),1)
    data_a@units=units
    data_a@long_name=long_name
    data_a_new = data_a(level|:,lon|:,time|:)

    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)

    ;-----------------------------------------------
    rc = regCoef(data_a_new&time,data_a_new)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data_a_new&level
    rc&lon=data_a_new&lon
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
return(rc)
end


begin
 diri="/media/sf_share/ncep/ncep1/"
 fils1=systemfunc("ls "+diri+"omega.mon.mean.nc")
 fils2=systemfunc("ls "+diri+"air.mon.mean.nc")
 f1=addfile(fils1, "r")
 f2=addfile(fils2, "r")
 time=f1->time
 lon=f1->lon
 lat=f1->lat
 level=f1->level
 date=cd_calendar(time, 1)

 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201701)
 latind1=ind(lat.eq.5)
 latind2=ind(lat.eq.-5)
 lonind1=ind(lon.eq.0)
 lonind2=ind(lon.eq.357.5)
 levind1=ind(level.eq.100)
 levind2=ind(level.eq.1000)
 
 omega=f1->omega(tind1:tind2,levind2:levind1,latind1:latind2,lonind1:lonind2)
 w=omega
 copy_VarMeta(omega, w)
 delete(omega)

 tem=f2->air(tind1:tind2,levind2:levind1,latind1:latind2,lonind1:lonind2)
 printMinMax(tem, 1)
 t=tem+273.15
 copy_VarMeta(tem, t)
 t@units="degK"
 printVarSummary(t)
 delete(tem)
 p=level(levind2:levind1)
 
 pt=pot_temp(p, t, 1, False)
 printMinMax(pt,1)
                                                          
 dtdp=center_finite_diff_n(pt, p*100,False , 0, 1)                                                         ;计算稳定度
 copy_VarCoords(pt, dtdp)


 
 Q=new((/dimsizes(pt&time),dimsizes(pt&level),dimsizes(pt&lat),dimsizes(pt&lon)/),float)
 do i=0,dimsizes(p)-1
 Q(:,i,:,:)=((p(i)/1000.0)^0.286)*w(:,i,:,:)*dtdp(:,i,:,:)*24.0*3600
 end do
 delete(i)
 copy_VarCoords(pt, Q)
 Q_new=dim_avg_n_Wrap(Q(time|:,lat|:,level|:,lon|:),(/0,1/))
 Q_time=dim_avg_n_Wrap(Q(time|:,lat|:,level|:,lon|:), 1)
 rc_Q=decade_trend(Q_time,"k/day","dabatic heating")
 printVarSummary(rc_Q)
   tval = onedtond(rc_Q@tval , dimsizes(rc_Q))
   df   = onedtond(rc_Q@nptxy, dimsizes(rc_Q)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
    copy_VarCoords(rc_Q,prob )
    prob1d=ndtooned(prob)
    i1d=ind(prob1d.lt.0.05)
    i=ind_resolve(i1d, dimsizes(prob))

    ilev=prob&lev(i(:,0))
    ilon=prob&lon(i(:,1))
    
 
 wks=gsn_open_wks("png", ""+diri+"Q_ncep1")
  gsn_define_colormap(wks,"BlRe")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="necp1"
    res@gsnRightString="1979-2017"
    ; res@lbLabelStrings=ispan(-6,6,1)
  
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.004
    res@cnMaxLevelValF = 0.004
    res@cnLevelSpacingF = 0.001
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    ; res@cnMonoFillColor=False
    
    res@vpHeightF=0.5
    res@vpWidthF=0.7
    res@vpXF=0.15
    res@vpYF=0.9
    res@cnFillColors=(/8,11,14,18,48,50,67,74,84,93/)
    res@cnLineLabelsOn=False
    ; res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res@cnLabelScalingMode="ScaleFactor"
       ; res@cnLabelScaleValueF=10^9
      
    res@lbOrientation="Horizontal"
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
     res@tmXBValues=(/0,90,180,270,357.5/)
     res@tmXBLabels=(/"0","90~S~o~N~E","180~S~o~N","90~S~o~N~W","0"/)
     res@tmXBMinorOn=True
 
   plot=gsn_csm_pres_hgt(wks,rc_Q,res)
    getvalues plot@contour                                                     
  "pmAnnoManagers" : am_ids1       
  end getvalues                  
           
  index1 = ind(NhlName(am_ids1).eq."right_axis")  
  if(.not.ismissing(index1)) then                                             
    NhlRemoveAnnotation(plot@contour,am_ids1(index1))                          
  end if 


     polyres                   = True
   polyres@gsMarkerIndex     = 16          ; polymarker style
   polyres@gsMarkerSizeF     = 3.         ; polymarker size

    dum1=gsn_add_polymarker(wks, plot, ilon, ilev, polyres)

 res@cnLevelSelectionMode = "ExplicitLevels"
 res@cnLinesOn=True
 res@cnLevels = fspan(-1.6,2.4,6)
 res@cnLineLabelsOn=True
 res@cnLineLabelBackgroundColor="Transparent"
 res@cnLineThicknessF = 3.
 res@cnFillOn        = False
 res@gsnContourNegLineDashPattern=1
 plot1 = gsn_csm_pres_hgt(wks, Q_new, res)
;-----------------------------------------------------
 ;去掉height字体
;----------------------------------------------------
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



 end