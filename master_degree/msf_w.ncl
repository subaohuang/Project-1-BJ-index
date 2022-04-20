
function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin

    ; data_new=doubletofloat(dim_standardize_n_Wrap(data(::-1,:,:),1,(/0,1/)))       ;标准化
    ; print(data_new)
    ; copy_VarMeta(data(::-1,:,:), data_new)
    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)
    ;-----------------------------------------------
    rc = regCoef_n(data&time,data,0,2)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data&level
    rc&lon=data&lon
    rc&lev@units = "hpa"
    rc&lon@units = "degrees_east"
     ; print(rc)
return(rc)
end
;;;-------
begin
 diri="/mnt/f/share/era/data/"
 fils1=systemfunc("ls "+diri+"msf_time_era_global.nc")
 f1=addfile(fils1, "r")

 fils2=systemfunc("ls "+diri+"msf_era_location1979-2017.nc")
 f2=addfile(fils2, "r")
 location=f2->locationxi_a

 msf=f1->msf(:,:)
 ud=f1->ud(:,:)

 w=f1->w(:,:)
 level=f1->level(:)
 lon=f1->lon(:)
 msf_time=f1->msf_time
 msf_time=(/msf_time/10^9/)

 ud_time=f1->ud_time
 w_time=f1->w_time
 msf&level=level
 msf&lon=lon
 msf&level@units="hpa"
 copy_VarCoords(msf, ud)
 copy_VarCoords(msf, w)
 copy_VarCoords(msf_time, ud_time)
 copy_VarCoords(msf_time, w_time)
 w_time=(/w_time*2/)
 w_clm=dim_avg_n_Wrap(w_time, 2)
 msf_clm=dim_avg_n_Wrap(msf_time, 2)


 w_rc=decade_trend(w_time,"units","w_rc")
 w_rc=(/120*w_rc/)
  
 location_clm=dim_avg_n_Wrap(location, 0)

  ;;;
   tval = onedtond(w_rc@tval , dimsizes(w_rc))
   df   = onedtond(w_rc@nptxy, dimsizes(w_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   
    copy_VarCoords(w_rc,prob )
    printVarSummary(prob)
    prob1d=ndtooned(prob)
    i1d=ind(prob1d.lt.0.05)
    i=ind_resolve(i1d, dimsizes(prob))

    ilev=prob&lev(i(:,0))
    ilon=prob&lon(i(:,1))
   delete(i)
  ;;;
  
;;;;---------------    
  wks=gsn_open_wks("pdf", ""+diri+"msf_w")

  wks2=gsn_open_wks("pdf", ""+diri+"w_trend&w")
  ;;;画图
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="Climatology"
    
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnFillPalette  = "GMT_polar"
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -5
    res@cnMaxLevelValF = 5
    res@cnLevelSpacingF = 1
  
    res@cnFillOn = True 

    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res@cnLabelScalingMode="ScaleFactor"
       ; res@cnLabelScaleValueF=10^9
      
    res@lbOrientation="Horizontal"
 
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
    res@tmYRMode="Automatic"
    

    
    res@vpHeightF=0.3
    res@vpWidthF=0.7
    
    plot=gsn_csm_pres_hgt(wks,w_clm,res)
 
    
   
  ;;---Resources for the polyline  添加下垫面的陆地
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 15.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=15
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 25)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)

    ployres@gsLineDashPattern=0
    
    
  dum1= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)
  ployres@gsLineDashPattern=16
  dum2= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)
  ployres@gsLineDashPattern=0
  dum3= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)
  
  ploylon=conform_dims(dimsizes(level), location_clm, 0)
  ployres@gsLineColor      = "green"         ; color of lines
  ployres@gsLineThicknessF = 4.0            ; line thickness   
  dum7=gsn_add_polyline(wks, plot ,ploylon ,level, ployres)

   ployres@gsLineThicknessF = 15.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
;;;添加msf的等值线

 res@cnLevelSelectionMode = "ExplicitLevels"
 res@cnLinesOn=True
 res@cnLevels = ispan(-6,6,1)
 res@cnLineLabelsOn=True
 res@cnLineLabelBackgroundColor="Transparent"
 res@cnLineThicknessF = 3.
 res@cnFillOn        = False
 res@gsnContourNegLineDashPattern=1

  plot1 = gsn_csm_pres_hgt(wks, msf_clm, res)

  
    

overlay(plot,plot1)
draw(plot)
frame(wks)
delete(plot)
delete(plot1)

;----------------------------------------------------------------------
;画随时间变化
;----------------------------------------------------------------------
   gsn_define_colormap(wks2,"GMT_polar")
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="era"
    res2@gsnRightString="1979-2017"
    ; res2@lbLabelStrings=ispan(-6,6,1)
  
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -0.5
    res2@cnMaxLevelValF = 0.5
    res2@cnLevelSpacingF = 0.1
    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    ; res2@cnMonoFillColor=False
    
    ; res2@cnFillColors=(/2,10,15,21,30,40,49,50,60,66,74,81,90,97/)
    res2@cnLineLabelsOn=False
    ; res2@cnMonoFillScale=True
    res2@cnInfoLabelOn=False

     ; res2@cnExplicitLineLabelsOn=True
     ; res2@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res2@cnLabelScalingMode="ScaleFactor"
       ; res2@cnLabelScaleValueF=10^9
      
    res2@lbOrientation="Horizontal"
    
    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    res2@tmXBValues=(/0,90,180,270,357.5/)
    res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(0,360,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7

    plot=gsn_csm_pres_hgt(wks2,w_rc,res2)
  ; ;;;通过显著性的点
  ;    polyres                   = True
  ;    polyres@gsMarkerIndex     = 16          ; polymarker style
  ;    polyres@gsMarkerSizeF     = 3.         ; polymarker size

  ;   dum1=gsn_add_polymarker(wks2, plot, ilon, ilev, polyres)

;;;gsn_contour_shade
   res_pro=True
   res_pro@gsnDraw=False
   res_pro@gsnFrame=False
   res_pro@cnLevelSelectionMode = "ExplicitLevels"
   res_pro@cnLinesOn=False
   res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing
 ; res_pro@cnLevels = ispan(-6,6,1)
   res_pro@cnLineLabelsOn=False
   res_pro@tmYRMode = "Automatic"
   res_pro@cnInfoLabelOn=False
 ; res_pro@cnLineLabelBackgroundColor="Transparent"
 ; res_pro@cnLineThicknessF = 3.
   ; res_pro@cnFillOn        = True
  
   iplot=gsn_csm_pres_hgt(wks2,prob,res_pro)
  

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
  iplot = gsn_contour_shade(iplot, 0,0.05, opt)
 ;;;

 res2@cnLevelSelectionMode = "ExplicitLevels"
 res2@cnLinesOn=True
 res2@cnLevels = ispan(-5,5,1)
 res2@cnLineLabelsOn=True
 res2@cnLineLabelBackgroundColor="Transparent"
 res2@cnLineThicknessF = 2.
 res2@cnFillOn        = False
 res2@gsnContourNegLineDashPattern=1

 plot1 = gsn_csm_pres_hgt(wks2, w_clm, res2)
;;;添加地形
 
  ployres@gsLineDashPattern=0
  dum4= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)
  ployres@gsLineDashPattern=16
  dum5= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)
  ployres@gsLineDashPattern=0
  dum6= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)

  ployres@gsLineColor      = "green"         ; color of lines
  ployres@gsLineThicknessF = 4.0            ; line thickness   
  dum8=gsn_add_polyline(wks, plot ,ploylon ,level, ployres)


 overlay(plot, iplot)
 overlay(plot,plot1)
 draw(plot)
 frame(wks2)


end
   