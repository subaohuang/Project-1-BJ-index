
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

begin
 diri="/mnt/f/share/era/data/"
 fils1=systemfunc("ls "+diri+"msf_time_era.nc")
 f1=addfile(fils1, "r")
 msf=f1->msf(:,:)
 ud=f1->ud(:,:)

 w=f1->w(:,:)
 level=f1->level(:)
 lon=f1->lon(:)
 msf_time=f1->msf_time
 msf&level=level
 msf&lon=lon
 msf&level@units="hpa"
 copy_VarCoords(msf, ud)
 copy_VarCoords(msf, w)
 

 msf_scale=msf/10^9
 copy_VarMeta(msf, msf_scale)

 ; msf_time_scale=msf_time/10^9
 ; copy_VarMeta(msf_time, msf_time_scale)

 printVarSummary(msf_scale)
 msf_rc=decade_trend(msf_time,"units","msf_rc")

  msf_rc_scale=120*msf_rc/(10^9)
 copy_VarMeta(msf_rc, msf_rc_scale)
    
 printVarSummary(msf_rc_scale)
   
   tval = onedtond(msf_rc_scale@tval , dimsizes(msf_rc_scale))
   df   = onedtond(msf_rc_scale@nptxy, dimsizes(msf_rc_scale)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   
    copy_VarCoords(msf_rc_scale,prob )
    printVarSummary(prob)
    prob1d=ndtooned(prob)
    i1d=ind(prob1d.lt.0.05)
    i=ind_resolve(i1d, dimsizes(prob))

    ilev=prob&lev(i(:,0))
    ilon=prob&lon(i(:,1))
  
    
    
  wks=gsn_open_wks("png", ""+diri+"mass_stream_era")
  wks2=gsn_open_wks("png", ""+diri+"msf_rc_era")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="1979-2017"
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


   plot=gsn_csm_pres_hgt_vector(wks,msf_scale,ud,w,res)


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

plot1 = gsn_csm_pres_hgt(wks, msf_scale, res)
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
delete(plot)
delete(plot1)
;----------------------------------------------------------------------
;画随时间变化
;----------------------------------------------------------------------
    gsn_define_colormap(wks2,"BlRe")
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="era"
    res2@gsnRightString="1979-2017"
    ; res2@lbLabelStrings=ispan(-6,6,1)
  
    res2@cnLinesOn=False
    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -0.6
    res2@cnMaxLevelValF = 0.6
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
      
    res2@lbOrientation="Vertical"
    
    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
    res2@tmXBValues=ispan(100,290,40)
    res2@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    res2@tmXBMinorOn=False
    res2@tmYRMode  = "Automatic"
   plot=gsn_csm_pres_hgt(wks2,msf_rc_scale,res2)

  printVarSummary(msf_rc_scale)
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

   getvalues wks2
    "wkColorMapLen" : len
  end getvalues

  print(len)
  delete(plot)


  plot=gsn_csm_pres_hgt(wks2,msf_rc_scale,res2)

   getvalues plot@contour                                                     
  "cnFillColors" :cnFillColors
  end getvalues      
  print(cnFillColors)


  NhlSetColor(wks2,46,1,1,1)
  NhlSetColor(wks2,53,1,1,1)
  ;  polyres                   = True
  ;  polyres@gsMarkerIndex     = 16          ; polymarker style
  ;  polyres@gsMarkerSizeF     = 3.         ; polymarker size

  ;   dum1=gsn_add_polymarker(wks2, plot, ilon, ilev, polyres)

  ;   getvalues plot@contour                                                     
  ; "pmAnnoManagers" : am_ids1
                                         
  ; end getvalues                  
   ; printVarSummary(plot)   

 res2@cnLevelSelectionMode = "ExplicitLevels"
 res2@cnLinesOn=True
 res2@cnLevels = ispan(-6,6,1)
 res2@cnLineLabelsOn=True
 res2@cnLineLabelBackgroundColor="Transparent"
 res2@cnLineThicknessF = 3.
 res2@cnFillOn        = False
 res2@gsnContourNegLineDashPattern=1
 plot1 = gsn_csm_pres_hgt(wks2, msf_scale, res2)

;-----------------------------------------------------
 ;去掉height字体
;----------------------------------------------------
 
 ; overlay(plot,plot1)
 overlay(plot, iplot)


 draw(plot)
 frame(wks2)


end
   