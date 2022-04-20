
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
 fils1=systemfunc("ls "+diri+"msf_time_era_global_10--10.nc")
 f1=addfile(fils1, "r")
 msf=f1->msf(:,:)
 ud=f1->ud(:,:)

 w=f1->w(:,:)
 level=f1->level(:)
 lon=f1->lon(:)
 msf_time=f1->msf_time
 ud_time=f1->ud_time
 w_time=f1->w_time
 msf&level=level
 msf&lon=lon
 msf&level@units="hpa"
 copy_VarCoords(msf, ud)
 copy_VarCoords(msf, w)
 copy_VarCoords(msf_time, ud_time)
 copy_VarCoords(msf_time, w_time)

 msf_annual=month_to_annual(msf_time(time|0:455,level|:,lon|:), 1)
 msf_annual!0="time"
 msf_annual_scale=msf_annual(level|:,lon|:,time|:)/10^9
 copy_VarCoords(msf_annual(level|:,lon|:,time|:), msf_annual_scale)
 printVarSummary(msf_annual_scale)
 delete(msf_annual)
 
 msf_scale=msf/10^9
 copy_VarMeta(msf, msf_scale)

 ud_annual_1=month_to_annual(ud_time(time|0:455,level|:,lon|:), 1)
 ud_annual_1!0="time"
 ud_annual=ud_annual_1(level|:,lon|:,time|:)
 delete(ud_annual_1) 
 w_annual_1=month_to_annual(w_time(time|0:455,level|:,lon|:), 1)
 w_annual_1!0="time"
 w_annual=w_annual_1(level|:,lon|:,time|:)
 delete(w_annual_1) 
 ; msf_time_scale=msf_time/10^9
 ; copy_VarMeta(msf_time, msf_time_scale)

 printVarSummary(msf_scale)
 msf_rc=decade_trend(msf_time,"units","msf_rc")

  msf_rc_scale=120*msf_rc/(10^9)
 copy_VarMeta(msf_rc, msf_rc_scale)
    
 printVarSummary(msf_rc_scale)
;;;趋势的显著性检验
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
   delete(i)
  ;;;
  
 
;;;;选择enino_modoki的年份
   cano_year=(/1997/)
   modoki1_year=(/1987,1990,1991,2002/)
   modoki2_year=(/1992,1994,2004/)

   printVarSummary(msf_annual_scale)
   msf_annual_scale&time=ispan(1979, 2016, 1)
   copy_VarCoords(msf_annual_scale, ud_annual)
   copy_VarCoords(msf_annual_scale, w_annual)

  
   msf_annual_scale_cano=msf_annual_scale(:,:,get1Dindex(msf_annual_scale&time,cano_year))
   msf_annual_scale_modoki1=dim_avg_n_Wrap(msf_annual_scale(:,:,get1Dindex(msf_annual_scale&time,modoki1_year)),2)
   msf_annual_scale_modoki2=dim_avg_n_Wrap(msf_annual_scale(:,:,get1Dindex(msf_annual_scale&time,modoki2_year)),2)
   
   ud_annual_scale_cano=ud_annual(:,:,get1Dindex(ud_annual&time,cano_year))
   ud_annual_scale_modoki1=dim_avg_n_Wrap(ud_annual(:,:,get1Dindex(ud_annual&time,modoki1_year)),2)
   ud_annual_scale_modoki2=dim_avg_n_Wrap(ud_annual(:,:,get1Dindex(ud_annual&time,modoki2_year)),2)

   w_annual_scale_cano=w_annual(:,:,get1Dindex(w_annual&time,cano_year))
   w_annual_scale_modoki1=dim_avg_n_Wrap(w_annual(:,:,get1Dindex(w_annual&time,modoki1_year)),2)
   w_annual_scale_modoki2=dim_avg_n_Wrap(w_annual(:,:,get1Dindex(w_annual&time,modoki2_year)),2)




;;;;---------------    
  wks=gsn_open_wks("pdf", ""+diri+"mass_stream_era_globle")

  wks2=gsn_open_wks("pdf", ""+diri+"msf_rc_era_globle")
  ;;;画图
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="Climatology"
    
    ; res@lbLabelStrings=ispan(-6,6,1)
    res@cnFillPalette  = "BlueDarkRed18"
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -4
    res@cnMaxLevelValF = 4
    res@cnLevelSpacingF = 0.5
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
 
    res@cnLineLabelsOn=True
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

    res@vcRefAnnoOrthogonalPosF  = -1.15
    res@vcRefAnnoPerimOn=False
    
    res@vpHeightF=0.3
    res@vpWidthF=0.7
    
    plot=gsn_csm_pres_hgt_vector(wks,msf_scale,ud,w,res)
 
    ;;;画38年的情况
    plot_aun=new(38,graphic)
   
    
   
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
  
  dum= gsn_add_polyline(wks, plot,  ploylon1,ploylevel1,  ployres)


   dum1=new(38,graphic)
   dum2=dum1
   dum3=dum1
   

 res_aun=res

;-------------------------------
;添加0线
;-----------------------
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn=True
res@cnLevels = 0
res@cnLineThicknessF = 4.
res@cnFillOn        = False

plot1 = gsn_csm_pres_hgt(wks, msf_scale, res)

plot1_anunal=new(38,graphic)

do i=0,37
   res_aun@gsnRightString=""+(/1979+i/)+""
   wks1=gsn_open_wks("png", ""+diri+"msf_w_anu/"+"msf_anu"+""+(/1979+i/)+"")   
   plot_aun(i)=gsn_csm_pres_hgt_vector( wks1,w_annual(:,:,i),ud_annual(:,:,i),w_annual(:,:,i),res_aun)   

   dum1(i)= gsn_add_polyline(wks1,plot_aun(i),  ploylon1,ploylevel1,  ployres)
   ployres@gsLineDashPattern=16
   dum2(i)= gsn_add_polyline(wks1,plot_aun(i),  ploylon2,ploylevel2,  ployres)
   ployres@gsLineDashPattern=0
   dum3(i)= gsn_add_polyline(wks1,plot_aun(i),  ploylon3,ploylevel3,  ployres)
   delete(res@gsnRightString)
   plot1_anunal(i) = gsn_csm_pres_hgt(wks1, msf_annual_scale(:,:,i), res)
   
  overlay(plot_aun(i), plot1_anunal(i))
  draw(plot_aun(i))
  frame(wks1)
  
 end do

overlay(plot,plot1)
draw(plot)
frame(wks)
delete(plot)
delete(plot1)
;;; 画3种时期的情况
   res_aun@lbLabelBarOn=False
   dum4=new(3,graphic)
   dum5=new(3,graphic)
   dum6=new(3,graphic)

   plot_com=new(3,graphic)
   wks2=gsn_open_wks("pdf", ""+diri+"msf_enino3种情况_10--10")
   res_aun@gsnRightString="canonical enino"
   plot_com(0)=gsn_csm_pres_hgt_vector( wks2,w_annual_scale_cano,ud_annual_scale_cano,w_annual_scale_cano,res_aun) 
   res_aun@gsnRightString="enino modoki1"
   plot_com(1)=gsn_csm_pres_hgt_vector( wks2,w_annual_scale_modoki1,ud_annual_scale_modoki1,w_annual_scale_modoki1,res_aun)
   res_aun@gsnRightString="enino modoki2"
   plot_com(2)=gsn_csm_pres_hgt_vector( wks2,w_annual_scale_modoki2,ud_annual_scale_modoki2,w_annual_scale_modoki2,res_aun)
   
    ;;;请添加0线
    plot1=new(3,graphic)
   
    plot1(0) = gsn_csm_pres_hgt(wks2, msf_annual_scale_cano, res)
    plot1(1) = gsn_csm_pres_hgt(wks2, msf_annual_scale_modoki1, res)
    plot1(2) = gsn_csm_pres_hgt(wks2, msf_annual_scale_modoki2, res)
  

    do i=0,2
    dum4(i)= gsn_add_polyline(wks2,plot_com(i),  ploylon1,ploylevel1,  ployres)
    end do
    ployres@gsLineDashPattern=16
    do i=0,2
    dum5(i)= gsn_add_polyline(wks2,plot_com(i),  ploylon2,ploylevel2,  ployres)
    end do
    ployres@gsLineDashPattern=0
    do i=0,2
    dum6(i)= gsn_add_polyline(wks2,plot_com(i),  ploylon3,ploylevel3,  ployres)
    end do
    do i=0,2
   overlay(plot_com(i),plot1(i))
 end do
   pres=True
   pres@gsnPanelLabelBar=True
   gsn_panel(wks2, plot_com, (/3,1/), pres)
exit()
;----------------------------------------------------------------------
;画随时间变化
;----------------------------------------------------------------------
   gsn_define_colormap(wks2,"BlueWhiteOrangeRed")
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

    plot=gsn_csm_pres_hgt(wks2,msf_rc_scale,res2)
  ; ;;;通过显著性的点
  ;    polyres                   = True
  ;    polyres@gsMarkerIndex     = 16          ; polymarker style
  ;    polyres@gsMarkerSizeF     = 3.         ; polymarker size

  ;   dum1=gsn_add_polymarker(wks2, plot, ilon, ilev, polyres)

;;;添加点，gsn_contour_shade
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
 res2@cnLevels = ispan(-6,6,1)
 res2@cnLineLabelsOn=True
 res2@cnLineLabelBackgroundColor="Transparent"
 res2@cnLineThicknessF = 3.
 res2@cnFillOn        = False
 res2@gsnContourNegLineDashPattern=1

 plot1 = gsn_csm_pres_hgt(wks2, msf_scale, res2)
;;;添加地形
 
  ployres@gsLineDashPattern=0
  dum1= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)
  ployres@gsLineDashPattern=16
  dum2= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)
  ployres@gsLineDashPattern=0
  dum3= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)

  


 overlay(plot, iplot)
 overlay(plot,plot1)
 draw(plot)
 frame(wks2)


end
   