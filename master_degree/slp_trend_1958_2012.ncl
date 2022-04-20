function decade_trend(data)
    local data_a, data_a_new, rc
    begin
    time=ispan(1,dimsizes(data&time),1)
    rc = regCoef_n(time,data,0,0)
    rc!0="lat"
    rc!1="lon"
    rc&lat=data&lat
    rc&lon=data&lon
    rc&lat@units = "degrees_north"
    rc&lon@units = "degrees_east"
    ; print(rc)
return(rc)
end
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end

;************************************************
; 画20cr,ncep,era的msf的变化趋势。
;************************************************
begin
    diri="/mnt/f/share/"
    fils1=systemfunc("ls "+diri+"ncep/20CR/prmsl.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"ncep/20CR/uwnd.mon.mean.nc")
    fils3=systemfunc("ls "+diri+"ncep/20CR/vwnd.mon.mean.nc")
    

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f3=addfile(fils3, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)

    slp_20cr=f1->prmsl(timeind,{20:-20},:)
    slp_20cr=(/slp_20cr/100.0/)
    uwnd_20cr=f2->uwnd(timeind,{850},{20:-20},:)
    vwnd_20cr=f3->vwnd(timeind,{850},{20:-20},:)

    delete(slp_20cr&time)
    delete(uwnd_20cr&time)
    delete(vwnd_20cr&time)
    slp_20cr&time=date(timeind)
    uwnd_20cr&time=date(timeind)
    vwnd_20cr&time=date(timeind)

    delete(fils1)
    delete(f1)
    delete(fils2)
    delete(f2)
    delete(fils3)
    delete(f3)
    delete(date)
    delete(timeind)
;;;


;;;
    fils1=systemfunc("ls "+diri+"ncep/ncep1/slp.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"ncep/ncep1/uwnd.mon.mean.nc")
    fils3=systemfunc("ls "+diri+"ncep/ncep1/vwnd.mon.mean.nc")

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    f3=addfile(fils3, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    
    slp_ncep1=f1->slp(timeind,{20:-20},:)
    uwnd_ncep1=f2->uwnd(timeind,{850},{20:-20},:)
    vwnd_ncep1=f3->vwnd(timeind,{850},{20:-20},:)

    delete(slp_ncep1&time)
    delete(uwnd_ncep1&time)
    delete(vwnd_ncep1&time)
    slp_ncep1&time=date(timeind)
    uwnd_ncep1&time=date(timeind)
    vwnd_ncep1&time=date(timeind)

    delete(fils1)
    delete(f1)
    delete(fils2)
    delete(f2)
    delete(fils3)
    delete(f3)
    delete(date)
    delete(timeind)
;;;

   

;;;计算趋势 

    rc1=decade_trend(slp_20cr({195801:197812},:,:))
    rc2=decade_trend(slp_20cr({197901:201212},:,:))
    rc3=decade_trend(slp_20cr({195801:201212},:,:))

    rc1=(/rc1*120/)
    rc2=(/rc2*120/)
    rc3=(/rc3*120/)

    rc1_u=decade_trend(uwnd_20cr({195801:197812},:,:))
    rc2_u=decade_trend(uwnd_20cr({197901:201212},:,:))
    rc3_u=decade_trend(uwnd_20cr({195801:201212},:,:))

    rc1_u=(/rc1_u*120/)
    rc2_u=(/rc2_u*120/)
    rc3_u=(/rc3_u*120/)

    rc1_v=decade_trend(vwnd_20cr({195801:197812},:,:))
    rc2_v=decade_trend(vwnd_20cr({197901:201212},:,:))
    rc3_v=decade_trend(vwnd_20cr({195801:201212},:,:))

    rc1_v=(/rc1_v*120/)
    rc2_v=(/rc2_v*120/)
    rc3_v=(/rc3_v*120/)

;;;
   
    rc4=decade_trend(slp_ncep1({195801:197812},:,:))
    rc5=decade_trend(slp_ncep1({197901:201212},:,:))
    rc6=decade_trend(slp_ncep1({195801:201212},:,:))

    rc4=(/rc4*120/)
    rc5=(/rc5*120/)
    rc6=(/rc6*120/)
 
    rc4_u=decade_trend(uwnd_ncep1({195801:197812},:,:))
    rc5_u=decade_trend(uwnd_ncep1({197901:201212},:,:))
    rc6_u=decade_trend(uwnd_ncep1({195801:201212},:,:))

    rc4_u=(/rc4_u*120/)
    rc5_u=(/rc5_u*120/)
    rc6_u=(/rc6_u*120/)

    rc4_v=decade_trend(vwnd_ncep1({195801:197812},:,:))
    rc5_v=decade_trend(vwnd_ncep1({197901:201212},:,:))
    rc6_v=decade_trend(vwnd_ncep1({195801:201212},:,:))

    rc4_v=(/rc4_v*120/)
    rc5_v=(/rc5_v*120/)
    rc6_v=(/rc6_v*120/)
;;;

;;;显著性检验
    prob1=rc_test(rc1)
    prob2=rc_test(rc2)
    prob3=rc_test(rc3)
   
    prob1_u=rc_test(rc1_u)
    prob2_u=rc_test(rc2_u)
    prob3_u=rc_test(rc3_u)

    prob4=rc_test(rc4)
    prob5=rc_test(rc5)
    prob6=rc_test(rc6)

    prob4_u=rc_test(rc4_u)
    prob5_u=rc_test(rc5_u)
    prob6_u=rc_test(rc6_u)
;;;

;;;没通过检验的风场不画上去
   rc1_u=where(prob1_u.le.0.05, rc1_u, rc1_u@_FillValue)
   rc2_u=where(prob2_u.le.0.05, rc2_u, rc2_u@_FillValue)
   rc3_u=where(prob3_u.le.0.05, rc3_u, rc3_u@_FillValue)
   
   rc4_u=where(prob4_u.le.0.05, rc4_u, rc4_u@_FillValue)
   rc5_u=where(prob5_u.le.0.05, rc5_u, rc5_u@_FillValue)
   rc6_u=where(prob6_u.le.0.05, rc6_u, rc6_u@_FillValue)
;;;





;;;

;;;画图

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/slp_850hpa风场的趋势_1958-2012")
    gsn_define_colormap(wks, "GMT_polar")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
    res@gsnLeftString="20CR"
    res@gsnRightString="1958-1978"
    
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnFillMode="RasterFill"
    res@cnRasterSmoothingOn=True
    res@gsnScalarContour        = True
 
    res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"transparent","transparent"/) ;填充位置相应的具体颜色
    res@mpOutlineOn = True; turn off outline 



    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(rc1&lat)
    res@mpMaxLatF             = max(rc1&lat)
    res@mpMinLonF             = min(rc1&lon)
    res@mpMaxLonF             = max(rc1&lon)
    res@mpCenterLonF          =avg(rc1&lon)


    res@tmYLMode="Explicit"
    res@tmYLValues=ispan(-20,20,10)
    res@tmYLLabels=(/"20~S~o~N~S","10~S~o~N~S","0~S~o~N","10~S~o~N~N","20~S~o~N~N"/)
    res@tmYLLabelFontHeightF=0.01

    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(0,360,30)
    res@tmXBLabelFontHeightF=0.01
    res@tmXBLabelFontThicknessF=0.02
    res@tmXBMinorOn=True
    res@tmXBLabels=(/"0","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
                "180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0"/)

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.8
    res@vpHeightF = 0.2

    
    res@vcRefMagnitudeF          = 0.8          ; define vector ref mag
    res@vcGlyphStyle="LineArrow"
    res@vcRefLengthF             = 0.03         ; define length of vec ref
    res@vcMinFracLengthF        = 0.3
    res@vcMinDistanceF           = 0.016          ; thin vector
    res@vcLineArrowColor="green"
    res@vcLineArrowThicknessF=1.5
    ; res@vcFillArrowMinFracWidthF=1.0  ;;;所有箭头调成一样的大小
    ; res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadMinFracXF=1.0
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上角

 plot=new(6,graphic)
    
    plot(0)=gsn_csm_vector_scalar_map(wks, rc1_u, rc1_v, rc1, res)
    res@gsnLeftString="NCEP1"
    plot(1)=gsn_csm_vector_scalar_map(wks, rc4_u, rc4_v, rc4, res)
    res@gsnRightString="1979-2012"
    res@gsnLeftString="20CR"
    plot(2)=gsn_csm_vector_scalar_map(wks, rc2_u, rc2_v, rc2, res)
    res@gsnLeftString="NCEP1"
    plot(3)=gsn_csm_vector_scalar_map(wks, rc5_u, rc5_v, rc5, res)
    res@gsnRightString="1958-2012"
    res@gsnLeftString="20CR"
    res@vcRefMagnitudeF          = 0.6          ; define vector ref mag
    plot(4)=gsn_csm_vector_scalar_map(wks, rc3_u, rc3_v, rc3, res)
    res@gsnLeftString="NCEP1"
    res@vcRefMagnitudeF          = 0.6         ; define vector ref mag
    plot(5)=gsn_csm_vector_scalar_map(wks, rc6_u, rc6_v, rc6, res)
    res@gsnRightString=""


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(6,graphic)

    iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob4,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob2,res_pro)
    iplot(3)=gsn_csm_contour(wks,prob5,res_pro)
    iplot(4)=gsn_csm_contour(wks,prob3,res_pro)
    iplot(5)=gsn_csm_contour(wks,prob6,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

    do i=0,5 
    iplot(i) = gsn_contour_shade(iplot(i), 0,0.05, opt)
    end do
    ;;;  
    do i=0,5
    overlay(plot(i), iplot(i))
    end do
    pres=True
    pres@gsnPanelLabelBar=True
    gsn_panel(wks, plot, (/3,2/), pres)
    end