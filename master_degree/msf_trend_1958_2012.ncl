function decade_trend(data)
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
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
begin
    diri="/mnt/f/share/"
    fils1=systemfunc("ls "+diri+"ncep/20CR/msf_time_20cr_global_6--6.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    
    msf_time_20cr=f1->msf_time(:,:,timeind)
    ud_time_20cr=f1->ud_time(:,:,timeind)
    w_time_20cr=f1->w_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)
     
    printVarSummary(ud_time_20cr)
    printVarSummary(msf_time_20cr)

    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)

    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)
;;;

;;;
    fils1=systemfunc("ls "+diri+"ncep/ncep1/msf_time_ncep1_global.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    
    msf_time_ncep1=f1->msf_time(:,:,timeind)
    ud_time_ncep1=f1->ud_time(:,:,timeind)
    w_time_ncep1=f1->w_time(:,:,timeind)
    msf_time_ncep1=(/msf_time_ncep1/10^9/)
    
    delete(msf_time_ncep1&time)
    msf_time_ncep1&time=date(timeind)


    delete(fils1)
    delete(f1)
    delete(date)
    delete(timeind)
;;;

   msf_clm_20cr=dim_avg_n_Wrap(msf_time_20cr, 2)
   msf_clm_ncep1=dim_avg_n_Wrap(msf_time_ncep1, 2)



;;;计算趋势 

    rc1=decade_trend(msf_time_20cr(:,:,{195801:197812}))

    rc2=decade_trend(msf_time_20cr(:,:,{197901:201212}))

    rc3=decade_trend(msf_time_20cr(:,:,{195801:201012}))

    rc1=(/rc1*120/)
    rc2=(/rc2*120/)
    rc3=(/rc3*120/)

    rc4=decade_trend(msf_time_ncep1(:,:,{195801:197812}))

    rc5=decade_trend(msf_time_ncep1(:,:,{197901:201212}))

    rc6=decade_trend(msf_time_ncep1(:,:,{195801:201212}))

    rc4=(/rc4*120/)
    rc5=(/rc5*120/)
    rc6=(/rc6*120/)





;;;

;;;显著性检验
    prob1=rc_test(rc1)
    prob2=rc_test(rc2)
    prob3=rc_test(rc3)

    prob4=rc_test(rc4)
    prob5=rc_test(rc5)
    prob6=rc_test(rc6)
;;;

;;;画图

    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/20cr&ncep1的msf的趋势")
    ; gsn_define_colormap(wks,"GMT_polar")
    gsn_define_colormap(wks, "NCV_blu_red")

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString="20CR"
    res2@gsnRightString="1958-1978"
   
    res2@cnLinesOn=False

    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -.06
    res2@cnMaxLevelValF = 0.06
    res2@cnLevelSpacingF = 0.01

    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
  
    res2@lbOrientation="Horizontal"
    res2@lbLabelBarOn=False

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
   
   plot=new(6,graphic)
    
    plot(0)=gsn_csm_pres_hgt(wks,rc1,res2)
    res2@gsnLeftString="NCEP1"
    plot(1)=gsn_csm_pres_hgt(wks,rc4,res2)
    res2@gsnRightString="1979-2012"
    res2@gsnLeftString="20CR"
    plot(2)=gsn_csm_pres_hgt(wks,rc2,res2)
    res2@gsnLeftString="NCEP1"
    plot(3)=gsn_csm_pres_hgt(wks,rc5,res2)
    res2@gsnRightString="1958-2012"
    res2@gsnLeftString="20CR"
    plot(4)=gsn_csm_pres_hgt(wks,rc3,res2)
    res2@gsnLeftString="NCEP1"
    plot(5)=gsn_csm_pres_hgt(wks,rc6,res2)
    res2@gsnRightString=""

;;;gsn_contour_shade
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@gsnRightString=""
  res_pro@gsnLeftString=""
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
    iplot(0)=gsn_csm_pres_hgt(wks,prob1,res_pro)
    iplot(1)=gsn_csm_pres_hgt(wks,prob4,res_pro)
    iplot(2)=gsn_csm_pres_hgt(wks,prob2,res_pro)

    iplot(3)=gsn_csm_pres_hgt(wks,prob5,res_pro)
    iplot(4)=gsn_csm_pres_hgt(wks,prob3,res_pro)
    iplot(5)=gsn_csm_pres_hgt(wks,prob6,res_pro)



  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

    iplot(0) = gsn_contour_shade(iplot(0), 0,0.05, opt)
    iplot(1) = gsn_contour_shade(iplot(1), 0,0.05, opt)
    iplot(2) = gsn_contour_shade(iplot(2), 0,0.05, opt)

    iplot(3) = gsn_contour_shade(iplot(3), 0,0.05, opt)
    iplot(4) = gsn_contour_shade(iplot(4), 0,0.05, opt)
    iplot(5) = gsn_contour_shade(iplot(5), 0,0.05, opt)
 ;;;

 res2@cnLevelSelectionMode = "ExplicitLevels"
 res2@cnLinesOn=True
 res2@cnLevels = ispan(-6,6,1)
 res2@cnLineLabelsOn=True
 res2@cnLineLabelBackgroundColor="Transparent"
 res2@cnLineThicknessF = 1.
 res2@cnFillOn = False
 res2@gsnContourNegLineDashPattern=1
 res2@gsnContourZeroLineThicknessF=2.0
 res2@gsnLeftString=""

 plot1=new(6,graphic)
 do i=0,4,2
 plot1(i) = gsn_csm_pres_hgt(wks, msf_clm_20cr, res2)
 end do

 do i=1,5,2
 plot1(i) = gsn_csm_pres_hgt(wks, msf_clm_ncep1, res2)
 end do
;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 8.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 10)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
 
    dum1=new(6,graphic)
    dum2=dum1
    dum3=dum1

    do i=0,5
    ployres@gsLineDashPattern=0
    dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
    end do 

    do i=0,5
    overlay(plot(i),iplot(i))
    overlay(plot(i),plot1(i))
    end do
    pres=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.8
    gsn_panel(wks, plot, (/3,2/), pres)
end



