undef("rmtrend") ;;;滑动趋势
function rmtrend(x,nave,year)
begin
length=dimsizes(x&year)
rc=new((/length-nave+1,dimsizes(x&lon)/),typeof(x))
rc_test=new((/length-nave+1,dimsizes(x&lon)/),typeof(x))

siglv=new(dimsizes(x&lon),float)

do t=nave,length  
   ssta_rc=regCoef_n(ispan(1,nave,1),x(t-nave:t-1,:),0,0)
   tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval                                   ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(ntime,nlon)
  
  rc_test(t-nave,:)=prob 
  rc(t-nave,:)=ssta_rc
end do 

rc!0="time"
rc!1="lon"
; rc&time=ispan(year+(nave-1)/2, (year+(nave-1)/2+length-nave), 1)
rc&time=ispan(year, year+length-nave, 1)
rc&lon=x&lon
rc@rc_test=rc_test
return (rc)
end

  undef("fit") ;;;得到回归的值
  function fit(slp_20cr_mean,ipindex1)
  begin
  rc_slp_20cr=regCoef_n(ipindex1, slp_20cr_mean ,0, 0)
  slp_20cr_re=new(dimsizes(slp_20cr_mean),typeof(slp_20cr_mean))
  do i=0,dimsizes(slp_20cr_mean&lon)-1,1
  slp_20cr_re(:,i)=rc_slp_20cr(i)*ipindex1+rc_slp_20cr@yintercept(i)
  end do 
  copy_VarCoords(slp_20cr_mean, slp_20cr_re)
  return(slp_20cr_re)
  end




begin
  diri="/mnt/f/share/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"ncep/20CR/uwnd.mon.mean.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.195801.and.date.le.201212)

  slp_20cr=f1->uwnd(timeind,{200:500},{-5:5},:)

  slp_20cr_annual=month_to_annual(slp_20cr,1)
  ; delete(msf_time_20cr&time)
  ; msf_mean_20cr&time=date(timeind)
  slp_20cr_mean=dim_avg_n_Wrap(slp_20cr_annual, (/1,2/))
  slp_20cr_mean&year=ispan(1,55,1)  
  printVarSummary(slp_20cr_mean)

  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)
;;;


;;;
  fils1=systemfunc("ls "+diri+"ncep/ncep1/uwnd.mon.mean.nc")
  ; fils1=systemfunc("ls "+diri+"ncep/ncep1/ud_ncep1_global.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.195801.and.date.le.201212)
 
  slp_ncep1_r=f1->uwnd(timeind,{1000},{5:-5},:)
   ; slp_ncep1_r=f1->ud_new({500:200},{-5:5},:,timeind)
  slp_ncep1=slp_ncep1_r(time|:,lat|:,lon|:)
   printVarSummary(slp_ncep1)
 
  slp_ncep1_annual=month_to_annual(slp_ncep1, 1)
  ; slp_ncep1_mean=dim_avg_n_Wrap(slp_ncep1_annual, (/1,2/))
  slp_ncep1_mean=dim_avg_n_Wrap(slp_ncep1_annual, 1)
  slp_ncep1_mean&year=ispan(1,55,1)  
   ; delete(msf_time_ncep1&time)
   ; msf_mean_ncep1&time=date(timeind)
  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)
;;;
ipindex1=asciiread("/mnt/f/share/era/data/GHCN_HAD_p_ud_timef1_1958_2012.ascii", 55, "float")
ipindex2=asciiread("/mnt/f/share/era/data/GHCN_HAD_p_ud_timef2_1958_2012.ascii", 55, "float")
;;;20cr
  ;;;回归，得到回归的值
  slp_20cr_re1=fit(slp_20cr_mean,ipindex1)
  slp_20cr_re2=fit(slp_20cr_mean,ipindex2)
  ;;用回归的值计算滑动趋势
  year=30
  slp_re_trend_20cr1=rmtrend(slp_20cr_re1,year,1958)
  slp_re_trend_20cr2=rmtrend(slp_20cr_re2,year,1958)
  slp_re_trend_20cr3=slp_re_trend_20cr1+slp_re_trend_20cr2
  copy_VarMeta(slp_re_trend_20cr1, slp_re_trend_20cr3)
  slp_trend_20cr=rmtrend(slp_20cr_mean,year,1958)

  slp_re_trend_20cr1=(/slp_re_trend_20cr1*year/)
  slp_re_trend_20cr2=(/slp_re_trend_20cr2*year/) 
  slp_re_trend_20cr3=(/slp_re_trend_20cr3*year/) 
  slp_trend_20cr=(/slp_trend_20cr*year/) 
;;

;;;ncep1
  ;;;回归，得到回归的值
  slp_ncep1_re1=fit(slp_ncep1_mean,ipindex1)
  slp_ncep1_re2=fit(slp_ncep1_mean,ipindex2)
  ;;用回归的值计算滑动趋势
 
  slp_re_trend_ncep11=rmtrend(slp_ncep1_re1,year,1958)
  slp_re_trend_ncep12=rmtrend(slp_ncep1_re2,year,1958)
  slp_re_trend_ncep13=slp_re_trend_ncep11+slp_re_trend_ncep12
  copy_VarMeta(slp_re_trend_ncep11, slp_re_trend_ncep13)
  slp_trend_ncep1=rmtrend(slp_ncep1_mean,year,1958)

  slp_re_trend_ncep11=(/slp_re_trend_ncep11*year/)
  slp_re_trend_ncep12=(/slp_re_trend_ncep12*year/) 
  slp_re_trend_ncep13=(/slp_re_trend_ncep13*year/) 
  slp_trend_ncep1=(/slp_trend_ncep1*year/) 
;;
;;;检验 
  ; rc_test1=slp_trend_ncep1@rc_test
  ; copy_VarCoords(slp_trend_ncep1, rc_test1)
  ; rc_test2=slp_re_trend_ncep12@rc_test
  ; copy_VarCoords(slp_trend_ncep1, rc_test2)

  ; rc_test3=slp_re_trend_ncep11@rc_test
  ; copy_VarCoords(slp_trend_ncep1, rc_test3)


;;;画图
wks  = gsn_open_wks ("pdf", ""+diri_pi+"EOF回归的uwnd_ncep1_1000hpa_u30年滑动趋势" )           ; send graphics to PNG file
 ; gsn_define_colormap(wks,"BlAqGrWh2YeOrReVi22")
 ; gsn_define_colormap(wks, "BlRe")
  gsn_define_colormap(wks, "BlueWhiteOrangeRed")

  res            = True                         ; plot mods desired
  res@gsnDraw    = False                        ; don't draw
  res@gsnFrame   = False                        ; don't advance frame
  res@gsnLeftString="20CR "
  res@gsnRightString="u trend"

  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@cnFillOn=True
  res@cnFillDrawOrder="PreDraw"
  res@cnLinesOn=False
  res@cnFillMode="RasterFill"
  res@cnRasterSmoothingOn=True
  res@cnInfoLabelOn=False
  res@cnLineLabelsOn=False

  res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
  res@cnMinLevelValF       = -2.0               ; min level
  res@cnMaxLevelValF       =  2.0             ; max leve
  res@cnLevelSpacingF      =0.25
  res@lbLabelStride=2

  ; res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  ; res@tmYLLabelsOn=True
  ; res@tmYUseLeft=True

  ; res@tmYLValues    = ispan(0,39,5)
  ; res@tmYLLabels    = ispan(1958,1998,5)

  res@tmXTOn=False
  res@tmYROn=False

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues=(/0,90,180,270,357.5/)
  res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
  res@tmXBMinorOn=True
  res@tmXBMinorValues=ispan(0,360,10)
  
  plot=new(4,graphic)
  plot(0)= gsn_csm_hov(wks, slp_trend_ncep1, res)   ; create plot
   res@gsnRightString="EOF1+EOF2 fit u trend"
  plot(1)= gsn_csm_hov(wks, slp_re_trend_ncep13, res)   ; create plot
   res@gsnRightString="EOF2 fit u trend"
  plot(2)= gsn_csm_hov(wks, slp_re_trend_ncep12, res)   ; create plot
   res@gsnRightString="EOF1 fit u trend"
  plot(3)= gsn_csm_hov(wks, slp_re_trend_ncep11, res)   ; create plot
 

   ;;;打点
   
  ;  ;;;打点
  ; res_pro=True
  ; res_pro@gsnDraw=False
  ; res_pro@gsnFrame=False
  ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
  ; res_pro@cnLinesOn=False
  ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  ; res_pro@cnLevelSpacingF     = 0.01       ; set contour spacing
  ; ; res_pro@cnLevels = ispan(-6,6,1)
  ; res_pro@cnLineLabelsOn=False
  ; res_pro@tmYRMode = "Automatic"
  ; res_pro@cnInfoLabelOn=False
  ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; ; res_pro@cnLineThicknessF = 3.
  ; ; res_pro@cnFillOn        = True

  ; iplot=new(3,graphic)
   
  ; iplot(0)=gsn_csm_contour(wks,rc_test1,res_pro)   ;;rc_test1是对原来进行检验
  ; iplot(1)=gsn_csm_contour(wks,rc_test3,res_pro)   ;;rc_test2是对ipo进行检验
  ; iplot(2)=gsn_csm_contour(wks,rc_test2,res_pro)   ;;rc_test3是对全球增温进行检验

  ; opt                      = True
  ; opt@gsnShadeFillType     = "pattern"
  ; opt@gsnShadeMid     = 17         ; stipple
  ; ; opt@gsnShadeHigh   = "green"
  ; opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  ; opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
  
  ; do i=0,2
  ; iplot(i) = gsn_contour_shade(iplot(i), 0,0.01, opt)
  ; end do
  ; ;;;  
  
 ;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 5.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=15
   ;---Attach the polyline
    ploylevel=1958
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 25)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
      
    dum1=new(4,graphic)
    dum2=dum1
    dum3=dum1
    do i=0,3  
    dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
    end do
    overlay(plot(0), iplot(0))
    overlay(plot(2), iplot(2))
    overlay(plot(3), iplot(1))
    pres=True
    pres@gsnPanelLabelBar=True
    gsn_panel(wks, plot, (/2,2/), pres)
  
end





