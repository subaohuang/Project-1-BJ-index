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

 begin
 diri="/mnt/f/share/JRA-55/"
 fils1=systemfunc("ls "+diri+"Q_time_JRA55.nc")
 f1=addfile(fils1, "r")
 time=f1->time
 date=cd_calendar(time, 1)
 timeind=ind(date.ge.197901.and.date.le.201212)
 Q=f1->Q_time(:,:,timeind)

year1=1979
year2=2012
Q_pi=lonPivot(Q(time|:,level|:,lon|:), 30)
    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    Q_djf=month_to_season(Q_pi, "DJF")
    delete(Q_djf&time)
    Q_djf&time=ispan(year1,year2,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;;
    delete(Q)
    ; Q_djf=(/month_to_annual(Q_pi, 1)/)
    Q_clm=dim_avg_n_Wrap(Q_djf, 0)
    ;;;计算线性趋势
    Q_rc=regCoef_n(ispan(1,dimsizes(Q_djf&time),1), Q_djf, 0, 0)
    copy_VarCoords(Q_djf(1,:,:), Q_rc)
    prob1=rc_test(Q_rc)

    Q_rc=(/Q_rc*10/)


  wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/Q_20CR的气候态")
 ;************************************************
; 画图气候态
;************************************************
;************************************************
; 画图气候态
;************************************************


    res = True
    res@gsnDraw=False
    res@gsnFrame=False

   
    res@gsnMaximize=False

    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.6
    res@cnMaxLevelValF = 0.6
    res@cnLevelSpacingF = 0.1
    res@cnFillOn = True 

    res@cnLineLabelsOn=False
    res@cnLineLabelInterval=2
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False


    res@lbLabelBarOn=True
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"

    res@tmXBValues=(/30,90,150,210,270,330,387.5/)
    res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
                    ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(30,390,10)
     
    res@tiYAxisString=""
    res@vpHeightF=0.3
    res@vpWidthF=0.7
   
      plot=gsn_csm_pres_hgt(wks, Q_rc, res)
    ; plot=gsn_csm_pres_hgt(wks, msf_trend(:,:), res)

 

    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -2
    res@cnMaxLevelValF = 2
    res@cnLevelSpacingF = 1
    res@cnLinesOn=True
    res@cnFillOn=False
    res@gsnContourPosLineDashPattern=0
    res@gsnContourNegLineDashPattern=1
    plot1=gsn_csm_pres_hgt(wks,Q_clm,res)

overlay(plot, plot1)

;;;打点
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  
  iplot=new(1,graphic)

  iplot(0)=gsn_csm_pres_hgt(wks,prob1,res_pro)


  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 1.0         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

  
  iplot= gsn_contour_shade(iplot, 0,0.1, opt)

;;;  

overlay(plot, iplot)

 ;;---Resources for the polyline  添加下垫面的陆地
      ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=1000
      ploylon1=fspan(9, 30, 25)+360.0
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylon6=fspan(30, 39, 25)


      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      ploylevel6=conform_dims(dimsizes(ploylon6), ploylevel, 0)
    

    str1 = unique_string("dum1")  ; "unique_string" will return a unique                                     ; string every time it is called fro                                 ;  within a single NCL session.
    str2 = unique_string("dum2")
    str3 = unique_string("dum3")
    str31 = unique_string("dum31")
    ployres@gsLineDashPattern=0
     
    plot@$str1$= gsn_add_polyline(wks,plot,  ploylon1,ploylevel1,  ployres)

    ployres@gsLineDashPattern=16

    plot@$str2$= gsn_add_polyline(wks,plot,  ploylon2,ploylevel2,  ployres)

    ployres@gsLineDashPattern=0
    plot@$str3$= gsn_add_polyline(wks,plot,  ploylon3,ploylevel3,  ployres)
    ployres@gsLineDashPattern=0
    plot@$str31$= gsn_add_polyline(wks,plot,  ploylon6,ploylevel6,  ployres)

draw(plot)
frame(wks)
end 


