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

    
;;;;
   
    yearstart = 197901
    yearend = 201512
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_ERAI_197901-201612.nc"/)

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     


    msf = f1->msf(timeind,{1000:100},:,:)
    
    msf = (/msf/10^10/)
    printVarSummary(msf)


    msf_mean = dim_avg_n_Wrap(msf(:,:,:,{125:160}), (/3/))



    season = (/"MAM","JJA","SON","DJF"/)

    title = (/"(a)","(b)","(c)","(d)"/)

    wks=gsn_open_wks("pdf", "local-msf-HC-trend")
    ; wks1=gsn_open_wks("pdf", "local-hc-clm")

    plot=new(4,graphic)
    ; plot1=new(4,graphic)
    str = ispan(1,4,1)+""
    iplot = plot
    plot1 = plot

    dum1=new(4,graphic)

do i = 0,3 

    ; msf_mean_annual = month_to_annual(msf_mean, 1)
    msf_mean_annual = month_to_season(msf_mean, season(i))
    msf_mean_clm = dim_avg_n_Wrap(msf_mean_annual, 0)


    dim = dimsizes(msf_mean_annual)

    ;;;计算趋势
    ; rc = msf_mean_annual(1,:,:)
    rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual(:,:,{-40:40}),0,0)

    copy_VarCoords(msf_mean_annual(1,:,{-40:40}), rc)

    rc  = (/rc*dim(0)/)


    prob := rc_test(rc)

    copy_VarCoords(msf_mean_annual(1,:,{-40:40}), prob)
    prob&level@units = "hPa"
    


;;;;;draw 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6.0
    res@cnMaxLevelValF = 6.0
    res@cnLevelSpacingF = 1.0

    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    


    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0
  
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
 
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(-40,40,10)
     
  
    res@tmYRMode="Automatic"                                                            
     
    res@vpHeightF=0.35
    res@vpWidthF=0.6

    ;;;;
    res@gsnLeftString = title(i) + " " + season(i) + " local HC trend"
    res@gsnLeftStringFontHeightF = 0.02
    plot(i) = gsn_csm_pres_hgt(wks,rc,res)



    ;;;画0-20的加粗直线
    ;---Resources for the polyline  添加下垫面的陆地
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 10.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=0
    ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(0, 20, 10)
    ; ploylon2=fspan(280, 290, 11)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ; ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)

    ; dum2=new(4,graphic) 
    dum1(i)= gsn_add_polyline(wks,plot(i), ploylon1,ploylevel1,  ployres)



  
    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1     ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

   
    iplot(i) = gsn_csm_contour(wks,prob,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.05      ; make dots larger

    iplot(i) = gsn_contour_shade(iplot(i), 0,0.1, opt)
    ;;;  
    overlay(plot(i), iplot(i))


    ;;;画气候态的等值线
    res@gsnRightString=""
    res@gsnLeftString=""
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLineLabelsOn=True
    res@cnLineLabelInterval = 1
     res@cnExplicitLineLabelsOn  = True

    res@cnLineLabelBackgroundColor="White"
    res@cnLineThicknessF = 2.
    res@cnFillOn = False
    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=1.0

    if i.eq.0.or.i.eq.2 then 
    res@cnLevels := ispan(-20,16,4)
    res@cnLineLabelStrings := ispan(-20,16,4) + ""
    end if 
    if i.eq.1 then 
    res@cnLevels := ispan(-28,16,4)
    res@cnLineLabelStrings := ispan(-28,16,4) + ""
    end if 
    if i.eq.3 then  
    res@cnLevels := ispan(-28,28,4)
    res@cnLineLabelStrings := ispan(-28,28,4) + ""
    end if

    
    plot1(i) = gsn_csm_pres_hgt(wks, msf_mean_clm, res)
    ;;;  
    overlay(plot(i), plot1(i))

end do 


pres  = True
pres@gsnPanelLabelBar = True
; drawNDCGrid(wks)


gsn_panel(wks, plot, (/2,2/), pres)

end 


