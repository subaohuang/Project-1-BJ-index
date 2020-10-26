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
    plot = new(1,graphic)
    iplot= plot
    plot1 = plot


    a=6371000
    g=9.8
    pi=get_pi("double")
    dfi=(pi*10)/180

    wks = gsn_open_wks("pdf", "201026-20CR-东太哈德来环流的趋势")
    
;;;;
   
    yearstart = 190001
    yearend = 201012
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_20CR_187101-201212.nc"/)

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    msf = f1->msf(timeind,{1000:100},:,:)
    
    msf = (/msf/10^10/)
    printVarSummary(msf)

  
    msf_mean = dim_avg_n_Wrap(msf(:,:,:,{240:270}), (/3/))
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,{180:270}), (/3/))



    msf_mean_annual = month_to_annual(msf_mean, 1)
    msf_mean_clm = dim_avg_n_Wrap(msf_mean_annual, 0)


    dim = dimsizes(msf_mean_annual)

    ;;;计算趋势
    ; rc = msf_mean_annual(1,:,:)
    rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual(:,:,{-60:60}),0,0)

    copy_VarCoords(msf_mean_annual(1,:,{-60:60}), rc)

    rc  = (/rc*dim(0)/)


    prob := rc_test(rc)

    copy_VarCoords(msf_mean_annual(1,:,{-60:60}), prob)
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
    res@cnMinLevelValF = -3.0
    res@cnMaxLevelValF = 3.0
    res@cnLevelSpacingF = 0.5

    res@lbLabelBarOn =True
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    


    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0
  
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
 
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-60,60,10)

  
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7

    ;;;;
    res@gsnLeftString = "Linear trend of Regional Hadley Circulation in 20CR"
    res@gsnLeftStringFontHeightF = 0.015
    plot = gsn_csm_pres_hgt(wks,rc,res)


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

   
    iplot=gsn_csm_contour(wks,prob,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.05      ; make dots larger

    iplot = gsn_contour_shade(iplot, 0,0.1, opt)
    ;;;  
    overlay(plot, iplot)


    ;;;画气候态的等值线
    res@gsnRightString=""
    res@gsnLeftString=""
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = ispan(-8,12,2)
    res@cnLineLabelsOn=True
    res@cnLineLabelBackgroundColor="White"
    res@cnLineThicknessF = 2.
    res@cnFillOn = False
    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=1.0

    plot1 = gsn_csm_pres_hgt(wks, msf_mean_clm, res)
    ;;;  
    overlay(plot, plot1)

    draw(plot)
    frame(wks)
end 



