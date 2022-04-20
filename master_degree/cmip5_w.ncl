begin
	diri="/mnt/f/share/era/data/"
	fils1=systemfunc("ls /mnt/f/share/cmip/wap_Amon_NorESM1-M_amip_r1i1p1_197901-200512.nc") 
	f1=addfile(fils1, "r")
	; date=cd_calendar(f1->time, 1)
	; timeind=ind(date.le.201212.and.date.ge.195801)
    level=f1->plev
    levelind=ind(level.eq.50000)
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.le.200512.and.date.ge.197901)
	p=f1->wap(:,levelind,{5:-5},:)
    p_mean=dim_avg_n_Wrap(p, 1)
    printVarSummary(p_mean)
    p_annual=month_to_annual(p_mean, 1)



    wks=gsn_open_wks("pdf", ""+diri+"cmip5 w的气候态")

    res            = True                         ; plot mods desired
    res@gsnDraw    = False                        ; don't draw
    res@gsnFrame   = False                        ; don't advance frame
    res@cnLineLabelsOn =False
    res@cnFillDrawOrder="PreDraw"
    res@cnLinesOn=False
    res@cnFillPalette="testcmap"
    ; res@cnFillPalette="BlueDarkRed18"
    ; res@cnFillPalette="testcmap"
    ; res@cnFillPalette="MPL_PRGn"
    ; res@cnFillPalette="WhViBlGrYeOrReWh"
    ; res@cnFillPalette="NCV_bright"
    res@cnFillOn=True
    res@lbLabelBarOn=True
    res@cnInfoLabelOn=False
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
    res@cnMinLevelValF       = -0.04               ; min level
    res@cnMaxLevelValF       =  0.08               ; max leve
    res@cnLevelSpacingF      =0.005
    res@lbLabelBarOn=False
    

    res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
    res@tmYLLabelsOn=True
    res@tmYUseLeft=True

    res@tmYLValues    = ispan(2,dimsizes(timeind)/12-1,5)    
    res@tmYLLabels    = ispan(1980,2005,5)
    ; res@tmYLValues    = ispan(1,dimsizes(timeind)/12-1,5)
    ; res@tmYLLabels=ispan(1980, 2010, 5)

    res@tmXTOn=False
    res@tmYROn=False

    res@tmXBValues=(/0,90,180,270/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
    res@tmXBMinorOn=True

    res@vpHeightF=0.25
    res@vpWidthF=0.4
    plot = gsn_csm_contour(wks, p_annual, res)   ; create plot
    

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = 0
    res@cnLineThicknessF = 4.
    res@cnFillOn        = False

    plot1 = gsn_csm_contour(wks, p_annual, res)
    
    overlay(plot, plot1)
  draw(plot)
  frame(wks)
  end 