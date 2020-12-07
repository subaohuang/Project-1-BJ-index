begin

	ftp = addfile("tp.assignment3.nc", "r")
    fu  = addfile("uwnd.assignment3.nc", "r")
	fv  = addfile("vwnd.assignment3.nc", "r")
	
	tp   = ftp->tp
    uwnd = fu->uwnd
    vwnd = fv->vwnd
    
    timeflp = ftp->time
    timef1  = cd_calendar(timeflp, 2)
    timef1_ind1 = ind(timef1.eq.20180627)           ;2018年6月27日降水时间地址
    timef1_ind2 = ind((20180615.le.timef1).and.(timef1.le.20180715))  ;6月15日-7月15日降水时间地址
    
    tp1   = tp(timef1_ind1,:,:)       
    tp1   = (/tp1*1000/)                            ;6月27日降水(单位mm)
    tp1@units = "mm"                                

    time  = timef1(timef1_ind2)
    x = ispan(1, 31, 1)
    tp2   = tp(time|timef1_ind2,{latitude|35},{longitude|127.5})                     
    tp2   = (/tp2*1000/)                            ;6月15日-7月15日降水(单位mm)
    tp2@units = "mm"

    timewind = cd_calendar(fu->time, 2)
    timewind_ind = ind(timewind.eq.20180627)        ;2018年6月27日风的时间地址

    u = uwnd(time|timewind_ind,{level|700},lat|:,lon|:) ;2018年6月27日u风
    v = vwnd(time|timewind_ind,{level|700},lat|:,lon|:) ;2018年6月27日v风

    ; printVarSummary(tp1)

    ; printMinMax(tp1,1)   

    ; ========================================================================
    pltName = "18304070.assignment3.q1"
	wks = gsn_open_wks("eps", pltName)
	
	; =========================================================================

	res 								= True

	res@gsnDraw 						= False
	res@gsnFrame 						= False

	res@gsnLeftString 					= ""
	res@gsnRightString 					= ""
	
	res@gsnAddCyclic					= False

	; res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	res@tiMainString = "                                             Picture(a)~C~Total precipitation(shading;mm) and Wind vector(vector:m/s)"
    res@gsnCenterString = "Time:2018/06/27 ; Location:East Asia ~C~  "
	
	res@mpMinLatF = 15
    res@mpMaxLatF = 45
    res@mpMinLonF = 100
    res@mpMaxLonF = 140

    res@tmYLMode   = "Explicit"
    res@tmYLValues = (/15,20,25,30,35,40,45/)
    res@tmYLLabels = (/15,20,25,30,35,40,45/) + "N"

    res@tmYROn = False
    res@tmXTOn = False
	
	res@mpDataBaseVersion = "HighRes"
	
	res@cnFillPalette = "WhiteBlueGreenYellowRed"
 
    ; res@cnLevelSelectionMode = "ManualLevels"
    ; res@cnMaxLevelValF       = 55
    ; res@cnMinLevelValF       = 0
    ; res@cnLevelSpacingF      =  5.
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/0,5,10,15,20,25,30,35,40,45,50,85,125/)
    
    res@tiMainFont = 22 ;设定字体可以让字体粗细一致

	res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    res@tiMainFontHeightF    = 0.025
    res@tiMainFontThicknessF = 1
    res@gsnCenterStringFontHeightF = 0.0225

    res@tiMainOffsetYF = 0.015
    res@gsnCenterStringOrthogonalPosF = 0.031

    res@vpXF = 0.1
    res@vpYF = 0.9
    res@vpHeightF = 0.7
    res@vpWidthF = 0.85

    res@pmLabelBarOrthogonalPosF = 0.1
    res@pmLabelBarParallelPosF   = 0.375
    ; ---------------------------------------------------------------------------
	rescn = True

    rescn@gsnDraw  = False
    rescn@gsnFrame = False

    rescn@gsnAddCyclic      = False

    rescn@gsnLeftString                   = " "
    rescn@gsnRightString                  = " "

	rescn@cnLevelSelectionMode = "ExplicitLevels"
    rescn@cnLevels = 30.
	
	rescn@cnLineColor = "black"
    rescn@cnLineThicknessF = 3.5
	
	rescn@cnLineLabelsOn = False
    rescn@cnInfoLabelOn  = False
	
	; ---------------------------------------------------------------------------
	resvc = True

    resvc@gsnDraw  = False
    resvc@gsnFrame = False

    resvc@gsnAddCyclic      = False
	
	resvc@vcGlyphStyle = "CurlyVector"
	
	resvc@vcLineArrowThicknessF = 1.5
    resvc@vcMinDistanceF        = 0.025
	
	resvc@vcRefLengthF    = 0.1
    resvc@vcRefMagnitudeF = 20
	
	; 参考矢量信息的X,Y方向的位置
	
	resvc@vcRefAnnoParallelPosF   = 1.05
    resvc@vcRefAnnoOrthogonalPosF = -0.225
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

    ; resvc@vcRefAnnoFontHeightF = 
	
	
	; ---------------------------------------------------------------------------
	resxy = True

    resxy@gsnDraw  = False
    resxy@gsnFrame = False

	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	resxy@gsnXYBarChart   = True

    resxy@tiMainString    = "               Picture(b)~C~Time series of precipitation"
    resxy@gsnCenterString = "Time:2018/6/15-2018/7/15 ; Location:35~S~o~N~N,127.5~S~o~N~E"
	
	resxy@tiXAxisString = "Time(day)"
    resxy@tiYAxisString = "Precipitation(mm)"

    resxy@tmXBMode   = "Explicit"
    resxy@tmXBValues = (/2,7,12,17,22,27/)
    resxy@tmXBLabels = (/"16 Jun.","21 Jun.","26 Jun.","1 Jul.","6 Jul.","11 Jul."/)
	resxy@tmXBMinorValues = x
 
	resxy@gsnXYBarChartColors = "skyblue"

    resxy@tiMainFont = 22

	resxy@tiMainFontHeightF = 0.025
    resxy@tiMainFontThicknessF = 1
    resxy@tiXAxisFontHeightF = 0.02
    resxy@tiYAxisFontHeightF = 0.02
    resxy@gsnCenterStringFontHeightF = 0.0225

    resxy@tiMainOffsetYF = 0.015
    resxy@gsnCenterStringOrthogonalPosF = 0.031
    ; =========================================================================
	resmark = True
    resmark@gsMarkerIndex = 14
    resmark@gsMarkerSizeF = 0.015
    resmark@gsMarkerThicknessF = 3.5
	resmark@gsMarkerColor = "White"
	
	; =========================================================================
	
	; 图 (a)
	
	; 降水地图
	plottp = gsn_csm_contour_map(wks, tp1 , res)
	; 降水=30 mm
	plot30mm = gsn_csm_contour(wks, tp1, rescn)
	; 700百帕风场
	plotvc = gsn_csm_vector(wks, u, v, resvc)
	; 添加标记
    lon0 = (/127.5/)
    lat0 = (/35.0/)
	dum = gsn_add_polymarker(wks, plottp, lon0 , lat0, resmark)
	
    overlay(plottp, plotvc)
    overlay(plottp, plot30mm)


    ; draw(plottp)
    ; frame(wks)
	; --------------------------------------------------------------------------
	; 图 (b)
    resxy@vpXF = 0.1
    resxy@vpYF = 0.9
    resxy@vpHeightF = 0.7
    resxy@vpWidthF = 0.9
	plotxy = gsn_csm_xy(wks, x, tp2, resxy)
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	plot = new(2,graphic)
	plot(0) = plottp
    plot(1) = plotxy

    respanel = True
    respanel@gsnPanelXWhiteSpacePercent = 5

    ; drawNDCGrid(wks)

    text1 = (/"Ren longhao"/)
    text2 = (/"18304070"/)
    res1 = True
    res1@txFontHeightF = 0.02
    res2 = True
    res2@txFontHeightF = 0.02
    gsn_text_ndc(wks, text2, 0.49, 0.82, res1)
    gsn_text_ndc(wks, text1, 0.49, 0.77, res2)
	

    gsn_panel(wks, plot, (/1,2/), respanel)

	; =========================================================================
	; 去掉空白空间
	
    delete (wks)            
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    system("rm " + pltName + ".eps -f") 
end