begin

	; 读取 .nc 文件

    path        = ""
    f_tp        = addfile(path+"tp.assignment3.nc", "r")
    f_uwnd      = addfile(path+"uwnd.assignment3.nc", "r")
    f_vwnd      = addfile(path+"vwnd.assignment3.nc", "r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	
	tp_time        = f_tp->time
    uwnd_time      = f_uwnd->time
    vwnd_time      = f_vwnd->time

    tp_timeARR      = cd_calendar(tp_time, -2)
    uwnd_timeARR    = cd_calendar(uwnd_time, -2)
    vwnd_timeARR    = cd_calendar(vwnd_time, -2)

    tp_timeInd      = ind(tp_timeARR.eq.20180627)       ;tp_timeInd = 514
    uwnd_timeInd    = ind(uwnd_timeARR.eq.20180627)     ;uwnd_timeInd = 1734
    vwnd_timeInd    = ind(vwnd_timeARR.eq.20180627)     ;vwnd_timeInd = 1734
	
	; 2018年6月15日至7月15日
	
	tp_timeInd_start       = ind(tp_timeARR.eq.20180615)       ;tp_timeInd_start = 502
    tp_timeInd_end         = ind(tp_timeARR.eq.20180715)       ;tp_timeInd_end = 532
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	tp         = f_tp->tp(time|514,{latitude|15:45},{longitude|100:140})
    tp         = tp*1000
    tp         = tp > 0          ;0到161.196
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	tp_start_to_end         = f_tp->tp(time|502:532,{latitude|35.0},{longitude|127.5})
    tp_start_to_end         = tp_start_to_end*1000
    tp_start_to_end         = tp_start_to_end > 0        ;0到88.1043
	
	; 提取 2018年6月27日 700百帕的风场
	
	uwnd       = f_uwnd->uwnd(time|1734,{level|700},{lat|15:45},{lon|100:140})        ;-9.2到25.425
    vwnd       = f_vwnd->vwnd(time|1734,{level|700},{lat|15:45},{lon|100:140})        ;-12.675到10.525
	
	; ========================================================================
	; 打开工作空间
	
    pltName     = "18304040.assignment3.q1"
	wks         = gsn_open_wks("eps", pltName)
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	res@gsnLeftString 					= ""
	res@gsnRightString 					= ""
	
	resxy	 							= res
	
	res@gsnAddCyclic					= False
	rescn								= res
	resvc 								= res
				
	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
    res@cnSmoothingOn                   = True
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	
    res@tiMainFont          = 0
	res@tiMainString        = "a)precipitation and wind field of 700hPa on 2018.06.27 in east Asia"
	res@tiMainFontHeightF   = 0.014
    res@tiMainOffsetYF      = 0.005

	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
    res@mpMinLatF       = 15
    res@mpMaxLatF       = 45
    res@mpMinLonF       = 100
    res@mpMaxLonF       = 140

    ; 去除右边和上边刻度

    res@tmXTOn      = False
    res@tmYROn      = False
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	res@mpDataBaseVersion           = "HighRes"
    res@mpDataResolution            = "Fine"
	
	; 色板
	; cnFillPalette
	
	res@cnFillPalette       = "WhiteBlue"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
	res@cnLevelSelectionMode        = "ExplicitLevels"
    res@cnLevels                    = ispan(0,165,15)
    res@lbTitleOn                   = True
    res@lbTitleString               = "precipitation (mm)"
    res@lbTitlePosition             = "Bottom"
    res@lbTitleFontHeightF          = 0.017
    res@lbTitleDirection            = "Across"
    res@lbTitleJust                 = "CenterCenter"

    res@pmLabelBarWidthF            = 0.55
    res@pmLabelBarParallelPosF      = 0.38
    res@pmLabelBarOrthogonalPosF    = 0.1
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
    rescn@cnLevelSelectionMode      = "ExplicitLevels"
    rescn@cnLevels                  = 30
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	
    rescn@cnLineThicknessF          = 1
    rescn@cnLineColor               = "blue"
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
    rescn@cnLineLabelsOn                = True
    rescn@cnInfoLabelOn                 = False
    rescn@lbLabelBarOn                  = False
    rescn@cnLabelMasking                = True
    rescn@cnLineLabelBackgroundColor    = -1
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	
	resvc@vcGlyphStyle          = "CurlyVector"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	resvc@vcLineArrowThicknessF        = 0.5
    resvc@vcMinDistanceF               = 0.01
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefLengthF         = 0.06
    resvc@vcRefMagnitudeF      = 3.0
    resvc@vcLineArrowColor     = "red"
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	
	resvc@vcRefAnnoOn               = True
    resvc@vcRefAnnoParallelPosF     = 1.0
    resvc@vcRefAnnoOrthogonalPosF   = -0.23
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

    resvc@vcRefAnnoString1            = "3 m/s"
	resvc@vcRefAnnoString2            = "wind speed"
    resvc@vcRefAnnoFontHeightF        = 0.016

	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.6
	resxy@vpHeightF 					= 0.6
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	
	resxy@tiMainString      = "b)time series of precipitation during 2018.06.15 to 2018.07.15 at 35~S~o~N~N,127.5~S~o~N~E"

	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString         = "date"
    resxy@tiYAxisString         = "precipitation (mm)"
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart        = True
	
	; X 轴
	; tmXBMode 及相应属性

    resxy@tmXBMode                  = "Explicit"
    resxy@tmXBValues                = tp_time((/503,508,513,518,523,528/))
    resxy@tmXBLabels                = (/"16 Jun","21 Jun","26 Jun","1 Jul","6 Jul","11 Jul"/)
    resxy@tmXBLabelFontHeightF      = 0.018
    resxy@tmXBMinorValues           = tp_time(502:532)

	; (可选) XY线及柱状的颜色
	
	resxy@gsnXYBarChartColors      = "blue"
    resxy@gsnXYBarChartBarWidth    = 10
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
    resxy@tiMainFont               = 0
	resxy@tiMainFontHeightF        = 0.014
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
    respm                       = True
    respm@gsMarkerColor         = "black"
    respm@gsMarkerIndex         = 5
    respm@gsMarkerSizeF         = 0.015
    respm@gsMarkerThicknessF    = 2
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
    plotbase        = gsn_csm_contour_map_ce(wks, tp, res)
	
	; 降水=30 mm
	plot_explicit   = gsn_csm_contour(wks, tp, rescn)
    overlay(plotbase, plot_explicit)
	
	; 700百帕风场
    plot_wind_field     = gsn_csm_vector(wks, uwnd, vwnd, resvc)
	overlay(plotbase, plot_wind_field)
	
	; 添加标记

    marker      = gsn_add_polymarker(wks, plotbase, 127.5, 35.0, respm)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	
	plot_xy        = gsn_csm_xy(wks, tp_time(502:532), tp_start_to_end, resxy)
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
    resP                            = True
    resP@gsnFrame                   = False
    resP@gsnPanelMainString         = "18304040 hzl"
    resP@gsnPanelMainFontHeightF    = 0.03
    resP@gsnPanelLabelBar           = False
    resP@gsnPanelXWhiteSpacePercent = 8

    gsn_panel(wks, (/plotbase, plot_xy/), (/1, 2/), resP)

    frame(wks)
	
	; =========================================================================
	; 去掉空白空间
	
    delete(wks)
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    system("rm " + pltName + ".eps -f")

end