begin

	; 读取 .nc 文件
    f1 = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc", "r")
    f2 = addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc", "r")
    f3 = addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc", "r")
    tp = f1->tp
    timetp = f1->time
    lon = f1->longitude
    lat = f1->latitude
    uwnd = f2->uwnd
    timew = f2->time
    vwnd = f3->vwnd

    ;printVarSummary(tp)
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
    timetp := cd_calendar(timetp, 2)
    timew := cd_calendar(timew, 2)
	; 2018年6月27日
    timetp_ind1 = ind(timetp.eq.20180627)
    timew_ind1 = ind(timew.eq.20180627)
	
	
	; 2018年6月15日至7月15日
	timetp_ind2 = ind(timetp.ge.20180615.and.timetp.le.20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
    ; 需要注意运算之后没了单位，要重新给他单位
	tp1 = tp(time|timetp_ind1,latitude|:,longitude|:)*1000.
	tp1!1 = "longitude"
    tp1&longitude = lon
    lon@unit = "degree_east"
    tp1!0 = "latitude"
    tp1&latitude = lat
    lat@unit = "degree_north"

	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2 = tp(time|timetp_ind2,{latitude|35},{longitude|127.5})*1000.
    tp2_time = timetp(timetp_ind2)
    ;print(tp2_time)
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	uwnd1 = uwnd(time|timew_ind1,{level|700},lat|:,lon|:)
    vwnd1 = vwnd(time|timew_ind1,{level|700},lat|:,lon|:)
	
	
	; printMinMax(tp1, 1)
 ;    printMinMax(tp2, 1)
 ;    printMinMax(uwnd1, 1)
 ;    printMinMax(vwnd1, 1)
 ;    printVarSummary(tp1)
 ;    printVarSummary(tp2)
 ;    printVarSummary(uwnd1)
 ;    printVarSummary(vwnd1)
	; ========================================================================
	; 打开工作空间
    wksName = "18304054.assignment3.q1"
	wks = gsn_open_wks("eps", wksName)
	
	
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
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@gsnCenterString = "(a) East Asia precipitation(mm) and 700hPa wind(ms~S~-1~N~) on Jun27,2018"
    ;res@gsnLeftString = "(a)"
    ;res@tiXAxisString = "precipitation (mm)"
    ;res@tiYAxisString = "wind speed (ms~S~-1~N~)"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
    res@mpMinLatF = 15
    res@mpMaxLatF = 45
    res@mpMinLonF = 100
    res@mpMaxLonF = 140
    res@mpGeophysicalLineColor = "Gray40"
    res@mpGeophysicalLineThicknessF = 0.8
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/15,25,35,45/)
    res@tmYLLabels = (/"15~S~o~N~N","25~S~o~N~N","35~S~o~N~N","45~S~o~N~N"/)
    res@tmYLMinorOn = True
    res@tmYLMinorValues = (/20,30,40/)
    res@tmXBMode = "Explicit"
    res@tmXBValues = (/100,110,120,130,140/)
    res@tmXBLabels = (/"100~S~o~N~E","110~S~o~N~E","120~S~o~N~E","130~S~o~N~E","140~S~o~N~E"/)
    res@tmXBMinorValues = (/105,115,125,135/)
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion = "HighRes"
	
	
	; 色板
	; cnFillPalette
    res@cnFillPalette = "CBR_wet"
    res@lbBoxEndCapStyle = "TriangleHighEnd"
    res@pmLabelBarOrthogonalPosF = 0.1
	
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnInfoLabelOn = False
    res@cnLineLabelsOn = False
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMaxLevelValF = 80
    res@cnMinLevelValF = 0
    res@cnLevelSpacingF = 10
	
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
    res@tmXTOn = False
    res@tmYROn = False
    res@tmXBLabelFontHeightF = 0.025
    res@tmYLLabelFontHeightF = 0.025
    res@gsnCenterStringFontHeightF = 0.025
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@cnLevelSelectionMode = "ExplicitLevels"
    rescn@cnLevels = 30
	
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "orange"
    rescn@cnLineThicknessF = 3
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn = False
    rescn@cnInfoLabelOn = False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle = "LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF = 1
    resvc@vcMinDistanceF = 0.025
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF = 0.1
    resvc@vcRefMagnitudeF = 20
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	;resvc@vcRefAnnoParallelPosF = 1.2
    resvc@vcRefAnnoOrthogonalPosF = -0.49
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF = 0.01
	resvc@vcRefAnnoString2On = False
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@gsnCenterString = "(b) Precipitation on Jun15 to Jul15, 2018 (35~S~o~N~N,127.5~S~o~N~E)"
    resxy@gsnCenterStringOrthogonalPosF = 0.04
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString = "Date"
    resxy@tiYAxisString = "Precipitation(mm)"
	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart = True
	
	
	; X 轴
	; tmXBMode 及相应属性
    ; resxy@tmXBMode = "Explicit"
    ; resxy@tmXBValues = tp2_time(1:26:6)
    ; resxy@tmXBLabels = (/"16Jun","21Jun","26Jun","1Jul","6Jul","11Jul"/)

    resxy@tmXBMode = "Explicit"
    resxy@tmXBValues = (/2,7,12,17,22,27/)
    resxy@tmXBLabels = (/"16Jun","21Jun","26Jun","1Jul","6Jul","11Jul"/)
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = (/1,3,4,5,6,8,9,10,11,13,14,15,16,18,19,20,21,23,24,25,26,28,29,30,31/)


	; (可选) XY线及柱状的颜色
    ;resxy@xyDashPattern = 0
    ;resxy@xyMarkerColor = "red"
    colors = read_colormap_file("WhiteBlue")
	resxy@gsnXYBarChartColors = colors(174,:)
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@gsnCenterStringFontHeightF = 0.025
    resxy@tiXAxisFontHeightF = 0.025
    resxy@tiYAxisFontHeightF = 0.025
    resxy@tmXBLabelFontHeightF = 0.02
    resxy@tmYLLabelFontHeightF = 0.02
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resmarker = True
    resmarker@gsMarkerIndex = 5
    resmarker@gsMarkerColor = "red"
    resmarker@gsMarkerSizeF = 15
    resmarker@gsMarkerThicknessF = 3
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plotpre = gsn_csm_contour_map_ce(wks, tp1, res)
	
	; 降水=30 mm
	plot30 = gsn_csm_contour(wks, tp1, rescn)
	
    
	; 700百帕风场
	plotvc = gsn_csm_vector(wks, uwnd1, vwnd1, resvc)

	; 添加标记
    lat1 = 35
    lon1 = 127.5
	plotmarker = gsn_add_polymarker(wks, plotpre, lon1, lat1, resmarker)
    overlay(plotpre, plot30)
    overlay(plotpre, plotvc)

	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	;plotxy = gsn_csm_xy(wks, tp2_time, tp2, resxy)
    x = ispan(1,31,1)
    plotxy = gsn_csm_xy(wks, x, tp2, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
    ; plots = new((/1,2/),graphic)
    ; plots(0,0) = plotpre
    ; plots(0,1) = plotxy
	respanel = True
    respanel@gsnPanelMainString = "18304054(Liang Kaixin)"
    respanel@gsnPanelMainFontHeightF = 0.02
    respanel@gsnPanelYWhiteSpacePercent = 2.

    gsn_panel(wks, (/plotpre,plotxy/), (/2,1/), respanel)
	
	
	; =========================================================================
	; 去掉空白空间
	
    delete (wks)
    system("convert -geometry 3000x3000 -density  600x600 -trim " + wksName + ".eps " + wksName + ".png")
    system("/bin/rm -f " + wksName + ".eps")

end