begin

	; 读取 .nc 文件
    f1       = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc","r")
    f2       = addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc","r")
    f3       = addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc","r")
    time1 = f1->time
    time2 = f2->time
    time3 = f3->time

	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
    timeARR1 = cd_calendar(time1,-2)
    timeARR2 = cd_calendar(time2,-2)
    timeARR3 = cd_calendar(time3,-2)

	; 2018年6月27日
	timeInd11   = ind(timeARR1.eq.20180627)
	timeInd12   = ind(timeARR2.eq.20180627)
	timeInd13   = ind(timeARR3.eq.20180627)
	
	; 2018年6月15日至7月15日
	timeInd21   = ind(timeARR1.ge.20180615.and.timeARR1.le.20180715)
    timeInd22   = ind(timeARR2.ge.20180615.and.timeARR2.le.20180715)
    timeInd23   = ind(timeARR3.ge.20180615.and.timeARR3.le.20180715)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
    tp1 = f1->tp
    tp11 = tp1
    tp11 = (/1000*tp1/)
    tp11@units = "mm"

	tp111 = tp11(time|timeInd11,latitude|:,longitude|:)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp222 = tp11(time|timeInd21,{latitude|35},{longitude|127.5})
	
	; 提取 2018年6月27日 700百帕的风场
	u1 = f2->uwnd(time|timeInd12,{level|700},lat|:,lon|:)
	v1 = f3->vwnd(time|timeInd13,{level|700},lat|:,lon|:) 
	
	; ========================================================================
	; 打开工作空间
	pltName = "18304088.assignment3.q1"
	wks 	= gsn_open_wks("pdf",pltName)
	
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
				
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@gsnCenterString = "(a) Total Precipitation in East Asia (mm) and 700 hPa Wind Field (ms~S~-1~N~) on June 27, 2018"
    res@gsnLeftString = ""
    res@gsnRightString = ""

	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF = 15.
	res@mpMaxLatF = 45.
	res@mpMinLonF = 100.
	res@mpMaxLonF = 140.
	res@mpCenterLonF       = 120.

	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion = "HighRes"
	
	; 色板
	; cnFillPalette
	res@cnFillPalette = "WhViBlGrYeOrRe"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnFillOn = True
    res@cnLinesOn = False

	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels = ispan(0,110,5)
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    res@tiMainFontHeightF = 0.05
    res@tiMainFontThicknessF = 2.
	res@gsnCenterStringFontHeightF = 0.02
	res@lbBoxEndCapStyle = "TriangleBothEnds"

	res@lbOrientation           = "Vertical"
	; res@lbOrientation = "Horizontal"
 ;    res@tmXTOn = False
	; res@tmYROn = False

	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
    res@tiMainOffsetYF = 0

    res@vpWidthF 						= 0.5
	res@vpHeightF 					= 0.3
	res@mpShapeMode = "FreeAspect"
	
		
	; res 								= True
	; res@gsnDraw 						= False
	; res@gsnFrame 						= False
	; res@gsnLeftString 					= ""
	; res@gsnRightString 					= ""
	
	; resxy	 							= res
	
	; res@gsnAddCyclic					= False
	; rescn								= res
	; resvc 								= res
				
	; res@vpWidthF 						= 0.7
				
	; res@cnFillOn 						= True		
	; res@cnLinesOn						= False
	
	; ; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; ; res: 图 (a) 底图
	
	; ; 标题
	; res@gsnCenterString                 = "precipitation in Jun 27 and wind in 700hPa"
 ;    res@gsnCenterStringFontHeightF      = 0.028
	
	; ; 地图范围
	; ; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
	; res@mpMinLatF          = 15.
 ;    res@mpMaxLatF          = 45.
 ;    res@mpMinLonF          = 100.
 ;    res@mpMaxLonF          = 140.
 ;    res@mpCenterLonF       = 120.
	
	; ; 高分辨率地图
	; ; mpDataBaseVersion
	
	; res@mpDataBaseVersion   = "HighRes"
    
	
	; ; 色板
	; ; cnFillPalette
	
	; res@cnFillPalette       = "precip3_16lev"
	
	; ; 等值线间距
	; ; cnLevelSelectionMode 及相应属性
	; res@cnLevelSelectionMode    = "ExplicitLevels"
 ;    res@cnLevels                = ispan(0, 200, 10)
 ;    res@lbBoxEndCapStyle        = "TriangleBothEnds"
 ;    res@lbOrientation           = "Vertical"
	
	
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
	rescn@cnLevelSelectionMode   = "ExplicitLevels"
    rescn@cnLevels               = 30.
    

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor            = "darkorchid2"
    rescn@cnLineDashPattern      = 0
    rescn@cnLineThicknessF       = 2.5
    
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn         = False
    rescn@cnInfoLabelOn          = False
    
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@cnFillOn = False
    rescn@cnLinesOn = True
	rescn@cnLevelSelectionMode = "ExplicitLevels"
	rescn@cnLevels = 30
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "Red"
	rescn@cnLineThicknessF = 3.5
	
	
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
	resvc@vcLineArrowThicknessF = 0.6
    resvc@vcMinDistanceF = 0.02
	resvc@vcLineArrowColor = "Black"

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF = 0.1
	resvc@vcRefMagnitudeF = 20.
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF = 1.15
	resvc@vcRefAnnoOrthogonalPosF = -0.2
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	
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
	resxy@gsnCenterString = "(b) Total precipitation from June 15 to July 15, 2018 at 35~S~o~N~N,127.5~S~o~N~E (mm)"
    resxy@gsnLeftString   = " "
    resxy@gsnRightString  = " "

	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString = "Time(day)"
	resxy@tiYAxisString = "Tp(mm)"
	
	;gsnXYBarChart 
	resxy@gsnXYBarChart = True
	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode   = "Explicit"
    resxy@tmXBValues = ispan(1,26,5)
    resxy@tmXBLabels = (/"June 16","June 21","June 26","July 1","July 6","July 11"/)
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = ispan(0,31,1)
	
	
	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors = "Red"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@gsnCenterStringFontHeightF = 0.02
    resxy@tiYAxisFontHeightF = 0.02
    resxy@tiXAxisFontHeightF = 0.02
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark = True
	x1 = 127.5
	y1 = 35
	resMark@gsMarkerIndex = 7
	resMark@gsMarkerColor = "Yellow"
	resMark@gsMarkerSizeF = 0.03
    resMark@gsMarkerThicknessF = 3.0

	; =========================================================================
	; 画图
	; plot = new(2,graphic)

	; 图 (a)
	
	; 降水地图
	plot1 = gsn_csm_contour_map_ce(wks,tp111,res)
	
	; 降水=30 mm
	plotcn = gsn_csm_contour(wks,tp111,rescn)
	overlay(plot1,plotcn)
	
	; 700百帕风场
	plotwind = gsn_csm_vector(wks,u1,v1,resvc)
	overlay(plot1,plotwind)

	; 添加标记
	plotaMarker = gsn_add_polymarker(wks, plot1, x1, y1, resMark)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plot2 = gsn_csm_y(wks, tp222, resxy)
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel = True
    resPanel@gsnPanelMainString = "Xu Xionghui(18304088)"
    resPanel@lbBoxEndCapStyle ="TriangleBothEnds"

    resPanel@gsnPanelXF = (/0.29,0.29/)
    resPanel@gsnPanelYF = (/0.8,0.35/)

	gsn_panel(wks, (/plot1,plot2/), (/2,1/), resPanel)

	; =========================================================================
	; 去掉空白空间
    delete (wks) 			
    ; system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    ; system("rm " + pltName + ".eps -f") 

end