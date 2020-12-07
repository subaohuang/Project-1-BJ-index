begin

	; 读取 .nc 文件

	a = addfile("tp.assignment3.nc","r")
	b = addfile("uwnd.assignment3.nc","r")
	c = addfile("vwnd.assignment3.nc","r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日

	timetp 		= a->time
	timea	    = cd_calendar(timetp,-2)
	timeIndtp0627 	= ind(timea.eq.20180627)
	print(timeIndtp0627)

	timeuv      = b->time
	timeb	    = cd_calendar(timeuv,-2)
	timeInduv0627 	= ind(timeb.eq.20180627)
	print(timeInduv0627)

	; 2018年6月15日至7月15日

	timeIndtp0615 	= ind(timea.eq.20180615)
	timeIndtp0715 	= ind(timea.eq.20180715)
	print(timeIndtp0615)
	print(timeIndtp0715)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp0 = a->tp(time|514,latitude|:,longitude|:)
	tp0 = (/tp0*1000/)
	tp0@units             = "mm"
	tp0&longitude@units   = "degrees_east"
    tp0&latitude@units    = "degrees_north"
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm

	tp1 = a->tp(time|502:532,{latitude|35},{longitude|127.5})
	tp1 = (/tp1*1000/)
	tp1@units             = "mm"

	; 提取 2018年6月27日 700百帕的风场
	u = b->uwnd(time|1734,level|1,lat|:,lon|:)
	v = c->vwnd(time|1734,level|1,lat|:,lon|:)	
	u&lon@units     = "degrees_east"
    u&lat@units     = "degrees_north"
    v&lon@units     = "degrees_east"
    v&lat@units     = "degrees_north"
	
	
	; ========================================================================
	; 打开工作空间

	pltName ="18304105.assignment3.q1"
	wks 	= gsn_open_wks("eps",pltName)

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

	res@gsnCenterString					="Precipitation in Jun 27 and Wind in 700hPa"
	res@gsnCenterStringFontHeightF		= 0.028
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
	res@mpMinLatF						=15.
	res@mpMaxLatF						=45.
	res@mpMinLonF						=100.
	res@mpMaxLonF						=140.
	res@mpCenterLonF       = 120.
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion     		 	= "HighRes"
	;res@mpDataResolution				= "FinestResolution"
	;res@mpDataSetName   				= "Earth..4"
	;res@mpOutlineBoundarySets     		= "National"  
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette				= "WhiteBlueGreenYellowRed"
	res@gsnSpreadColorEnd			= 220	

	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
	res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels             = ispan(-10, 130, 10)
	
	res@lbBoxEndCapStyle        = "TriangleBothEnds"
    res@lbOrientation           = "Vertical"
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	
	;res@tmXBLabelFontHeightF 			= 0.02
	;res@tmYLLabelFontHeightF 			= 0.02
	
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
	rescn@cnLineColor				="Red"
	rescn@cnLineDashPattern      	= 0
	rescn@cnLineThicknessF			=3.

	
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn			=False
	rescn@cnInfoLabelOn				=False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle

	resvc@vcGlyphStyle 	 	= "CurlyVector"
 	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF =2.

	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF 	 	= 0.1	
	resvc@vcRefMagnitudeF 	= 25.
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF 	=1.0
	resvc@vcRefAnnoOrthogonalPosF 	=0.02

	resvc@vcRefAnnoPerimOn              = False
    resvc@vcRefAnnoBackgroundColor      = -1
    resvc@vcRefAnnoString2              = "m/s"
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF          =0.02
	
	
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
	
	resxy@tiMainString					="total Precipitation in Jun 15-Jul 15"
	resxy@gsnCenterStringFontHeightF    = 0.03
	;resxy@gsnXYBarChart					= True
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString 				="Day"
	resxy@tiYAxisString 				="Precipitation (mm)"
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 				=True
	
	; X 轴
	; tmXBMode 及相应属性

	resxy@tmXBOn                        = True
	resxy@tmXBMode						="Explicit"
	resxy@tmXBValues                    = (/1,6,11,16,21,26/)
	resxy@tmXBLabels                    = (/"Jun 16","Jun 21","Jun 26","Jul 1","Jul 6","Jul 11"/)
	
	;resxy@tmXBMinorValues               = fspan(0,30,31)
	
	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors  =  "cornflowerblue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF  =0.02
	resxy@tiXAxisFontHeightF =0.02
	resxy@tiYAxisFontHeightF =0.02
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark						= True
	resMark@gsMarkerIndex 		= 5
	resMark@gsMarkerThicknessF 	= 3.
	resMark@gsMarkerColor 		= "White"
	resMark@gsMarkerSizeF 		= 0.024

	; =========================================================================
	; 画图
	plot 						=new(2,graphic)	
	; 图 (a)

	; 降水地图
	plot(0)    = gsn_csm_contour_map_ce(wks,tp0,res)
	
	; 降水=30 mm
	plotcn	   = gsn_csm_contour(wks,tp0,rescn)


	; 700百帕风场
	plotOverlay = gsn_csm_vector(wks,u,v,resvc)


	; 添加标记
	plotM		= gsn_add_polymarker(wks,plot(0),127.5,35,resMark)

	overlay(plot(0),plotOverlay)
	overlay(plot(0),plotcn)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	time 						= ispan(0,30,1)
	plot(1) 					= gsn_csm_xy(wks,time,tp1,resxy)	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)

    resPanel                           = True
    resPanel@gsnPanelMainFontHeightF   = 0.022	
    resPanel@gsnPanelMainString        = "18304105.zhuyixiao" 
    gsn_panel(wks,plot,(/2,1/),resPanel)
	
	; =========================================================================
	; 去掉空白空间
	delete(wks)
	system("convert -geometry 1600x1600 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
	system("/bin/rm " + pltName + ".eps -f")


end