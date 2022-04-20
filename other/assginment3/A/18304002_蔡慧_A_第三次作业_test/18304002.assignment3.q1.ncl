begin

	; 读取 .nc 文件
	ft 				= addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc","r")
	fU 				= addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc","r")
	fV 				= addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc","r")



	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time 	=	ft->time
	timearr = 	cd_calendar(time, -2)
	timeind1= 	ind(timearr.eq.20180627)

	
	
	; 2018年6月15日至7月15日
	timeind2= 	ind(timearr.ge.20180615.and.timearr.le.20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1 		= 	ft->tp(timeind1,:,:)
	tp11 		=	1000*tp1
	copy_VarMeta(tp1, tp11)
	; delete(units)
	tp11@units 	=	 "mm"

	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2 		= 	ft->tp(timeind2,{35},{127.5})
	tp22 		= 	1000*tp2	
	copy_VarMeta(tp2, tp22)
	; delete(units)
	tp22@units 	=	 "mm"
	
	; 提取 2018年6月27日 700百帕的风场
	u 			=	fU->uwnd(timeind1,{700},:,:)
	v 			=	fV->vwnd(timeind1,{700},:,:)
	
	
	; ========================================================================
	; 打开工作空间
	wksName 		= "18304002.assignment3.q1"
	wks 			= gsn_open_wks("eps",wksName)
	
	
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
	res@tiMainString 					= "20180627 precipitation & 700hPa wind in East Asia"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF 						=	15
	res@mpMaxLatF 						=	45
	res@mpMaxLonF 						= 	145
	res@mpMinLonF 						= 	100
	
	; 高分辨率地图
	; mpDataBaseVersion
	; res@mpDataBaseVersion     		 	= 	 "HighRes"
	res@mpDataBaseVersion     		 	= 	 "LowRes"

	res@mpDataResolution				=	 "FinestResolution"
	res@mpDataSetName   				=	 "Earth..4"
	res@mpOutlineBoundarySets     		= 	 "National" 
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette 					=	"BlAqGrYeOrRe"
	res@mpGeophysicalLineThicknessF 	=	 0.8
	res@mpGeophysicalLineColor 			=	"Gray40"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode 			=	"ManualLevels"
	res@cnMaxLevelValF 					= 	 4
	res@cnMinLevelValF 					= 	 0
	res@cnLevelSpacingF 				=    1	
	; res@cnLevels						=	 30
	res@cnLineThicknessF				= 	0.3	 
	
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tiXAxisFontAspectF 	 			= 0.015
	res@tiYAxisFontAspectF 	 			= 0.015
	res@tiMainFontHeightF 		 		= 0.017
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF 					= 0.005
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 			=	"ManualLevels"
	rescn@cnMaxLevelValF				=	 60
	rescn@cnMinLevelValF				=	 0
	rescn@cnLevelSpacingF				=	 15
	
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor 					= 	"Red"
	rescn@cnLineThicknessF 				= 	0.1
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn 				= 	False
	rescn@cnInfoLabelOn   				= 	False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle 					= 	"LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF 	 				= 	0.1
	resvc@vcRefMagnitudeF 				= 	20
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF 		= 	0.55
	resvc@vcRefAnnoOrthogonalPosF 		=	0	
	
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
	resxy@tiMainString 	 				=	"time series(20180615-20180715) of precipitation at 35~S~o~N~N,127.5~S~o~N~E"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString 				= 	"time series/day"
	resxy@tiYAxisString 				= 	"the amount of precipitation/mm"
	resxy@tiXAxisOffsetYF 				= 	-0.02
	resxy@tiYAxisOffsetXF 				= 	0.015
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 				= 	True
	
	
	; X 轴
	; tmXBMode 及相应属性
	date 								= 	ispan(1,31,1)
	resxy@tmXBMode 						= 	"Explicit"
	resxy@tmXBMinorOn	   				=  	True	
	resxy@tmXBValues 			 		=	ispan(2,27,5)
	resxy@tmXBMinorValues 				= 	date
	; resxy@tmXBDataLeftF 				= 	0
	; resxy@tmXBDataRightF 				= 	26
	resxy@tmXBLabelsOn 					= 	True
	resxy@tmXBLabels 					= 	(/"0616","0621","0626","0701","0706","0711"/)	
	; resxy@tmXBTickStartF		 		= 	1	
	; resxy@tmXBTickEndF 	 				=	26
	; resxy@tmXBTickSpacingF 				=	5
	
	
	; (可选) XY线及柱状的颜色
	resxy@xyLineColor 					=	"Blue"

	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiXAxisFontHeightF 			= 0.015
	resxy@tiYAxisFontHeightF 			= 0.015
	resxy@tiMainFontHeightF 			= 0.017
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resmarker 									=	True
	resmarker@gsMarkerIndex 	 				=	13
	resmarker@gsMarkerColor 	 				=	"white"
	resmarker@gsMarkerSizeF 	 				= 	0.015
	; resmarker@gsMarkerThicknessF 				= 	0.1
	; =========================================================================
	; 画图
	plot 								= 	new(2,graphic)
	; 图 (a)
	
	; 降水地图
	plot(0) 							=	gsn_csm_contour_map_ce(wks, tp11, res)
	
	; 降水=30 mm
	plot2 								= 	gsn_csm_contour(wks, tp11, rescn)
	
	; 700百帕风场
	plot3 								= 	gsn_csm_vector(wks, u, v, resvc)
	
	; 添加标记
	plotmarker 							= 	gsn_add_polymarker(wks, plot3, 127.5, 35, resmarker)
	
	overlay(plot(0), plot2)
	overlay(plot(0), plot3)

	; --------------------------------------------------------------------------
	; 图 (b)
	plot(1) 							=	gsn_csm_xy(wks,date , tp22, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	respanel 							= 	True
	respanel@gsnPanelMainString 		=	"CAI HUI(18304002)"
	respanel@tiMainFontHeightF 			= 	0.03
	respanel@gsnPanelXF 				= 	0.1
	respanel@gsnPanelYF 				=	0.02

	gsn_panel(wks, plot, (/1,2/), respanel)
	
	; draw(plot)
	; frame(wks)
	; =========================================================================
	; 去掉空白空间
	delete(wks)
	system("convert -geometry 1200x1200 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
	system("/bin/rm -f " + wksName + ".eps")
	

end