 begin

	; 读取 .nc 文件
	f_uwnd 				= addfile("uwnd.assignment3.nc", "r")
	f_vwnd 				= addfile("vwnd.assignment3.nc", "r")
	f_tp 				= addfile("tp.assignment3.nc"  , "r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time_tp 			= f_tp->time
	timeARR_tp			= cd_calendar(time_tp,-2)
	timeInd_tp1 		= ind(timeARR_tp.eq.20180627)
	
	time_wind   		= f_uwnd->time
	timeARR_wind 		= cd_calendar(time_wind,-2)
	timeInd_wind 		= ind(timeARR_wind.eq.20180627)
	
	
	; 2018年6月15日至7月15日
	timeInd_tp2 		= ind(timeARR_tp.ge.20180615 .and. timeARR_tp.le. 20180715) 

	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1 				= f_tp->tp(timeInd_tp1,:,:)	
	tp1					= (/tp1*1000./)
	printMinMax(tp1, 0)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2 				=f_tp->tp(timeInd_tp2,{35},{127.5}) 
	tp2					= (/tp2*1000./)
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	u 					= f_uwnd->uwnd(timeInd_wind,{700},:,:)
	v 					= f_vwnd->vwnd(timeInd_wind,{700},:,:)
	wind 				= wind_speed(u, v)
	printMinMax(wind, 0)
	
	
	; ========================================================================
	; 打开工作空间
	wksName 			= "18304020.assignment3.q1"
	wks  				= gsn_open_wks("eps", wksName)
	
	
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
	res@tiMainString 					= "(a)20180627 Precipitation in East Asia (unit: mm) and 700 hPa wind (unit: ms~S~-1~N~)"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF 						= 15
	res@mpMaxLatF 						= 45
	res@mpMinLonF 						= 100
	res@mpMaxLonF 						= 140
	
	res@pmTickMarkDisplayMode			= "Always"
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataResolution 				= "CoarseResolution"
	res@mpDataBaseVersion 				= "HighRes"
	res@mpDataSetName 					= "Earth..4"
	res@mpOutlineOn 			 		= True
	res@mpOutlineBoundarySets     		= "National" 
	res@mpFillOn 						= False
	res@mpGeophysicalLineColor			= "Gray20"
	res@mpGeophysicalLineThicknessF 	= 0.8
	
	
	; 色板
	; cnFillPalette
	res@cnFillOn 						= True
	res@cnFillPalette 					= "precip2_17lev"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode 			= "ManualLevels"
	; res@cnLevels 						= ispan(0,100,10)
	res@cnMaxLevelValF					= 100
	res@cnMinLevelValF 					= 0
	res@cnLevelSpacingF 				= 10
	
	
	res@cnLinesOn 						= False
	res@cnInfoLabelOn 					= False
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF 			= 0.018
	res@tmYLLabelFontHeightF 			= 0.018
	res@tiMainFontHeightF 				= 0.02
	res@tiMainFont 						= "times-bold"
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF 					= -0.01
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 			= "ExplicitLevels"
	rescn@cnLevels 						= 30
	
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLinesOn 					= True
	rescn@cnLineColor 					= "Black"
	rescn@cnLineThicknessF 				= 3.

	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn 				= False
	rescn@cnInfoLabelOn 				= False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle 					= "LineArrow"	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcMinDistanceF 				= 0.025
	resvc@vcLineArrowThicknessF 		= 1.5
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF	 				= 0.08
	resvc@vcRefMagnitudeF 				= 25
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoSide        			= "Right"
	resvc@vcRefAnnoOrthogonalPosF 		= -0.01
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.55
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString 					= "(b)20180615-20180715 precipitation (35~S~o~N~N, 127.5~S~o~N~E)"

	
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString 				= "time"
	resxy@tiYAxisString 				= "precipitation(mm)"
	
	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 				= True
	
	
	; X 轴
	; tmXBMode 及相应属性
	time_x 								= time_tp(timeInd_tp2)
	resxy@tmXBMode 						= "Explicit"
	resxy@tmXBValues 					= time_x(::5)
	resxy@tmXBLabels 					= (/"0615","0620","0625","0630","0705","0710","0715"/)
	resxy@tmXBMinorOn 					= True
	resxy@tmXBMinorValues 				= time_x
	 		
	
	
	; (可选) XY线及柱状的颜色
	colors 	 							= read_colormap_file("MPL_Blues")
	resxy@gsnXYBarChartColors 			= colors(12,:)
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF 			= 0.02
	resxy@tiMainFont 					= "times-bold"
	resxy@tiMainOffsetYF 				= 0.02
	
	resxy@tiXAxisFontHeightF 		= 0.024
	resxy@tiYAxisFontHeightF 			= 0.024
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark						= True
	resMark@gsMarkerIndex 		= 15
	resMark@gsMarkerThicknessF 	= 2.
	resMark@gsMarkerColor 		= "Red"
	resMark@gsMarkerSizeF 		= 0.02
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot 						= gsn_csm_contour_map(wks, tp1, res)
	
	
	; 降水=30 mm                                                                                                                                                                                                                         
	plotcn 						= gsn_csm_contour(wks, tp1, rescn)
	overlay(plot, plotcn)
	; 700百帕风场
	plotvc 						= gsn_csm_vector(wks, u, v, resvc)
	overlay(plot, plotvc)
	
	; 添加标记
	
	
	plotM						= gsn_add_polymarker(wks,plot,127.5,35.,resMark)

	
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plotxy 						= gsn_csm_xy(wks, time_x, tp2, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	; drawNDCGrid(wks)
	resPanel 							= True
	resPanel@gsnPanelMainString 		= "dengxiaopeng(18304020)"
	resPanel@gsnPanelMainFontHeightF 	= 0.02
	resPanel@gsnPanelMainPosYF			= 0.97
	resPanel@gsnPanelDebug 				= True
	gsn_panel(wks, (/plot,plotxy/), (/2,1/), resPanel)
	
	
	
	; =========================================================================
	; 去掉空白空间
	delete(wks)
	system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
	system("/bin/rm -f " + wksName + ".eps")

end