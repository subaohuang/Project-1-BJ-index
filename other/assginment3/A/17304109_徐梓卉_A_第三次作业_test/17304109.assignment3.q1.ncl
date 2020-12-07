begin

	; 读取 .nc 文件
	tpf  = addfile("tp.assignment3.nc", "r")
	uf   = addfile("uwnd.assignment3.nc", "r")
	vf   = addfile("vwnd.assignment3.nc", "r")
	;printFileVarSummary(vf, "vwnd")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	time        = tpf->time
	timeConvert = cd_calendar(time, 2)
	; 2018年6月27日
	timeInd1    = ind(timeConvert.eq.20180627)
	
	; 2018年6月15日至7月15日
	timeInd2    = ind(timeConvert.ge.20180615 .and. timeConvert.le.20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tpf1    = tpf->tp(time|timeInd1,{latitude|15:45},{longitude|100:140})
	;printVarSummary(tpf1)
	tp1         = 1000*tpf1
	copy_VarCoords(tpf1, tp1)
	;printMinMax(tp1, 0)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tpf2    = tpf->tp(time|timeInd2,{latitude|35},{longitude|127.5})
	tp2         = 1000*tpf2
	copy_VarCoords(tpf2, tp2)
	time2 	= ispan(1, 31, 1)
	;print(tp2)
	; 提取 2018年6月27日 700百帕的风场
	uwnd 	= uf->uwnd(time|timeInd1,{level|700},{lat|15:45},{lon|100:140})
	vwnd 	= vf->vwnd(time|timeInd1,{level|700},{lat|15:45},{lon|100:140})
	
	
	; ========================================================================
	; 打开工作空间
	pltName = "17304109.assignment3.q1"
    wks 	= gsn_open_wks("eps", pltName)
	
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
	
	res@tmXTOn 							= False
	res@tmYROn 							= False
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@tiMainString = "Precipitation(mm) and wind field of 700 hpa over East Asia, June 27,2018"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF 					= 15
	res@mpMaxLatF 	 				= 45
	res@mpMinLonF 	 				= 100
	res@mpMaxLonF	 				= 140
	res@mpGeophysicalLineColor 		= "gray50"
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion 			= "HighRes"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette 				= "wind_17lev"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode 		= "ManualLevels"
	res@cnMaxLevelValF				= 80.
	res@cnMinLevelValF				= 0.
	res@cnLevelSpacingF				= 5.
	res@cnLineThicknessF			= 2.
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF 		= 0.02
	res@tmYLLabelFontHeightF 		= 0.02
	res@tiMainFontHeightF    		= 0.025
	;res@gsnCenterStringFontHeightF =
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF  			= 0.00
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 		= "ExplicitLevels"
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor 				= "Black"
	rescn@cnLevels 					= 30
	rescn@cnLineThicknessF          = 4
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnInfoLabelOn 	       	= False
	rescn@cnLineLabelsOn 	     	= False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle 				= "CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF 	= 1.5
	resvc@vcMinDistanceF 			= 0.015
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF 				= 0.1
	resvc@vcRefMagnitudeF 			= 20
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoOn 				= True
	resvc@vcRefAnnoParallelPosF 	= 0.188
	resvc@vcRefAnnoOrthogonalPosF 	= -1.36
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	;resvc@vcRefAnnoFontHeightF 		= 0.02
	;resvc@vcMonoLineArrowColor 		= True
	resvc@vcLineArrowColor          = "Gray20"
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 					= 0.7
	resxy@vpHeightF 				= 0.4
	
	resxy@tmXTOn 					= False
	resxy@tmYROn 					= False
	
	resxy@trYMinF 					= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString 				= "35~S~o~N~N 127.5~S~o~N~E, Jun 15 - Jul 15, 2018"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString 			= "Time"
	resxy@tiYAxisString 			= "Precipitation"
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 			= True
	
	; X 轴
	; tmXBMode 及相应属性
	resxy@tmXBMode 					= "Explicit"
	resxy@tmXBValues 	 			= (/2,7,12,17,22,27/)
	resxy@tmXBLabels 				= (/"0616","0621","0628","0701","0706","0711"/)
	resxy@tmXBMinorValues 			= time2


	; (可选) XY线及柱状的颜色
	colors 							= read_colormap_file("wind_17lev")
	index 							= (/15/)
	resxy@gsnAboveYRefLineColor 	= colors(index,:)
	resxy@xyMarkLineMode 			= "Lines"
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF 		= 0.025
	resxy@tiXAxisFontHeightF 		= 0.02
	resxy@tiYAxisFontHeightF 		= 0.02
	resxy@tiMainOffsetYF  			= 0.00
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resmark 						= True
	resmark@gsMarkerIndex 			= 6
	resmark@gsMarkerColor 			= "Black"
	resmark@gsMarkerSizeF 			= 0.015
	resmark@gsMarkerThicknessF 		= 6
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot   = gsn_csm_contour_map(wks, tp1, res)
	
	; 降水=30 mm
	plot1  = gsn_csm_contour(wks, tp1, rescn)
	
	; 700百帕风场
	plot2  = gsn_csm_vector(wks, uwnd, vwnd, resvc)
	
	
	overlay(plot, plot1)
	overlay(plot, plot2)
	; 添加标记
	plotM  = gsn_add_polymarker(wks, plot, 127.5, 35, resmark)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plot_b = gsn_csm_xy(wks, time2, tp2, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel = True
	resPanel@gsnPanelMainString = "17304109 xuzihui"
	resPanel@gsnPanelXF 		= 0
	resPanel@gsnPanelYF 		= 0
	;resPanel@gsnPanelXF
	;resPanel@gsnPanelYF
	plotf 	 = new(2,graphic)
	plotf(0) = (plot)
	plotf(1) = (plot_b)
	gsn_panel(wks, plotf, (/2,1/), resPanel)
	
	
	
	; =========================================================================
	; 去掉空白空间
	frame(wks)
	delete (wks)           ; 必须先关闭 workstation
    ;system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    ;system("rm " + pltName + ".eps -f")

end