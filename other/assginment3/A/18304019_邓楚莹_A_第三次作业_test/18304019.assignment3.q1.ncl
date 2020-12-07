begin
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

	; 读取 .nc 文件
    f = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc", "r")
    u = addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc", "r")
    v = addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc", "r")
    ;print(f)
    time  = f->time
    timeARR		= cd_calendar(time,-2)
    time2 = u->time
    timeARR2    = cd_calendar(time2,-2)
    

	;print(timeARR)
	 tp1 = f->tp
     tp11 = tp1
     tp11(:,:,:) = (/1000*tp1/)
     tp11@units = "mm"

	;printVarSummary(tp11)
    u1 = u->uwnd
	v1 = v->vwnd	
	
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	timeInd1 	= ind(timeARR.eq.20180627) ;降水的时间
	timeInd11 	= ind(timeARR2.eq.20180627);风场的时间
	
	; 2018年6月15日至7月15日
	timeInd2 	= ind(timeARR.ge.20180615 .and. timeARR.le.20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp2 = tp11(time|timeInd1,latitude|:,longitude|:)
    ;printVarSummary(tp2)
    
    ;printMinMax(tp2, True)
	
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp3 = tp11(time|timeInd2,{latitude|35},{longitude|127.5})
	;print(tp3)

	
	
	; 提取 2018年6月27日 700百帕的风场
	u2 = u1(time|timeInd11,{level|700},lat|:,lon|:)
	v2 = v1(time|timeInd11,{level|700},lat|:,lon|:)
 
	;printVarSummary(u2)
	;printVarSummary(v2)
	
	
	; ========================================================================
	; 打开工作空间
	plottype = "eps"
	pltName  = "18304019.assignment3.q1_Test"
	wks = gsn_open_wks(plottype, pltName)
	
	
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
	resMark                             = res
	resPanel                            = res

	
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@tiMainString      = "18304019(Deng Chuying)"
	res@gsnCenterString = "(a) Total precipitation of east Asia(mm) and Wind field of 700 hPa(m/s) in June 27,2018"
	res@gsnLeftString   = " "
	res@gsnRightString  = " "
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF = 15
	res@mpMaxLatF = 45
	res@mpMinLonF = 100
	res@mpMaxLonF = 140
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion = "HighRes"
	
	
	; 色板
	; cnFillPalette
	;res@cnFillPalette = "WhiteGreen"
    res@cnFillPalette = "WhiteBlueGreenYellowRed"

	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnFillOn  = True
	res@cnLinesOn = False
	;printMinMax(wspd, True)
	res@cnLevelSelectionMode = "ManualLevels"
	res@cnMinLevelValF       = 0
	res@cnMaxLevelValF       = 110
	res@cnLevelSpacingF      = 5
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF	   = 0.015
	res@tmYLLabelFontHeightF	   = 0.015
	res@tiMainFontHeightF 		   = 0.02
	res@gsnCenterStringFontHeightF = 0.018
	res@pmLabelBarOrthogonalPosF   = 0.15
	res@lbBoxEndCapStyle 		   = "TriangleBothEnds"
	res@tmXTOn                     = False
	res@tmYROn                     = False
	
	; (可选) 主标题Y方向移动位置
	res@tiMainOffsetYF = 0
	
	;图像位置
	res@vpXF		= 0.1
    res@vpYF		= 0.86
    res@vpWidthF	= 0.8
    res@vpHeightF	= 0.3
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnFillOn  = False
	rescn@cnLinesOn = True
	rescn@cnLevelSelectionMode = "ExplicitLevels"
	rescn@cnLevels             = 30
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "black"
	rescn@cnLineThicknessF = 3.0
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn = False
	rescn@cnInfoLabelOn  = False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle = "LineArrow"
	resvc@mpMinLatF = 15
	resvc@mpMaxLatF = 45
	resvc@mpMinLonF = 100
	resvc@mpMaxLonF = 140
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF = 0.5
	resvc@vcLineArrowColor  = "black"
    resvc@vcMinDistanceF 	= 0.025
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF = 0.1
	resvc@vcRefMagnitudeF = 20
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoOrthogonalPosF = -0.25
	resvc@vcRefAnnoParallelPosF   = 1.2
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.6
	resxy@vpHeightF 					= 0.3
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@gsnCenterString = "(b) Total precipitation from June 15 to July 15, 2018 at 35~S~o~N~N,127.5~S~o~N~E"
	resxy@gsnLeftString   = " "
	resxy@gsnRightString  = " "
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString = "Date"
	resxy@tiYAxisString = "Total precipitation(mm)"
	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart = True
	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode   = "Explicit"
    resxy@tmXBValues = (/1,6,11,16,21,26/)
    resxy@tmXBLabels = (/"June 16","June 21","June 26","July 1","July 6","July 11"/)
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = ispan(0,31,1)
	
	
	; (可选) XY线及柱状的颜色
	;resxy@gsnXYBarChartColors   = "Blue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	resxy@tiYAxisFontHeightF       = 0.015
	resxy@tiXAxisFontHeightF       = 0.015
	resxy@gsnCenterStringFontHeightF = 0.02
	
	;图像位置
	resxy@vpXF		= 0.2
    resxy@vpYF		= 0.4
    
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	x		= 127.5
    y		= 35
	resMark@gsMarkerIndex = 12
	resMark@gsMarkerColor = "black"
	resMark@gsMarkerThicknessF = 2.0
	resMark@gsMarkerSizeF = 0.015
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plotBase1       = gsn_csm_contour_map_ce(wks,tp2,res)
	

	; 降水=30 mm，黑色粗线
	plotOverlay1    = gsn_csm_contour(wks,tp2,rescn)
	overlay(plotBase1, plotOverlay1)
	
	; 700百帕风场
	plotOverlay2    = gsn_csm_vector(wks,u2,v2,resvc)
	overlay(plotBase1, plotOverlay2)
	; 添加标记，星星标志
	plotMarker = gsn_add_polymarker(wks, plotBase1, x, y, resMark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plotBase2      = gsn_csm_y(wks, tp3, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel 							= True
	resPanel@gsnPanelMainString 		= "Deng Chuying(18304019)"
	resPanel@gsnPanelMainFontHeightF 	= 0.02
	resPanel@gsnPanelLabelBar 			= True
	resPanel@lbLabelFontHeightF 		= 0.015
	resPanel@lbBoxEndCapStyle           ="TriangleBothEnds"
    resPanel@gsnPanelYWhiteSpacePercent = 1.0

	gsn_panel(wks,(/plotBase1,plotBase2/),(/2,1/),resPanel)
	
	draw(wks)
	frame(wks)
	
	
	
	; =========================================================================
	; 去掉空白空间
    delete (wks) 			
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    system("rm " + pltName + ".eps -f") 
end