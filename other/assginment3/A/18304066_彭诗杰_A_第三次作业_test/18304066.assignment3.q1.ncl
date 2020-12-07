begin

	; 读取 .nc 文件
    tp_file = addfile("tp.assignment3.nc","r")
    u_file = addfile("uwnd.assignment3.nc","r")
    v_file = addfile("vwnd.assignment3.nc","r")

	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)




	time_tp = tp_file->time
    timeARR_tp = cd_calendar(time_tp, 2)
    time_u = u_file->time
    timeARR_u = cd_calendar(time_u, 2)
    time_v = v_file->time
    timeARR_v = cd_calendar(time_v, 2)

    ; 2018年6月27日
	
    ind_tp1 = ind(timeARR_tp .eq. 20180627)
    ind_u1 = ind(timeARR_u .eq. 20180627)
    ind_v1 = ind(timeARR_v .eq. 20180627)
    ; print(ind_tp1)
	
	; 2018年6月15日至7月15日
	
	ind_tp2 = ind(timeARR_tp .ge. 20180615 .and. timeARR_tp .le. 20180715)
    ind_u2 = ind(timeARR_u .ge. 20180615 .and. timeARR_u .le. 20180715)
    ind_v2 = ind(timeARR_v .ge. 20180615 .and. timeARR_v .le. 20180715)
    ; print(ind_tp2)


	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp = tp_file->tp
    tp1 = tp(ind_tp1,:,:)
    ; printMinMax(tp1, 1)
    tp1 = tp1*1000
    tp1@units = "mm"

    ; printMinMax(tp1, 1)
    ; printVarSummary(tp1)
	; print(tp1)
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
    tp = tp_file->tp
    ; printVarSummary(tp)
    tp2 = tp(time|ind_tp2,{latitude|35},{longitude|127.5})
    
    tp2 = tp2*1000
    
    tp2@units = "mm"


    date_x = timeARR_tp(ind_tp2)
    test_x = ispan(0, 30, 1)
    ; print(date_x)
	; printVarSummary(tp2)

	
	; 提取 2018年6月27日 700百帕的风场
	u = u_file->uwnd
    u1 = u(time|ind_u1,{level|700},lat|:,lon|:)
    v = v_file->vwnd
    v1 = v(time|ind_v1,{level|700},lat|:,lon|:)
    ; printVarSummary(v1)  
    ; wind = wind_speed(u1, v1)
    ; printMinMax(wind, 0)	
	
	; ========================================================================
	; 打开工作空间
    wksName = "18304066.assignment3.q1"
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
	res@tiMainString = "wind_speed at 700hpa and precipitation"
    res@gsnCenterString = "27th August 2018"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF = 15
    res@mpMaxLatF = 45
    res@mpMinLonF = 100
    res@mpMaxLonF = 140
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion      ="HighRes"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette = "precip3_16lev"
	res@lbBoxEndCapStyle = "TriangleHighEnd"


	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0
    res@cnMaxLevelValF = 50
    res@cnLevelSpacingF = 5



	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tiMainFontHeightF = 0.025
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    res@gsnCenterStringFontHeightF = 0.02


	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode = "ExplicitLevels"
    rescn@cnLevels = (/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "black"
    rescn@cnLineThicknessF = 4.
	
	
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
	resvc@vcLineArrowThicknessF = 2.
    resvc@vcMinDistanceF = 0.01
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF = 0.03
	resvc@vcRefMagnitudeF = 15


	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF = 0.183
    resvc@vcRefAnnoOrthogonalPosF = -0.47


	

	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF




	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.53
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString = "precipitation series"
	resxy@gsnCenterString = "15th Jun to 15th Jul 2018"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString = "Time"
    resxy@tiYAxisString = "Total precipitation(mm)"

	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart = True
	

	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode = "Explicit"
    ; resxy@tmXBMajorLine = (/1,6,11,16,21,26/)
    resxy@tmXBValues = (/1,6,11,16,21,26/)
    resxy@tmXBLabels = (/"0616","0621","0626","0701","0706","0711"/)
    resxy@tmXBLabelFontHeightF = 0.02

    resxy@tmXBMinorValues = test_x
    

    
	


	
	; (可选) XY线及柱状的颜色


    resxy@gsnAboveYRefLineColor = "deepskyblue"





	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF

    resxy@tiMainFontHeightF = 0.03
	resxy@tiYAxisFontHeightF = 0.02

    resxy@tiMainOffsetYF  = 0.002
    resxy@tiYAxisOffsetXF = 0.015
    resxy@tiXAxisOffsetXF = -0.016	
    resxy@gsnCenterStringFontHeightF = 0.02
    resxy@gsnCenterStringOrthogonalPosF = 0.04
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    resmark = True
    resmark@gsMarkerIndex = 12
    resmark@gsMarkerColor = "white"
	resmark@gsMarkerSizeF = 0.024
    resmark@gsMarkerThicknessF=3.
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
    plot = gsn_csm_contour_map_ce(wks, tp1, res)	
	
	; 降水=30 mm

    plotcn = gsn_csm_contour(wks, tp1, rescn)
    	
	overlay(plot, plotcn)
	; 700百帕风场
	plotvc = gsn_csm_vector(wks, u1, v1, resvc)
     
	overlay(plot, plotvc)
	; 添加标记
	plotM = gsn_add_polymarker(wks, plot, 127.5, 35, resmark)
	
	
    
	; --------------------------------------------------------------------------
	; 图 (b)
	
	plotxy = gsn_csm_xy(wks, test_x, tp2, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)

    resPanel = True
    resPanel@gsnPanelMainString = "pengshijie(18304066)"




    plots = new(2,graphic)
    plots(0) = plot
    plots(1) = plotxy
    gsn_panel(wks, plots, (/1,2/), resPanel)



	
	
	; =========================================================================
	; 去掉空白空间
	delete(wks)
    system("convert -geometry 3000x3000 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png") 
    system("rm "+wksName+".eps -f")

end