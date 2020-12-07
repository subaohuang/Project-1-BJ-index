begin

	; 读取 .nc 文件
	tpfile = "tp.assignment3.nc"
	ufile  = "uwnd.assignment3.nc"
	vfile  = "vwnd.assignment3.nc"
	tpf    = addfile(tpfile, "r")
	uf     = addfile(ufile, "r")
	vf 	   = addfile(vfile, "r")
	; print(tpf)
	; print(uf)
	; print(vf)
	

	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	timetp 			= cd_calendar(tpf->time, 2)
	; printVarSummary(timetp)
	timeuf 			= cd_calendar(uf->time, 2)
	timevf 			= cd_calendar(vf->time, 2)

	; 2018年6月27日
	date_tp 		= ind(timetp .eq. 20180627)
	; print(date_tp)
	date_uf 		= ind(timeuf .eq. 20180627)
	date_vf 		= ind(timevf .eq. 20180627)

	; 2018年6月15日至7月15日
	month_tp = ind(timetp .ge. 20180615 .and. timetp .le. 20180715)
	; print(month_tp)

	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp2 			= tpf->tp
	tp22 			= tp2
	tp22(:,:,:)		= (/tp2*1000/)
	tp22@units 		= "mm"
	; printVarSummary(tp2)
	; printMinMax(tp2, 0)
	; printVarSummary(tp22)
	; printMinMax(tp22, 0)

	tp_date = tp22(date_tp,:,:)
	; printVarSummary(tp_date)
	; printMinMax(tp_date, 0)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp_m 	= tp22(month_tp,{35},{127.5})
	; printVarSummary(tp_m)
    ; print(tp_m)
    ; month   = timetp(month_tp)
    ; print(month)
	
	
	; 提取 2018年6月27日 700百帕的风场
	u 	= uf->uwnd({level|700},time|date_uf,lat|:,lon|:)
	; printVarSummary(u)
	v 	= vf->vwnd({level|700},time|date_vf,lat|:,lon|:)


	
	; ========================================================================
	; 打开工作空间
	wkstype = "png"
	wksName = "18304015.assignment3.q1"
	wks 	= gsn_open_wks(wkstype, wksName)


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
	res@gsnCenterString = "(a) precipitation(mm) and 700hPa wind field on 2018.06.27"
	

	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF 		= 15.
	res@mpMaxLatF 		= 45.
	res@mpMinLonF 		= 100.
	res@mpMaxLonF 		= 140.
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion = "HighRes"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette = "precip2_17lev"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLevels 			 = ispan(0, 120, 6)
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
    ; res@tmXBLabelFontHeightF    = 0.024
    ; res@tmYLLabelFontHeightF    = 0.024
    ; res@tiMainFontHeightF       = 0.024
    ; res@gsnCenterStringFontHeightF  = 0.024
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
    ; res@tiMainOffsetYF = -0.02
	res@tmXTOn 		= False
	res@tmYROn 		= False
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 	= "ExplicitLevels"
	rescn@cnLevels 				= (/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLinesOn 			= True
	rescn@cnLineColor 			= "Red"
	rescn@cnLineThicknessF 		= 4
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn 		= False
	rescn@cnInfoLabelOn 		= False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle 			= "LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF 	= 2
	resvc@vcMinDistanceF 			= 0.02

	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF 		= 0.1
	resvc@vcRefMagnitudeF 	= 25
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF        = 1.15
    resvc@vcRefAnnoOrthogonalPosF      = -0.222 
    ; resvc@vcRefAnnoString2On         = False
    
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF         = 0.01
	
	
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
	resxy@gsnCenterString 	= "(b) 6.15-7.15 precipitation time series of 35~S~o~N~N,127.5~S~o~N~E in 2018"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString 	= "time"
	resxy@tiYAxisString 	= "precipitation(mm)"
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 	= True
	
	
	; X 轴
	; tmXBMode 及相应属性
	resxy@tmXBMode 			        = "Explicit"
    resxy@tmXBValues                = ispan(1, 30, 5)
    resxy@tmXBLabels                = (/"6.16","6.21","6.26","7.1","7.6","7.11",""/)
    resxy@tmXBMinorOn               = True
    resxy@tmXBMinorValues           = ispan(0, 31, 1)
 
	
	
	; (可选) XY线及柱状的颜色
	; resxy@gsnXYBarChartColors       = "Blue"
    ; resxy@gsnXYBarChartFillLineThicknessF = 0.03
    resxy@xyLineThicknesses         = 2.5
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@gsnCenterStringFontHeightF        = 0.022
    ;resxy@tiXAxisFontHeightF       = 0.026
    ;resxy@tiYAxisFontHeightF       = 0.026
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resmark                     = True
    resmark@gsMarkerIndex       = 12
    resmark@gsMarkerColor       = "red"
    resmark@gsMarkerSizeF       = 0.02
    resmark@gsMarkerThicknessF  = 6
    ; resmark@tfPolyDrawOrder     = "PreDraw"
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot1 	= gsn_csm_contour_map_ce(wks, tp_date, res)
	
	; 降水=30 mm
	plot2 	= gsn_csm_contour(wks, tp_date, rescn)
	
	; 700百帕风场
	plot3 	= gsn_csm_vector(wks, u, v, resvc)	
	
	; 添加标记
	plotM = gsn_add_polymarker(wks, plot1, 127.5, 35, resmark)


	overlay(plot1, plot3)
	overlay(plot1, plot2)

	; --------------------------------------------------------------------------
	; 图 (b)
	plot4 	= gsn_csm_y(wks, tp_m, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel 			        = True
    resPanel@gsnPanelMainString = "ChenZiqing(18304015)"
    ; resPanel@

	plot 				  = new(2,graphic)
  	plot(0) 			  = plot1
  	plot(1) 			  = plot4
	
	gsn_panel(wks,plot,(/2,1/),resPanel)
	
	
	; =========================================================================
	; 去掉空白空间
	system("convert -trim " + wksName + ".png " + wksName + ".png")


end