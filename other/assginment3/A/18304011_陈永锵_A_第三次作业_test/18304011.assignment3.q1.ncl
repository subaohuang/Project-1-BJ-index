begin

	; 读取 .nc 文件
    file_path  = "/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc"
    file_path2 = "/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc"
    file_path3 = "/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc"
    file1      = addfile(file_path, "r")
    file2      = addfile(file_path2, "r")
    file3      = addfile(file_path3, "r")
    tp         = file1->tp ; Total precipitation
    time1      = file1->time
    uwnd       = file2->uwnd
    vwnd       = file3->vwnd
    time2      = file2->time
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
    time_tp   = cd_calendar(time1, 2)
    time_wind = cd_calendar(time2, 2)

	; 2018年6月27日
	ind_tp_20180627    = ind(time_tp.eq.20180627)
    ind_wind_20180627  = ind(time_wind.eq.20180627) 
		
	; 2018年6月15日至7月15日
    ind_tp_20180615_0715   = ind(time_tp.ge.20180615.and.time_tp.le.20180715)
    ind_wind_20180615_0715 = ind(time_wind.ge.20180615.and.time_wind.le.20180715)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
    tp_20180627       = tp(ind_tp_20180627,:,:)
    tp_20180627       = (/tp_20180627*1000./)	
    tp_20180627@units = "mm"

	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
    tp_20180615_0715       = tp(ind_tp_20180615_0715,{35},{127.5})
    tp_20180615_0715       = (/tp_20180615_0715*1000/)
    tp_20180615_0715@units = "mm"	
	
	; 提取 2018年6月27日 700百帕的风场
    u           = uwnd(ind_wind_20180627,{700},:,:)
    v           = vwnd(ind_wind_20180627,{700},:,:)
    wind_700hPa = wind_speed(u, v) 

	; ========================================================================
	; 打开工作空间
    type    = "eps"
    wksName = "18304011.assignment3.q1.test"
    wks     = gsn_open_wks(type, wksName)	
	
	; =========================================================================
	; 定义绘图属性
	
    res                = True

    res@gsnDraw        = False
    res@gsnFrame       = False
    res@gsnLeftString  = ""
    res@gsnRightString = ""

    res@tmYROn         = False
    res@tmXTOn         = False 

    resxy              = res

    res@gsnAddCyclic   = False

    rescn              = res
    resvc              = res

    res@vpWidthF       = 0.8
    res@cnFillOn       = True		
    res@cnLinesOn      = False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@tiMainString = "(a)East Asian Precipitation(mm) and 700-hPa Wind Field(ms~S~-1~N~), June 27, 2018 "
		
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF = 15.
    res@mpMaxLatF = 45.
    res@mpMinLonF = 100.
    res@mpMaxLonF = 145.
		
	; 高分辨率地图
	; mpDataBaseVersion
    res@mpDataBaseVersion           = "HighRes"
    res@mpDataResolution            = "Coarse"
    res@mpGeophysicalLineThicknessF = 0.8
    res@mpGeophysicalLineColor      = "#3E2723"	
	; 色板
	; cnFillPalette
	res@cnFillPalette = "precip2_17lev";"precip3_16lev"
		
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMaxLevelValF       = 100.
    res@cnMinLevelValF       = 0.
    res@cnLevelSpacingF      = 5.
    
    ; tp_test = tp(ind_tp_20180627,{15:45},{100:145})
    ; tp_test = (/tp_test*1000./)
    ; printMinMax(tp_test, False) 
    ; result : min=-1.66533e-13   max=161.196
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
    res@tiMainFont               = "times-bold"
    res@tiMainFontHeightF        = 0.023
    res@pmLabelBarParallelPosF   = 0.32
    res@pmLabelBarOrthogonalPosF = 0.12
    res@tmXBLabelFontHeightF     = 0.02
    res@tmYLLabelFontHeightF     = 0.02

	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
    res@tiMainOffsetYF = 0.015

	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@cnFillOn             = False
    rescn@cnLinesOn            = True
    rescn@cnLevelSelectionMode = "ExplicitLevels"
    rescn@cnLevels             = (/30./)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
    rescn@cnLineColor      = "#ffff80";"#292421"
    rescn@cnLineThicknessF = 3.
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
    rescn@cnLineLabelsOn             = True
    rescn@cnInfoLabelOn              = False
    rescn@cnLineLabelFontHeightF     = 0.015
    rescn@cnLineLabelFont            = "helvetica-bold"
    rescn@cnLineLabelBackgroundColor = -1
    rescn@cnLineLabelPlacementMode   = "Constant" 	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle = "LineArrow"	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
    resvc@vcLineArrowThicknessF = 2.
    resvc@vcLineArrowColor      = "#082E54"
    resvc@vcMinDistanceF        = 0.025	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
    resvc@vcRefLengthF    = 0.05
    resvc@vcRefMagnitudeF = 25.
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
    resvc@vcRefAnnoOrthogonalPosF = -0.23
    resvc@vcRefAnnoParallelPosF   = 1.083
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF = 0.015
		
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
    resxy@vpWidthF  = 0.7
    resxy@vpHeightF = 0.3    
    resxy@tmXTOn    = False
    resxy@tmYROn    = False    
    resxy@trYMinF   = 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString = "(b)Precipitation at 35~S~o~N~N, 127.5~S~o~N~E  from June 15 to July 15, 2018"	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString = "Time series"
    resxy@tiYAxisString = "precipitation(mm)"	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart = True	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode             = "Explicit"
    resxy@tmXBValues           = (/1,6,11,16,21,26/)
    resxy@tmXBLabels           = (/"June 16", "June 21", "June 26","July 1", "July 6","July 11"/)
    resxy@tmXBLabelFont        = "helvetica"
    resxy@tmXBLabelFontHeightF = 0.022
    resxy@tmYLLabelFontHeightF = 0.022
    resxy@tmXBMinorValues      = ispan(1, 31, 1)
	
	; (可选) XY线及柱状的颜色
    resxy@gsnYRefLine               = 0.
    resxy@gsnAboveYRefLineBarColors = "Blue4"
    resxy@gsnXYBarChartBarWidth     = 0.8	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
    resxy@tiMainFontHeightF  = 0.024
    resxy@tiMainFont         = "times-bold"
    resxy@tiMainOffsetYF     = 0.013
    resxy@tiYAxisFont        = "times-bold"
    resxy@tiXAxisFont        = "times-bold"
    resxy@tiXAxisFontHeightF = 0.024
    resxy@tiYAxisFontHeightF = 0.024
    resxy@tiXAxisFontColor   = "#292421"
    resxy@tiYAxisFontColor   = "#292421"
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark                    = True
    resMark@gsMarkerIndex      = 5
    resMark@gsMarkerThicknessF = 7.
    resMark@gsMarkerColor      = "#0000FF"
    resMark@gsMarkerSizeF      = 0.022	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	;降水地图
	plot = gsn_csm_contour_map(wks, tp_20180627, res)
	
	; 降水=30 mm
	overlay(plot,gsn_csm_contour(wks, tp_20180627, rescn))

	; 700百帕风场
    overlay(plot, gsn_csm_vector(wks, u, v, resvc))  
    
	; 添加标记
    plotMark = gsn_add_polymarker(wks,plot,127.5,35,resMark)
    
	; --------------------------------------------------------------------------
	; 图 (b)
	plot2 = gsn_csm_xy(wks, ispan(1, 31, 1), tp_20180615_0715, resxy)

	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
    plot3    = new(2,graphic)
    plot3(0) = plot
    plot3(1) = plot2

    resPanel                            = True
    resPanel@gsnPanelYWhiteSpacePercent = 5.
    resPanel@gsnPanelMainString         = "18304011 ChenYongqiang"
    resPanel@gsnPanelMainFont           = "times-bold"
    resPanel@gsnPanelMainFontHeightF    = 0.02
    resPanel@gsnPanelMainPosYF          = 0.97 
	; resPanel@gsnPanelDebug = True

    gsn_panel(wks, plot3, (/2,1/), resPanel)	
	
	; =========================================================================
	; 去掉空白空间
	system("convert -geometry 4000x4000 -density 1200x1200 -trim " + wksName + ".eps " + wksName + ".png")

end