begin

	; 读取 .nc 文件

    path = "D:\课程\ncl\assignment3\"

    f_pre = addfile(path + "tp.assignment3.nc", "r")
    f_u = addfile(path + "uwnd.assignment3.nc", "r")
    f_v = addfile(path + "vwnd.assignment3.nc", "r")

    time_pre = f_pre->time
    time_u = f_u->time
    time_v = f_v->time

    pre_a = f_pre->tp
    uwnd_a = f_u->uwnd
    vwnd_a = f_v->vwnd

    ; printVarSummary(pre_a)
    ; printVarSummary(uwnd_a)
    ; printVarSummary(vwnd_a)

    ; printVarSummary(f_pre)
    ; printVarSummary(f_v)
    ; printVarSummary(f_u)	
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time_convert_pre = cd_calendar(time_pre, -2)
    time_convert_u = cd_calendar(time_u, -2)
    time_convert_v = cd_calendar(time_v, -2)

    ; print(time_convert_v)
    ; print(time_convert_pre)
	time_convert_pre_627Ind = ind(time_convert_pre .eq. 20180627)
    ; print(time_convert_pre_627Ind)
    time_convert_u_627Ind = ind(time_convert_u .eq. 20180627)
    ; print(time_convert_u_627Ind)
    time_convert_v_627Ind = ind(time_convert_v .eq. 20180627)



	; 2018年6月15日至7月15日
	time_convert_pre_615715Ind = ind(time_convert_pre .ge. 20180615 .and. time_convert_pre .le. 20180715)
    time_convert_u_615715Ind = ind(time_convert_u .ge. 20180615 .and. time_convert_v .le. 20180715)
    time_convert_v_615715Ind = ind(time_convert_v .ge. 20180615 .and. time_convert_v .le. 20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	pre_627 = pre_a(time_convert_pre_627Ind, :, :)
    ; printMinMax(pre_627, 0)
    pre_627 = 1000*pre_627
    ; printVarSummary(pre_627)
    printMinMax(pre_627, 0)
	
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	pre_615715 = pre_a(time_convert_pre_615715Ind, {35}, {127.5})
    pre_615715 = 1000*pre_615715

    ; printVarSummary(pre_615715)
    ; printMinMax(pre_615715, 0)
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	uwnd_627 = uwnd_a(time_convert_u_627Ind, {700}, :, :)
    vwnd_627 = vwnd_a(time_convert_v_627Ind, {700}, :, :)

    wind_627 = wind_speed(uwnd_627, vwnd_627)

    ; printVarSummary(wind_627)
    ; printMinMax(wind_627, 0)
	
	
	
	; ========================================================================
	; 打开工作空间
	pltName = "18304038.assignment3.q1"
	wks = gsn_open_wks("eps", pltName)
	
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
    ; resMark                             = res
				
	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	res@cnLineLabelsOn                  = False
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图

    ; res@cnFillOn = True
    ; res@cnLinesOn = False
	
	; 标题
    res@tiMainString = "(a)"
	res@gsnCenterString = "2018 0627 Precipitation(mm) and 700hpa wind vector(ms~S~-1~N~) of East-Asia"

	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF = 15.
    res@mpMaxLatF = 45.
    res@mpMinLonF = 100.
    res@mpMaxLonF = 140.
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion = "HighRes"
    res@mpDataSetName="Earth..4"
    res@mpDataResolution = "Coarse"
    res@mpGeophysicalLineThicknessF = 0.4
	
	; 色板
	; cnFillPalette

	res@lbBoxEndCapStyle = "TriangleBothEnds"
    ; res@lbLabelBarOn = False ;去掉色标
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    cmap = read_colormap_file("cmocean_balance")
    res@cnFillPalette = cmap(135:0, :)

	res@cnLevelSelectionMode = "ManualLevels"
    res@cnMaxLevelValF = 170.
    res@cnMinLevelValF = -10.
    res@cnLevelSpacingF = 2.

	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
    res@tmXTOn = False
    res@tmYROn = False

	res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@tiMainFontHeightF = 0.024
    res@gsnCenterStringFontHeightF = 0.02

	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
    res@tiMainOffsetYF = 0.

    ; plot = gsn_csm_contour_map_ce(wks, pre_627, res)
    ; draw(plot)
    ; frame(wks)
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性

	rescn@cnLevelSelectionMode = "ExplicitLevels"
    rescn@cnLevels = 30.

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "Green4"
    rescn@cnLineThicknessF = 2.5
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
    rescn@cnLineLabelsOn = False
    rescn@cnInfoLabelOn = False
    rescn@cnFillOn = False
	
	; plot = gsn_csm_contour_map_ce(wks, pre_627, res)
 ;    plotcn = gsn_csm_contour(wks,pre_627,rescn)
 ;    overlay(plot, plotcn)

 ;    draw(plot)
 ;    frame(wks)
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle = "LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcMinDistanceF        = 0.025
    resvc@vcLineArrowThicknessF = 1.3
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefLengthF          = 0.08
    resvc@vcRefMagnitudeF       = 25.
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
    resvc@vcRefAnnoOn = True ;去掉矢量图标
	resvc@vcRefAnnoParallelPosF = 1 
    resvc@vcRefAnnoOrthogonalPosF = -0.44
    resvc@vcRefAnnoString2On = False
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF = 0.01

    resvc@vcMonoLineArrowColor        = False
    resvc@vcLevelPalette                  = "cmocean_thermal"
    resvc@lbLabelBarOn = False


    ; resMark@gsMarkerIndex = 6
    ; resMark@gsMarkerThicknessF  = 3.
    ; resMark@gsMarkerColor       = "red"
    ; resMark@gsMarkerSizeF       = 0.015


    ; plot = gsn_csm_contour_map_ce(wks, pre_627, res)
    ; plotM = gsn_add_polymarker(wks, plot, 127.5, 35., resMark)
    ; plotcn = gsn_csm_contour(wks,pre_627,rescn)
    ; plotvc = gsn_csm_vector(wks, uwnd_627, vwnd_627, resvc)
    ; overlay(plot, plotcn)
    ; overlay(plot, plotvc)
    ; draw(plot)
    ; frame(wks)
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.53
	
    resxy@tmXBOn                        = True
    resxy@tmYLOn                        = True
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString = "(b)"
    resxy@gsnCenterString = "2018 0615-0715 35~S~o~N~N 127.5~S~o~N~E Precipitation Time Series"

	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString = "Time Series(YMD)"
    resxy@tiYAxisString = "Precipitation(mm)"
	
	
	; gsnXYBarChart
    resxy@gsnYRefLine           = 0.
	resxy@gsnXYBarChart  = True
    resxy@gsnAboveYRefLineColor = "red"              ; above ref line fill red
    resxy@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue
	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode = "Explicit"
    ; resxy@tmXBValues = time_convert_pre(time_convert_pre_615715Ind)
    resxy@tmXBValues = ispan(1, sizeof(time_convert_pre(time_convert_pre_615715Ind))/4 - 5, 5)
    time_need = new(sizeof(time_convert_pre(time_convert_pre_615715Ind))/4, "integer")
    time_need = time_convert_pre(time_convert_pre_615715Ind)
    printVarSummary(time_need)
    time_labels = time_need((/1, 6, 11, 16, 21, 26/))
    ; resxy@tmXBLabels = time_convert_pre(time_convert_pre_615715Ind)((/1, 6, 11, 16, 21, 26/))
    resxy@tmXBLabels = time_labels
    ; resxy@tmXBLabelDirection = ""
    ; resxy@tmXBLabelFontAspectF = 90.
    resxy@tmXBLabelAngleF = 45.
    ; resxy@
	resxy@tmXBMinorValues = ispan(0, 30, 1)
	; printVarSummary(time_convert_pre(time_convert_pre_615715Ind))
    ; print(sizeof(time_convert_pre(time_convert_pre_615715Ind))/4)
	; (可选) XY线及柱状的颜色
	; resxy@gsnScale = True    
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	resxy@tiMainFontHeightF  = 0.018
    resxy@tiMainFontThicknessF = 2.5
    resxy@tiMainFont = 0
    resxy@tiXAxisFontHeightF = 0.024
    resxy@tiYAxisFontHeightF = 0.024
    resxy@tiMainOffsetYF = 0.01
    resxy@gsnCenterStringOrthogonalPosF = 0.03
    resxy@gsnCenterStringFontHeightF = 0.02

    x = new(sizeof(time_convert_pre(time_convert_pre_615715Ind))/4, "integer")
    ; print(sizeof(time_convert_pre(time_convert_pre_615715Ind)))
    x = ispan(0, sizeof(time_convert_pre(time_convert_pre_615715Ind))/4 - 1, 1)
    ; plot = gsn_csm_xy(wks, x, pre_615715, resxy)

    ; draw(plot)
    ; frame(wks)
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF

	resMark = True
    resMark@gsMarkerIndex = 6
    resMark@gsMarkerThicknessF  = 3.
    resMark@gsMarkerColor       = "red"
    resMark@gsMarkerSizeF       = 0.015


    ; plot = gsn_csm_contour_map_ce(wks, pre_627, res)
    ; plotM = gsn_add_polymarker(wks, plot, 127.5, 35., resMark)
    ; plotcn = gsn_csm_contour(wks,pre_627,rescn)
    ; plotvc = gsn_csm_vector(wks, uwnd_627, vwnd_627, resvc)
    ; overlay(plot, plotcn)
    ; overlay(plot, plotvc)
    ; draw(plot)
    ; frame(wks)
    
	
	; =========================================================================
	; 画图
	
	; 图 (a)

	; 降水地图
    plot = gsn_csm_contour_map_ce(wks, pre_627, res)	
	
	; 降水=30 mm
    plotcn = gsn_csm_contour(wks,pre_627,rescn)
	
	; 700百帕风场
    plotvc = gsn_csm_vector(wks, uwnd_627, vwnd_627, resvc)
	
	; 添加标记
	plotM = gsn_add_polymarker(wks, plot, 127.5, 35., resMark)

    overlay(plot, plotcn)
    overlay(plot, plotvc)
	
	; --------------------------------------------------------------------------
	; 图 (b)
    plotxy = gsn_csm_xy(wks, x, pre_615715, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
    resPanel                            = True
    resPanel@gsnPanelMainString         = "huang yong tao(18304038)" ;图板标题
    resPanel@gsnPanelMainFontHeightF    = 0.02 ;标题字体大小
    resPanel@gsnPanelLabelBar           = False ;图板色标
    ; resPanel@lbLabelFontHeightF         = 0.02 ;色标字体大小
    resPanel@gsnPanelYWhiteSpacePercent = 5. ;图之间相对位置
    resPanel@gsnPanelMainFontColor = "Red" ;标题颜色

    plotBase = new(2, graphic)
    plotBase(0) = plot
    plotBase(1) = plotxy

    gsn_panel(wks, plotBase, (/1,2/), resPanel)
    delete (wks)
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    system("rm " + pltName + ".eps -f")
	
	
	; =========================================================================
	; 去掉空白空间
	

end