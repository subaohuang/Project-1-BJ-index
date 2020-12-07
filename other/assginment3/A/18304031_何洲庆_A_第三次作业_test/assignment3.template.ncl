begin

	; 读取 .nc 文件
    path="/home/ys17-19/lsh/Project/other/assginment3/"

    f_1= addfile(path + "tp.assignment3.nc", "r")
    f_u = addfile(path + "uwnd.assignment3.nc", "r")
    f_v = addfile(path + "vwnd.assignment3.nc", "r")

    time_1 = f_1->time
    time_u = f_u->time
    time_v = f_v->time

    f_a = f_1->tp
    u_a = f_u->uwnd
    v_a = f_v->vwnd
	
	; printVarSummary(u_a)
	; printVarSummary(v_a)
	printVarSummary(f_a)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time_11=cd_calendar(time_1,-2)
	time_u1=cd_calendar(time_u,-2)
	time_v1=cd_calendar(time_v,-2)
    
    ;printVarSummary(time_1)

	time_627_1=ind(time_11.eq.20180627)
	time_627_u=ind(time_u1.eq.20180627)
	time_627_v=ind(time_v1.eq.20180627)
    
    printVarSummary(time_627_u)
	;printVarSummary(time_627_1)

	
	; ; 2018年6月15日至7月15日
	; ; time_12=cd_calendar(time_1, -2)
	; ; time_u2=cd_calendar(time_u,-2)
	; ; time_v2=cd_calendar(time_v,-2)

	time_615_1=ind(time_11.ge.20180615.and.time_11.le.20180715)
	time_615_u=ind(time_u1.ge.20180615.and.time_u.le.20180715)
	time_615_v=ind(time_v1.ge.20180615.and.time_v.le.20180715)

	;printVarSummary(time_615_11)


	
	
	; ; ========================================================================
	; ; 提取 2018年6月27日的降水数据
	; ; 单位需要由 m 转为 mm
	data_627_m=f_a(time_627_1,:,:)
	data_627_m=1000*data_627_m
	;printMinMax(data_627_mm, 0)
	;printVarSummary(data_627_m)
	;print("___________")
	
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	data_615715_m=f_a(time_615_1,{35},{127.5})
	data_615715_m=1000*data_615715_m
	
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	;print("___________________")
	;print(time_627_u)
	;printVarSummary(u_a)
	data_627_u=u_a(time_627_u, {700}, :, :)
	data_627_v=v_a(time_627_v, {700}, :, :)

	;printVarSummary(data_627_v)

	
	
	
	
	; ========================================================================
	; 打开工作空间
	pltName="18304031_hezhouqing"
	wks=gsn_open_wks("eps",pltName)
	
	
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
	res@tiMainString                    ="(a)"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF=15.
    res@mpMaxLatF=45.
    res@mpMinLonF=100.
    res@mpMaxLonF=135.
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"
	
	
	; 色板
	; cnFillPalette
    res@cnFillPalette = "BlueYellowRed"

	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode = "ManualLevels"
    res@cnMaxLevelValF = 220.
    res@cnMinLevelValF = -10.
    res@cnLevelSpacingF = 2.

	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm

	rescn@cnLinesOn = True
	rescn@cnFillOn = False
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode ="ExplicitLevels"
	rescn@cnLevels             =30.
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor         ="red"
	rescn@cnLineThicknessF    =1.5
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn       =False
	rescn@cnInfoLabelOn        =False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle         ="LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF = 1.5
	resvc@vcMinDistanceF        = 0.025
    
	
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF          =0.08
	resvc@vcRefMagnitudeF       =15
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF  =1
	resvc@vcRefAnnoOrthogonalPosF=0.02
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF =0.025
	
	
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
	resxy@tiMainString = "(b)"

	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString   ="Time"
    resxy@tiYAxisString   = "Precipitation(mm)"
	
	
	
	
	; gsnXYBarChart

	resxy@gsnXYBarChart         = True
    resxy@gsnAboveYRefLineColor = "green"              
    resxy@gsnBelowYRefLineColor = "red"             
	
	
	; X 轴
	; tmXBMode 及相应属性

	; resxy@tmXBMode ="Explicit"
	; resxy@xyXStyle ="Irregular"
	; resxy@xyXIrregularPoints =(/2018616/2018621/2018626/2018701/2018706/2018711/)
	
	
	resxy@tmXBMode 		= "Explicit"
	resxy@tmXBValues 	= time_1((/503,508,513,518,523,528/))
	resxy@tmXBLabels 	= (/"2018616","2018621","2018626","2018701","2018706","2018711"/)
	resxy@tmXBMinorValues = time_1(502:532:1)
	; (可选) XY线及柱状的颜色
	
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF  = 0.02
    resxy@tiXAxisFontHeightF = 0.024
    resxy@tiYAxisFontHeightF = 0.024
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark = True
    resMark@gsMarkerIndex       = 4
    resMark@gsMarkerThicknessF  = 2.
    resMark@gsMarkerColor       = "black"
    resMark@gsMarkerSizeF       = 0.1


	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	plot     =gsn_csm_contour_map_ce(wks, data_627_m, res)
	; 降水地图
	plotcn   =gsn_csm_contour(wks, data_627_m, rescn)
	
	; ; 降水=30 mm
	
	
	; ; 700百帕风场
    plotvc   =gsn_csm_vector(wks, data_627_u, data_627_v, resvc)
	
	; ; 添加标记
	plotM    = gsn_add_polymarker(wks, plot, 127.5, 35., resMark)

    overlay(plot, plotcn)
    overlay(plot, plotvc)
    
	
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plotxy   =gsn_csm_xy(wks, time_1(502:532), data_615715_m, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	;gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	;gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
	resPanel                            = True
    resPanel@gsnPanelMainString         = "assignment3 hezhouqing 18304031" 
    resPanel@gsnPanelMainFontHeightF    = 0.03
    resPanel@gsnPanelLabelBar           = False 
    resPanel@gsnPanelMainFontColor      = "blue"

    plotBase = new(2, graphic)
    plotBase(0) = plot
    plotBase(1) = plotxy

	
	
	; ; =========================================================================
	; ; 去掉空白空间
	gsn_panel(wks, plotBase, (/2,1/), resPanel)
    
     delete (wks)
     system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
     system("rm " + pltName + ".eps -f")
	

end