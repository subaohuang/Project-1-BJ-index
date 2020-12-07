begin

	; 读取 .nc 文件

    path        = "C:\Users\Halo\Desktop\assignment3\"
    f_tp        = addfile(path+"tp.assignment3.nc", "r")
    f_uwnd      = addfile(path+"uwnd.assignment3.nc", "r")
    f_vwnd      = addfile(path+"vwnd.assignment3.nc", "r")

   ; print(f_tp)
    ;print(f_uwnd)
	;print(f_vwnd)

    ;printVarSummary(f_uwnd)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time1 =  f_tp->time;仅仅提取时间维
	timeARR1=cd_calendar(time1, 2)
	time2=f_uwnd->time
	timeARR2=cd_calendar(time2, 2)
    time1_ind=ind(timeARR1.eq.20180627)
    time2_ind=ind(timeARR2.eq.20180627)
	; 2018年6月15日至7月15日
	
	time3_ind=ind(timeARR1.ge.20180615.and.timeARR1.le.20180715)
	print(time1_ind)
	print(time2_ind)
	print(time3_ind)
		; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	tp1      = f_tp ->tp(time|514,latitude|:,longitude|:);提取某时间的数据
	tp1       =1000*(tp1)
	;tp1@units  ="mm"；;确定单位
	printMinMax(tp1,0)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp3      = f_tp ->tp(time|502:532,{latitude|35},{longitude|127.5});提取某时间的数据
	tp3       =1000*(tp3)
	;tp3@units  ="mm"；;确定单位
	printMinMax(tp3,0)

	
	
	; 提取 2018年6月27日 700百帕的风场

	uwnd1=f_uwnd->uwnd(time|1734,{level|700},lat|:,lon|:);
    vwnd1=f_vwnd->vwnd(time|1734,{level|700},lat|:,lon|:);

	
	printMinMax(uwnd1, 0)
	printMinMax(vwnd1, 0)
	
	; ========================================================================
	; 打开工作空间
	
	pltName="18304055.assignment3.q1"
    wks = gsn_open_wks("eps", pltName);使用eps裁剪
	
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
	res@tiMainString="2018.06.27 Precipitation and Wind of 700hpa of East Asia"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF   	=15
	res@mpMaxLatF       =45
	res@mpMinLonF       =100
	res@mpMaxLonF       =140
	
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"
	res@mpDataResolution="Fine"
	
	; 色板
	; cnFillPalette
	res@cnFillPalette           = "MPL_Blues"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 		= "ExplicitLevels";等值线间距取值方式
    res@cnLevels = ispan(0,170,5)

	res@lbBoxEndCapStyle = "RectangleEnds "

	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF = 0.02
	res@tmYLLabelFontHeightF = 0.02
	
	res@tiMainFontHeightF=0.02
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF = 0.025
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode    = "ExplicitLevels"
	rescn@cnLevels                = 30
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "blue"
	rescn@cnLineThicknessF = 4
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnLineLabelsOn=False
	rescn@cnInfoLabelOn=False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle     = "CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowColor  = "black"
	resvc@vcLineArrowThickness=1
 	resvc@vcMinDistanceF    = 0.025

	;resvc@vcMonoLineArrowColor 		= False
    ;resvc@vcLevelPalette	 				= "wgne15";颜色代表矢量大小

	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF=0.05; 参考矢量的长度
    resvc@vcRefMagnitudeF=10; 参考矢量的长度(NDC单位)表示的数值


	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
    resvc@vcRefAnnoParallelPosF   = 1.2
	resvc@vcRefAnnoOrthogonalPosF = -0.2
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF = 0.01
	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	 resxy@vpWidthF 				    = 0.7
	 resxy@vpHeightF 					= 0.4
	
	 resxy@tmXTOn 						= False
	 resxy@tmYROn 						= False
	
	 resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题

	
	resxy@tiMainString="Precipitation of 20180615 to 20180715 Precipitation in East Asia"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString = "Time(day)"
	resxy@tiYAxisString = "Precipitation(mm)"
	
	
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart = True
	
	; X 轴
	; tmXBMode 及相应属性

	resxy@tmXBMode        = "Explicit"
    resxy@tmXBValues      = (/503,508,513,518,523,528/)
    resxy@tmXBLabels      = (/"06.16","06.21","06.26","07.01","07.06","07.11"/)
	resxy@tmXBMinorValues =ispan(502,532,1)
	
	; (可选) XY线及柱状的颜色
	

	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF  = 0.02
	resxy@tiXAxisFontHeightF = 0.02
	resxy@tiYAxisFontHeightF = 0.02
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
	resMark				= True
	resMark@gsMarkerIndex 		= 1
	resMark@gsMarkerThicknessF 	= 3.
	resMark@gsMarkerColor 		= "red"
	resMark@gsMarkerSizeF 		= 0.025
	resText				= True
	resText@txFontHeightF 		= 0.025
	resText@txFontColor 		= "red"
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plotBase    = gsn_csm_contour_map_ce(wks,tp1, res)
	
	; 降水=30 mm
	
	plotcn         = gsn_csm_contour(wks,tp1,rescn)
	
	; 700百帕风场
	plotOverlay = gsn_csm_vector(wks, uwnd1,vwnd1, resvc)
	


    overlay(plotBase,plotcn)
    overlay(plotBase,plotOverlay)
	; 添加标记
	plotMark        =gsn_add_polymarker(wks, plotBase, 127.5, 35.0, resMark)
	plotText		= gsn_add_text(wks,plotBase,"35~S~o~N~N,127.5~S~o~N~E",133.5,33.0,resText)

	

	; --------------------------------------------------------------------------
	; 图 (b)
	plotBase1 = gsn_csm_xy(wks, ispan(502,532,1), tp3, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
    resPanel                            = True
    resPanel@gsnPanelYWhiteSpacePercent = 5
    resPanel@gsnPanelMainString         = "18304055 liang zhiliang "
    resPanel@gsnPanelMainFontHeightF    = 0.02

 
    resPanel@lbLabelFontHeightF         = 0.02
 
	
	gsn_panel(wks,(/plotBase, plotBase1/), (/2,1/), resPanel)
	; =========================================================================
	; 去掉空白空间
    delete (wks)            
    system("convert -geometry 1600x1600 -density 1000x1000 -trim " + pltName + ".eps " + pltName + ".png")
    system("rm " + pltName + ".eps -f")
	

end