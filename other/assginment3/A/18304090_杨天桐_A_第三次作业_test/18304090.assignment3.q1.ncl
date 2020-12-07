; load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl" ；使用time函数
begin

	; 读取 .nc 文件
	path 								= "/home/ys17-19/lsh/Project/other/assginment3/"
	f1 									= addfile(path+ "tp.assignment3.nc", "r")
	f2 									= addfile(path+ "uwnd.assignment3.nc", "r")
	f3 									= addfile(path+ "vwnd.assignment3.nc", "r")
		
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	timeTP 								= f1 ->time	
	YYYYMMDD_TP							= cd_calendar(timeTP, -2)        ;将时间数据转化为YYMMDD格式

	; 2018年6月27日
	timeTP0627Ind 						= ind(YYYYMMDD_TP .eq. 20180627) ;提取20180627时间位置

	; 2018年6月15日至7月15日
	timeTPStart							= 20180615
	timeTPLast 							= 20180715
	timeTPInd 							= ind(YYYYMMDD_TP .ge. timeTPStart .and. YYYYMMDD_TP .le. timeTPLast) ;提取0615~0715时间段的位置
	timeTPWant 							= timeTP(timeTPInd)              ;0615~0715的时间数据
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	rain 								= f1 ->tp
	rain0627 							= rain(timeTP0627Ind,:,:) * 1000
	rain0627!0 							= "lat"
	rain0627&lat   					    = rain&latitude
	rain0627!1 							= "lon"
	rain0627&lon 						= rain&longitude
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	rainWant 							= rain(timeTPInd,{35},{127.5}) * 1000
	; print(rainWant)
	
	; 提取 2018年6月27日 700百帕的风场
	timeWind 							= f2 ->time
	YYYYMMDD_Wind 						= cd_calendar(timeWind, -2)
	timeWindInd 						= ind(YYYYMMDD_Wind .eq. 20180627)   ;提取风场数据的时间位置

	uwnd 								= f2 ->uwnd
	vwnd 								= f3 ->vwnd
	uwndWant							= uwnd(timeWindInd,{700},:,:)
	vwndWant 							= vwnd(timeWindInd,{700},:,:)

	; ========================================================================
	; 打开工作空间
	wkstype 							= "eps"
	wksName 							= "18304090.assignment3.q1.test"
	wks 								= gsn_open_wks(wkstype, wksName)
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	res2          						= res
    ; res2@gsnMaximize 	= True
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
	res@gsnCenterString					= "Rain(mm) and 700 hPa wind(ms~S~-1~N~) in East Asia (180627)"
	res@gsnCenterStringFontHeightF		= 0.024
	
	; 地图范围
	; mpMinlatF, mpMaxlatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF						= 15
	res@mpMaxLatF						= 45
	res@mpMinLonF 						= 100
	res@mpMaxLonF 						= 140
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	res@mpDataBaseVersion 				= "HighRes"
	
	; 色板
	; cnFillPalette
	res@cnFillPalette 					= "WhiteBlueGreenYellowRed"
	res@pmLabelBarParallelPosF 			= 0.4
    res@pmLabelBarOrthogonalPosF 		= 0.1
	res@lbBoxEndCapStyle 				= "TriangleBothEnds"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode 			= "ManualLevels"
	res@cnMinLevelValF 					= 0
	res@cnMaxLevelValF 					= 110
	res@cnLevelSpacingF 				= 10
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
    res@tmXBLabelFontHeightF            = 0.024
    res@tmYLLabelFontHeightF            = 0.024
    res@tiXAxisFontHeightF              = 0.024
    res@tiYAxisFontHeightF              = 0.024
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	res@tmXTOn                          = False
    res@tmYROn                          = False
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 			= "ExplicitLevels"
	rescn@cnLevels 						= (/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor 					= "purple"
	rescn@cnLineThicknessF 				= 3
	
	
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
	resvc@vcLineArrowThicknessF 		= 1.5
	resvc@vcMinDistanceF 			    = 0.025
	
	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF 					= 0.1
	resvc@vcRefMagnitudeF 				= 15
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF 		= 1.05
	resvc@vcRefAnnoOrthogonalPosF 		= -0.23
	resvc@vcRefAnnoPerimSpaceF 			= 0.1
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	resvc@vcRefAnnoFontHeightF 			= 0.01
	
	
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
	resxy@tiMainString 					= "Rain in 180615-180715 (35~S~o~N~N,127.5~S~o~N~E)"
	resxy@tiMainFontHeightF 			= 0.024
	resxy@tiMainFontThicknessF 			= 2
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString 				= "date"
	resxy@tiYAxisString 				= "rain((mm)"
	 
	; gsnXYBarChart
	resxy@gsnXYBarChart 				= True
	
	; X 轴tmXBMode 及相应属性
 	; restick                             = True
    ; restick@ttmFormat                   = "%N%D"
    ; time_axis_labels(timeTPWant,resxy,restick) ; call the formatting procedure 格式化时间数据函数
    dateMajorWant                       = (/20180616,20180621,20180626,20180701,\
                                          20180706,20180711/)
    dateMajorInd                        = get1Dindex(YYYYMMDD_TP,dateMajorWant)
    dateMajor                           = timeTP(dateMajorInd)

    resxy@tmXBMode                      = "Explicit"
    resxy@tmXBValues                    = dateMajor
    resxy@tmXBLabels                    = (/"0616","0621","0626","0701","0706","0711"/)

    resxy@tmXBLabelsOn                  = True
    resxy@tmXBMinorOn                   = True
    resxy@tmXBMinorValues               = timeTPWant
	
	; (可选) XY线及柱状的颜色
    resxy@gsnXYBarChartColors           = "deepskyblue"
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiXAxisFontHeightF 			= 0.024
	resxy@tiYAxisFontHeightF 			= 0.024
    resxy@tmXBLabelFontHeightF          = 0.024
    resxy@tmYLLabelFontHeightF          = 0.024

    plotb 								= gsn_csm_xy(wks, timeTPWant, rainWant, resxy)
    ;！！！显示副刻度日期！！！
    ;获取原图位置及大小
    getvalues plotb 
        "vpXF"      : vpxf
        "vpYF"      : vpyf
        "vpHeightF" : vpheightf
        "vpWidthF"  : vpwidthf
    end getvalues



    res2@vpXF      = vpxf
    res2@vpYF      = vpyf
    res2@vpHeightF = vpheightf
    res2@vpWidthF  = vpwidthf

    ;创建表示前一个图中次要刻度线位置的值,我们将使用这些值作为第二个绘图中的主要刻度线值。
    labels                          = cd_calendar(timeTPWant, -5)
    dateMinorLabels                 = tostring(labels(:,2))                   ;年月日的日部分
    values                          = cd_calendar(timeTPWant, -5)
    values1                         = values(:,2)
    n                               = values1 - 1           				  ;用于区分主刻度与副刻度 
 
    dateMinorLabels                 = where((n%5),dateMinorLabels,"")		  ;主刻度判断值为0，副刻度判断值为1
    res2@tmXBMode                   = "Explicit"
    res2@tmXBValues                 = timeTPWant
    res2@tmXBLabels                 = dateMinorLabels
    res2@tmXBLabelFontHeightF       = 0.01
    res2@tmXBMajorOutwardLengthF    = 0.0
    res2@tmXBMajorLengthF           = 0.0

    res2@tmXBLabelFontColor         = "red"

    res2@tmYROn                     = False            						; Turn off right tickmarks.
    res2@tmXTOn                     = False            						; Turn off top tickmarks.
    res2@tmYLOn                     = False            						; Turn off left tickmarks.
    res2@tmXBLabelDeltaF            = 0.6              						; Move label away from tickmarks.
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark 							= True
	resMark@gsMarkerIndex 				= 3
	resMark@gsMarkerColor 				= "magenta"
	resMark@gsMarkerSizeF 				= 0.025
	resMark@gsMarkerThicknessF 			= 5
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plota 								= gsn_csm_contour_map_ce(wks, rain0627, res)
	
	; 降水=30 mm
	plot1 								= gsn_csm_contour(wks, rain0627, rescn) 			
	overlay(plota, plot1)
	; 700百帕风场
	plot2 								= gsn_csm_vector(wks, uwndWant, vwndWant, resvc)
	overlay(plota, plot2)
	; 添加标记
	plot3 								= gsn_add_polymarker(wks, plota,127.5 , 35, resMark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
    plotb 								= gsn_csm_xy(wks, timeTPWant, rainWant, resxy)
	blank 								= gsn_blank_plot(wks, res2)
	overlay(plotb, blank)

	; =========================================================================
	; 绘画展板
	;gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	;gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
	resPanel							= True
	resPanel@gsnPanelMainString 		= "18304090 yangtiantong" 
	resPanel@gsnPanelXF 				= 0.1
	resPanel@gsnPanelYF 				= 0.7			
	; gsn_panel(wks, (/plota,plotb/), (/2,1/), resPanel)
    draw(plota)
    frame(wks)
	
	; =========================================================================
	; 去掉空白空间

	delete(wks)
    system("convert -geometry 3000x3000 -density 600x600 -trim " +wksName + ".eps " + wksName + ".png")   
    system("rm " + wksName + ".eps -f") 

end