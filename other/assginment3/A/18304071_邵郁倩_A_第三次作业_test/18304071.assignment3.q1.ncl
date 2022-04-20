begin

	; 读取 .nc 文件

    f1      =addfile("tp.assignment3.nc", "r")
	f2      =addfile("uwnd.assignment3.nc", "r")
	f3      =addfile("vwnd.assignment3.nc", "r")

    tp   =f1->tp
    uwnd =f2->uwnd
    vwnd =f3->vwnd

    ; printVarSummary(tp)
    ; printVarSummary(uwnd)
    ; printVarSummary(vwnd)

	; ========================================================================
	
    ; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
    time1 =f1->time
    time1_1 =cd_calendar(time1, -2)  

    time2 =f2->time
    time2_1 =cd_calendar(time2, -2)

    time3 =f3->time
    time3_1 =cd_calendar(time3, -2)

	; 2018年6月27日
	time1_1_1 = ind(time1_1.eq.20180627)
    time2_1_1 = ind(time2_1.eq.20180627)
    time3_1_1 = ind(time3_1.eq.20180627)	

	; 2018年6月15日至7月15日
	time1_1_2 =ind(time1_1.ge.20180615.and.time1_1.le.20180715)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1 =tp(time|time1_1_1,latitude|:,longitude|:)
	tp1_1 =tp1*1000
    copy_VarMeta(tp1, tp1_1)
    tp1_1@unit ="mm"
    ; printVarSummary(tp1_1)	

	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2 =tp(time|time1_1_2,{latitude|35},{longitude|127.5})
    tp2_1 =tp2*1000
    copy_VarMeta(tp2, tp2_1)
    tp2_1@unit ="mm"
	; printVarSummary(tp2_1)

	; 提取 2018年6月27日 700百帕的风场
	uwnd1 =uwnd(time|time2_1_1,{level|700},lat|:,lon|:)
	vwnd1 =vwnd(time|time3_1_1,{level|700},lat|:,lon|:)
	; printVarSummary(uwnd1)
	; printVarSummary(vwnd1)

	; ========================================================================
	; 打开工作空间
	pltName ="18304071.assignment3.q1"
	wks =gsn_open_wks("eps",pltName)
	
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
	resvc                               = res

	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图

	; 标题
	res@gsnCenterString                 ="20180627 EastAsian pre"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF                    =15
    res@mpMaxLatF                    =45
	res@mpMinLonF                    =100
    res@mpMaxLonF                    =140
	
	; 高分辨率地图
	res@mpDataBaseVersion                  ="HighRes"
	
	; 色板
	res@cnFillPalette                       ="BlueDarkRed18"
	
	; 等值线间距
	; res@cnLevelSelectionMode             ="AutomaticLevels"
    res@cnLevelSelectionMode               ="ManualLevels"
	res@cnMaxLevelValF                     =70
    res@cnMinLevelValF                     =0
	res@cnLevelSpacingF                    =10
	
    res@mpGeophysicalLineColor             ="gray"

    ; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	res@tmXBLabelFontHeightF               =0.02 
    ; 主标题: tiMainFontHeightF
    ; res@tiMainFontHeightF                  =0.02
    ; res@tiMainFontThicknessF               =0.02
    res@tiXAxisFontHeightF                 =0.02
    res@tiYAxisFontHeightF                 =0.02
	; 副标题: gsnCenterStringFontHeightF
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode             = "ExplicitLevels"
    rescn@cnLevels                         = 30

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor                      ="black"
	rescn@cnLineThicknessF                 =1.8
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn                   =False
	rescn@cnInfoLabelOn                    =False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	resvc@vcGlyphStyle                     ="LineArrow"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF            =0.8
	resvc@vcMinDistanceF                   =0.025

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF                     =0.1
	resvc@vcRefMagnitudeF                  =20
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoOrthogonalPosF          =-0.05
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	
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
	resxy@tiMainString                  ="20180615-20180715 35N127.5E pre"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString                 =""
    resxy@tiYAxisString                 =""
	
	; gsnXYBarChart
	resxy@gsnXYBarChart                 =True
	
	; X 轴
	; tmXBMode 及相应属性
    timexy =ispan(1, 31, 1)

    resxy@tmXBMode                       ="Explicit"
    resxy@tmXBValues                     =ispan(2, 28, 5)
	resxy@tmXBLabels                     =(/"0616","0621","0626","0701","0706","0711"/)
    resxy@tmXBMinorOn                    =True
    resxy@tmXBMinorValues                =timexy

    ; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors            ="SlateBlue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF              =0.02
    resxy@tiMainFontThicknessF           =0.02
	resxy@tiXAxisFontHeightF             =0.02
	resxy@tiYAxisFontHeightF             =0.02

	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    resmarker                          =True
	resmarker@gsMarkerIndex            =12
    resmarker@gsMarkerColor            ="white"
	resmarker@gsMarkerSizeF            =0.02
	resmarker@gsMarkerThicknessF       =4
	; =========================================================================
	; 画图
	
	; 图 (a)
	; 降水地图
	plota =gsn_csm_contour_map(wks, tp1_1, res)
	
	; 降水=30 mm
	plotcn =gsn_csm_contour(wks, tp1_1, rescn)
	
	; 700百帕风场
	plotvc =gsn_csm_vector(wks, uwnd1, vwnd1, resvc)
	
	; 添加标记
	plotmarker =gsn_add_polymarker(wks, plota, 127.5, 35, resmarker)
	
    overlay(plota,plotcn)   
    overlay(plota,plotvc)	
	; --------------------------------------------------------------------------
	; 图 (b)
    plotb =gsn_csm_xy(wks, timexy, tp2_1, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)

    resPanel                            = True 
    resPanel@gsnPanelMainString         = "18304071 shaoyuqian"
    resPanel@gsnPanelMainFontHeightF    = 0.02
	
    plot =new(2, graphic)
	plot(0) =plota
    plot(1) =plotb
	gsn_panel(wks, plot, (/2,1/), resPanel)


	; =========================================================================
	; 去掉空白空间

    delete (wks)          
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    ; system("rm " + pltName + ".eps -f") 


end