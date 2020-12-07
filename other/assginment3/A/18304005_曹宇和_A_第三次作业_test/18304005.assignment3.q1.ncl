begin

	; 读取 .nc 文件
	path 	="/home/ys17-19/lsh/Project/other/assginment3/"
	f1 		=addfile(path+"tp.assignment3.nc", "r")
	f2 		=addfile(path+"uwnd.assignment3.nc", "r")
	f3 		=addfile(path+"vwnd.assignment3.nc", "r")
	
	; print(f1)
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time1 	=f1->time
	time2 	=f2->time
	time3 	=f3->time

	; print(time1)
	
	time_tp 	= cd_calendar(time1, 2)
	time_uwnd	= cd_calendar(time2, 2)
	time_vwnd	= cd_calendar(time3, 2)
	
	time0627	= ind(time_tp .eq. 20180627)
	; print(time0627)
	
	; 2018年6月15日至7月15日
	
	time06_0715 	=ind(time_tp .ge. 20180615 .and. time_tp .le. 20180715)
	; print(time06_0715)
	
	time_tp1 		= cd_calendar(time1, 0)
	time_tp_m 		=tointeger(time_tp1(502:532,1))
	time_tp_d		=tointeger(time_tp1(502:532,2))
	time06_0715a 	=sprinti("%0.2i", time_tp_m)+sprinti("%0.2i", time_tp_d)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	tp0627 	= f1->tp(time0627,:,:)
	tp0627 	=tp0627*1000
	; printMinMax(tp0627, True)

	; m1=f1->latitude
	; m2=f1->longitude
	; m11=ind(m1 .eq. 35)
	; m22=ind(m2 .eq. 127.5)
	; print(m11)
	; print(m22)



	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm

	tp06_0715 		= f1->tp(502:532,100,150)
	tp06_0715		=1000*tp06_0715
	
	
	
	; 提取 2018年6月27日 700百帕的风场

	time0627_u	= ind(time_uwnd .eq. 20180627)
	time0627_v	= ind(time_vwnd .eq. 20180627)
	; print(time0627_u)
	uwnd 	=f2->uwnd(time0627_u,1,:,:)
	vwnd 	=f3->vwnd(time0627_v,1,:,:)

	
	
	; ========================================================================
	; 打开工作空间
	
	wks 	=gsn_open_wks("png", "18304005.assignment3.q1")
	
	
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
	res@gsnCenterString            = "East Asia tp(mm) and 700hPa wind(ms~S~-1~N~)"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF 	=15.
	res@mpMaxLatF 	=45.
	res@mpMinLonF 	=100.
	res@mpMaxLonF 	=140.
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	res@mpDataBaseVersion  	="HighRes"
	
	; 色板
	; cnFillPalette
	
	res@cnFillPalette 		="wind_17lev"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode 	="ManualLevels"
	res@cnMaxLevelValF 			=200.
	res@cnMinLevelValF 			=0.
	res@cnLevelSpacingF 		=20.
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF 		=0.015
	res@tmYLLabelFontHeightF 		=0.015
	res@gsnCenterStringFontHeightF 	=0.02
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
	rescn@cnLevelSelectionMode 	="ExplicitLevels"
	rescn@cnLevels 				=30


	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	
	rescn@cnLineThicknessF 		=4
	rescn@cnLineColor  			="black"
	

	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnLineLabelsOn 	=False
	rescn@cnInfoLabelOn 	=False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle 		="LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	resvc@vcLineArrowThicknessF = 1.5
	resvc@vcMinDistanceF 		= 0.025
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefLengthF 		=0.1
	resvc@vcRefMagnitudeF 	=18.

	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	
	resvc@vcRefAnnoParallelPosF 	=1.14
	resvc@vcRefAnnoOrthogonalPosF 	=-0.24
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
	
	resxy@tiMainString 	="35~S~o~N~N,~N~127.5~S~o~N~E 0615-0715 precipitation(mm)"
	resxy@tiMainFontHeightF 	=0.02


	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString				= "time(day)"
	resxy@tiYAxisString 			= "precipitation(mm)"

	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 			= True
	
	
	; X 轴
	; tmXBMode 及相应属性

	t=ispan(1,31,1)
	resxy@tmXBMode 		="Explicit"
	resxy@tmXBValues 	=t(1::5)
	resxy@tmXBLabels 	=time06_0715a(1::5)
	resxy@tmXBMinorOn 	=True
	resxy@tmXBMinorValues = t 
	

	



	; (可选) XY线及柱状的颜色

	resxy@gsnAboveYRefLineColor 	= "Blue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiXAxisFontHeightF 		= 0.02
	resxy@tiYAxisFontHeightF 	 	= 0.02
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
	resMark 					=True
	resMark@gsMarkerIndex 		=8
	resMark@gsMarkerSizeF 		=0.025
	resMark@gsMarkerColor 		="Red"
	resMark@gsMarkerThicknessF 	=4.

	
	; =========================================================================
	; 画图

	resPanel 							= True
	resPanel@gsnPanelMainString 		= "Cao Yuhe(18304005)"
	resPanel@gsnPanelMainFontHeightF 	= 0.025
	resPanel@gsnPanelLabelBar 			= False
	resPanel@lbLabelFontHeightF 		= 0.025
	resPanel@gsnPanelYWhiteSpacePercent = 5.
	
	; 图 (a)
	
	; 降水地图
	
	
	; 降水=30 mm
	
	
	; 700百帕风场
	
	
	; 添加标记
	
	plot 	=new(2,graphic)

	plot1 	=gsn_csm_contour_map_ce(wks, tp0627, res)
	plotcn 	=gsn_csm_contour(wks, tp0627, rescn)
	plotvc	=gsn_csm_vector(wks, uwnd, vwnd, resvc)
	plot2 	=gsn_add_polymarker(wks, plot1, 127.5, 35, resMark)

	
	overlay(plot1,plotcn)
	overlay(plot1,plotvc)
	
	plot(0) 	=plot1
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	
	plot_z 	=gsn_csm_xy(wks, ispan(1,31,1), tp06_0715, resxy)
	
	plot(1) 	=plot_z
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
	gsn_panel(wks,plot,(/2,1/),resPanel)
  	 
	
	
	; =========================================================================
	; 去掉空白空间
	; delete (wks)
	; system("convert -geometry 3000x3000 -density 600x600 -trim " + "18304005.assignment3.q1" + ".eps " + "18304005.assignment3.q1" + ".png")   
	; system("rm " + "18304005.assignment3.q1" + ".eps -f") 



end