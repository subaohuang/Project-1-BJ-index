begin

	; 读取 .nc 文件
	fp  =  addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc", "r")
	fu  =  addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc","r")
	fv  =  addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc","r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time  =  fp->time
	timea  =  cd_calendar(time, -2)
	timeind1  =  ind(timea .eq. 20180627)
	
	
	; 2018年6月15日至7月15日
	timeind2  =  ind(timea .ge. 20180615 .and. timea .le. 20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1  =  fp->tp(timeind1,:,:)
	tpmm  =  tp1*1000.
	copy_VarCoords(tp1, tpmm)
	

	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2  =  fp->tp(time|timeind2,{latitude|35.},{longitude|127.5})
	tpmm2  =  tp2*1000.

	
	; 提取 2018年6月27日 700百帕的风场
	uw  =  fu->uwnd(timeind1,{700},:,:)
	vw  =  fv->vwnd(timeind1,{700},:,:)
	
	
	; ========================================================================
	; 打开工作空间
	pltname  =  "18304050.assignment3"
	wks  =  gsn_open_wks("eps", pltname)
	
	
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
	res@tiMainString  =  "(a)"
	res@gsnCenterString  =  "Precipitation in Jun 27,2018 (wnd700hPa)"
	res@gsnCenterStringFontHeightF  =  0.02
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF  =  15
	res@mpMaxLatF  =  45
	res@mpMinLonF  =  100
	res@mpMaxLonF  =  140
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion  =  "HighRes"
	res@mpGeophysicalLineColor  =  "gray60"
	res@mpGeophysicalLineThicknessF  =  1.

	; 色板
	; cnFillPalette
	res@cnFillPalette  =  "precip3_16lev"
	;res@gsnSpreadColorStart   =  2
	;res@gsnSpreadColorEnd     =  10
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode  =  "ManualLevels"
	res@cnMaxLevelValF  =  100.
	res@cnMinLevelValF  =  0.
	res@cnLevelSpacingF =  10.
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF  =  0.025
	res@tmYLLabelFontHeightF  =  0.025
	res@tiMainFontHeightF  =  0.03
	res@gsnCenterStringFontHeightF  =  0.028
	

	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF  =  0.
	res@lbBoxEndCapStyle  =  "TriangleBothEnds"
	res@lbOrientation  =  "Vertical"
	res@vpWidthF  =  0.7

	;res@pmLabelBarOrthogonalPosF  = 0.05
	;res@pmLabelBarParallelPosF  =  0.6
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode  =  "ExplicitLevels"
	rescn@cnLevels  =  30.

	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor  =  "black"
	rescn@cnLineThicknessF  =  2.
	rescn@cnLineDrawOrder  =  "PostDraw"
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn  =  False
	rescn@cnInfoLabelOn  =  False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle  =  "CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF  =  1.2
	resvc@vcMinDistanceF  =  0.025  ; ??
	resvc@vcLineArrowColor  =  "green4"

	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF  =  0.06
	resvc@vcRefMagnitudeF  = 15
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF  =  0.1
	resvc@vcRefAnnoOrthogonalPosF  =  -1.
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	resvc@vcRefAnnoFontHeightF  =  0.012
	resvc@vcRefAnnoString2On  =  False
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 			   			= 0.7
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString  =  "(b)"
	resxy@gsnCenterString  =  "Precipitation in Jun 15-Jul 15, 2018, 35~S~o~N~N,127.5~S~o~N~E"
	resxy@gsnCenterStringFontHeightF  =  0.028
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString  =  "time"
	resxy@tiYAxisString  =  "total_precipitation (mm)"
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart  =  True
	
	
	; X 轴
	; tmXBMode 及相应属性
	resxy@tmXBMode  =  "Explicit"
	resxy@tmXBValues  =  (/1,6,11,16,21,26/)
	resxy@tmXBLabels  =  (/"Jun16","Jun21","Jun26","Jul01","Jul06","Jul11"/)
	resxy@tmXBMinorValues  =  fspan(0,30,31)

	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors  =  "royalblue"
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF  =  0.03
	resxy@tiXAxisFontHeightF  =  0.025
	resxy@tiYAxisFontHeightF  =  0.025
	resxy@tmXBLabelFontHeightF  =  0.02
	resxy@tmYLLabelFontHeightF  =  0.02
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resmark  =  True
	resmark@gsMarkerIndex  =  8
	resmark@gsMarkerColor  =  "indianred1"
	resmark@gsMarkerSizeF  =  0.015
	resmark@gsMarkerThicknessF  =  5.
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot  =  gsn_csm_contour_map_ce(wks, tpmm,res)
	
	; 降水=30 mm
	plotcn  =  gsn_csm_contour(wks,tpmm,rescn)
	overlay(plot, plotcn)
	
	; 700百帕风场
	plotvc  =  gsn_csm_vector(wks, uw, vw, resvc)
	overlay(plot, plotvc)
	
	; 添加标记
	plotmarker  =  gsn_add_polymarker(wks, plot,127.5 , 35., resmark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)

	x1  =  fspan(0,30,31)
	plotxy  =  gsn_csm_xy(wks, x1, tpmm2, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	respanel  =  True
	respanel@gsnPanelXF  = 0.01
	respanel@gsnPanelYF  = 0.02
	respanel@gsnPanelMainString  =  "18304050liwenyi"

	plot2 = new(2,graphic)
	plot2(0)  =  plot
	plot2(1)  =  plotxy
	gsn_panel(wks,plot2,(/2,1/),respanel)
	
	
	; =========================================================================
	; 去掉空白空间
    delete (wks)
	system("convert -geometry 1000x1000 -density 600x600 -trim " + pltname + ".eps " + pltname + ".png")
	system("rm " + pltname + ".eps -f")
end
