begin

	; 读取 .nc 文件
	f=addfile("tp.assignment3.nc", "r")
	f1=addfile("uwnd.assignment3.nc","r")
	f2=addfile("vwnd.assignment3.nc","r")
	tp 			=f ->tp
	uwnd		=f1 ->uwnd
	vwnd 		=f2 ->vwnd
	time 		=f ->time
	time1		=f1 ->time
	time2 		=f2 ->time		
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	tt=cd_calendar(time, 2)
	tt1=cd_calendar(time1, 2)
	tt2=cd_calendar(time2, 2)
	; 2018年6月27日
	id=ind(tt .eq. 20180627)
	id1=ind(tt1 .eq. 20180627)
	id2=ind(tt2 .eq. 20180627)
	;print(id1)
	;print(id2)
	; 2018年6月15日至7月15日
	index=ind(tt .ge. 20180615 .and. tt .le. 20180715)
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	water1=tp(id,:,:  )
	water1=water1*1000
	;;printVarSummary(water1)
	;printMinMax(water1, -1)
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	water2=tp(time|index,{latitude|35}, {longitude|127.5} )
	water2=water2*1000
    ;print(water2)
	;printMinMax(water2, -1)
	; 提取 2018年6月27日 700百帕的风场
    u=uwnd(time|1734,{level|700},lat|:,lon|:)
    v=vwnd(time|1734,{level|700},lat|:,lon|:)
    ;; =========================================================================
    pltName="18304102.assignment3.q1"
    wks=gsn_open_wks("eps", "18304102.assignment3.q1")
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
	res@gsnCenterString="Rainfall(20180627)-wind(700hpa)"
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF=15
	res@mpMaxLatF=45
	res@mpMinLonF=100
	res@mpMaxLonF=140
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"
	; 色板
	; cnFillPalette
	res@cnFillPalette="cmocean_deep"
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode="ManualLevels"
	res@cnMinLevelValF=0
	res@cnMaxLevelValF=160
	res@cnLevelSpacingF=8
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF=0.025
	res@tmYLLabelFontHeightF=0.025
	res@gsnCenterStringFontHeightF=0.03
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	;调色板位置及格式
    res@lbBoxEndCapStyle="TriangleBothEnds"
	;res@pmLabelBarOrthogonalPosF=
    res@pmLabelBarParallelPosF=0.38
	res@pmLabelBarOrthogonalPosF=0.1
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode="ExplicitLevels"
	rescn@cnLevels=30.
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor="black"
	rescn@cnLineThicknessF=3.
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn=False
	rescn@cnInfoLabelOn=False
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle="CurlyVector"
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF=1.5
	resvc@vcMinDistanceF=0.025

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF=0.1
	resvc@vcMaxMagnitudeF=20
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF=1.05
	resvc@vcRefAnnoOrthogonalPosF=-0.23
	
	
	; (可选) 参考矢量信息的字体大小
	; res@vcRefAnnoFontHeightF
	; res@vcRefAnnoFontHeightF

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
	resxy@tiMainString="20180615-20180715_rainfall"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="Time(day)"
	resxy@tiYAxisString="Rainfall(mm)"
	; gsnXYBarChart
	resxy@gsnXYBarChart=True
	; X 轴
	; tmXBMode 及相应属性/)
    x=ispan(1,31,1)
	resxy@tmXBMode="Explicit"
	resxy@tmXBValues=(/2,7,13,18,23,28/)
	resxy@tmXBLabels=(/"6.16","6.21","6.26","7.1","7.6","7.11"/)
    resxy@tmXBMinorValues=x
    ;resxy@tmXBMinorValues=x
	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors="Red"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF=0.03
	resxy@tiXAxisFontHeightF=0.03
	resxy@tiYAxisFontHeightF=0.03
    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	x1=127.5
	y1=35
	resmark=True
	resmark@gsMarkerIndex=13
	resmark@gsMarkerColor="Red"
	resmark@gsMarkerSizeF=0.02
	; =========================================================================
	; 画图
	
	; 图 (a)
	plot=new(2, graphic)
	; 降水地图
	plot(0)=gsn_csm_contour_map(wks, water1, res)
	
	; 降水=30 mm
	plot1=gsn_csm_contour(wks, water1, rescn)
	overlay(plot(0), plot1)
	; 700百帕风场
	plot2=gsn_csm_vector(wks, u, v, resvc)
	overlay(plot(0),plot2)
	; 添加标记
	plotm=gsn_add_polymarker(wks, plot(0), x1, y1, resmark)
	; --------------------------------------------------------------------------
	; 图 (b)
	plot(1)=gsn_csm_xy(wks, x,water2, resxy)
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	respanel=True
    respanel@gsnPanelMainString="ZhouXuan-18304102"
    respanel@gsnPanelMainFontHeightF=0.025
	; respanel@gsnPanelXF=
	; respanel@gsnPanelYF=
	gsn_panel(wks, plot, (/2,1/), respanel)
	; =========================================================================
	; 去掉空白空间
	delete (wks) 			
	system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
	system("rm " + pltName + ".eps -f") 
end