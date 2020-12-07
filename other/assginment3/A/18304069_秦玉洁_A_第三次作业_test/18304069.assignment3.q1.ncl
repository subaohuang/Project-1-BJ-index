begin

	; 读取 .nc 文件
     f1=addfile("D:\NCLexercise\tp.assignment3.nc","r")
     f2=addfile("D:\NCLexercise\uwnd.assignment3.nc","r")
     f3=addfile("D:\NCLexercise\vwnd.assignment3.nc","r")
     time=f1->time
     time1=f2->time
     precipitation=f1->tp
     uwnd=f2->uwnd
     vwnd=f3->vwnd
     ;printVarSummary(precipitation)
     ;print(uwnd)
     ;print(vwnd)

   
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	timeConvert=cd_calendar(time, -2)
	time1_ind=ind(timeConvert.eq.20180627);降雨时间
	;print(time1_ind)
	
	timeConvert1=cd_calendar(time1,-2)
	time3_ind=ind(timeConvert1.eq.20180627);风场时间
	;print(time3_ind)
	; 2018年6月15日至7月15日
	time2_ind=ind(timeConvert.ge.20180615.and.timeConvert.le.20180715);降雨时间
	;print(time2_ind)
  
     
   
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm

	tp1=precipitation(time1_ind,:,:)
	tp1=(/tp1*1000./)
	;printMinMax(tp1, 1)
	;提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2=precipitation(time|time2_ind,{latitude|35},{longitude|127.5})
	tp2=(/tp2*1000/)
	;print(tp2)
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	uwnd1=uwnd(time|time3_ind,{level|700},lat|:,lon|:)
	vwnd1=vwnd(time|time3_ind,{level|700},lat|:,lon|:)
	; printVarSummary(uwnd1)
	; printVarSummary(vwnd1)
	; printVarSummary(uwnd1)
	
	
	; ========================================================================
	; 打开工作空间
	wksName="18304069.assignment3.q1"
	wks=gsn_open_wks("eps", wksName)
	
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
	res@cnFillOn=True
	; 标题
	res@tiMainString="precipitation of east Asia at 2018-06-27 and wind at 700hPa"
	
	
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
	
	res@cnFillPalette="precip2_17lev"

	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode="ManualLevels"
	res@cnMinLevelValF=0
	res@cnMaxLevelValF=100
	res@cnLevelSpacingF=10
	
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF=0.02
	res@tmYLLabelFontHeightF=0.02
    res@tiMainFontHeightF=0.03
  
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF=0.02
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	rescn=True

	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode="ExplicitLevels"
	rescn@cnLevels=30
	
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor="Black"
	rescn@cnLineThicknessF=2.5
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn=False
	rescn@cnInfoLabelOn=False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle="LineArrow"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	resvc@vcLineArrowThicknessF=0.02
	resvc@vcMinDistanceF       =0.02
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF=0.05
	resvc@vcRefMagnitudeF=20
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF=1.0
	resvc@vcRefAnnoOrthogonalPosF=0.08
	
	
	; (可选) 参考矢量信息的字体大小
	resvc@vcRefAnnoFontHeightF=0.02

	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.55
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	
	resxy@tiMainString="precipitation at 2018-06-15 (35~S~o~N~N 127.5~S~o~N~E)"
	resxy@tiMainFont="helvetica-bold"

	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="time"
	resxy@tiYAxisString="precipitation(mm)"
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart=True
	
	
	; X 轴
	; tmXBMode 及相应属性
	x=time(time2_ind)
    resxy@tmXBMode="Explicit"
	resxy@tmXBValues=x(::5)
	resxy@tmXBLabels=(/"0615","0620","0625","0630","0705","0710","0715"/)
	resxy@tmXBMinorOn=True
	resxy@tmXBMinorValues=x
	

    ; (可选) XY线及柱状的颜色
	resxy@xyLineColor="Blue"
	colors=read_colormap_file("cmocean_ice")
	resxy@gsnXYBarChartColors=colors(200,:)
	; (可选) 字体大小s
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF=0.03
	resxy@tiXAxisFontHeightF=0.02
	resxy@tiYAxisFontHeightF=0.02
	resxy@tiMainFontThicknessF=5
	resxy@tiMainOffsetYF 	= 0.02
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark                   =True
	resMark@gsMarkerColor     ="Red"
	resMark@gsMarkerIndex     =15
	resMark@gsMarkerSizeF     =0.02
	resMark@gsEdgeThicknessF  =2
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plotBase 				= gsn_csm_contour_map_ce(wks,tp1,res)
  
	
	; 降水=30 mm
	plotOverlay             = gsn_csm_contour(wks, tp1, rescn)
	overlay(plotBase, plotOverlay)
	; 700百帕风场
	plotOverlay1            = gsn_csm_vector(wks, uwnd1, vwnd1, resvc)
	overlay(plotBase, plotOverlay1)
	; 添加标记
	
	plotMarker              =gsn_add_polymarker(wks, plotBase, 127.5, 35, resMark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plot 					= gsn_csm_xy(wks,x,tp2,resxy)

	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel=True
	resPanel@gsnPanelXF=(/0.10,0.60/)
	resPanel@gsnPanelYF=(/0.66,0.66/)
	;resPanel@gsnPanelXF=-0.2
	;resPanel@gsnPanelDebug=True
	resPanel@tiMainString="18304069qinyujie"
	gsn_panel(wks,(/plotBase,plot/),(/1,2/),resPanel)
	; =========================================================================
	; 去掉空白空间
	delete(wks)
    system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
	system("/bin/rm -f " + wksName + ".eps")
end