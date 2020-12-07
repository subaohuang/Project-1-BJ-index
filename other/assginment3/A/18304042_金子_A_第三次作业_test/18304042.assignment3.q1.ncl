begin

	; 读取 .nc 文件

	ftp 			= addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc","r")
	fu 				= addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc","r")
	fv 				= addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc","r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日

	time_tp 		= cd_calendar(ftp->time, -2)
	time_u 			= cd_calendar(fu->time, -2)
	time_v 			= cd_calendar(fv->time, -2)

	ind_tp1 		= ind(time_tp.eq.20180627)
	ind_u   	    = ind(time_u.eq.20180627)
	ind_v   	    = ind(time_v.eq.20180627)
	
	; 2018年6月15日至7月15日
	
	ind_tp2 		= ind(time_tp.ge.20180615.and.time_tp.le.20180715)
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	tp  			= ftp->tp(:,:,:)
	tp1 			= ftp->tp(ind_tp1,:,:)*1000
	copy_VarCoords(tp(ind_tp1,:,:), tp1)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	tp2  			= ftp->tp(ind_tp2,{35},{127.5})*1000
	copy_VarCoords(tp(ind_tp2,0,0), tp2)
	
	; 提取 2018年6月27日 700百帕的风场
	
	u 				= fu->uwnd(ind_tp1,{700},:,:)
	v 				= fv->vwnd(ind_tp1,{700},:,:)
	
	; ========================================================================
	; 打开工作空间
	
	wksName 		= "18304042.assignment3.q1"
	wks 			= gsn_open_wks("eps",wksName )
	
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

	res@gsnLeftString		= "mm"
	res@gsnRightString		= "m/s"
	
	; 标题
	
	res@tiMainString 		= " (a) AAR & wind 700hPa (EA 20180607) "
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
	res@mpMinLatF   		= 15
	res@mpMaxLatF 			= 45
	res@mpMinLonF 			= 100
	res@mpMaxLonF 			= 140
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	; res@mpDataBaseVersion 	= "HighRes"
	res@mpDataBaseVersion  =  "HighRes"
	res@mpGeophysicalLineColor  =  "gray60"
	res@mpGeophysicalLineThicknessF  =  1.
	
	; 色板
	; cnFillPalette
	
	res@cnFillPalette 		= "gui_default"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
	res@cnLevelSelectionMode 			= "ManualLevels"
	res@cnMaxLevelValF      			= 40
	res@cnMinLevelValF      			= 0
	res@cnLevelSpacingF     			= 5
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	
	res@tmXBLabelFontHeightF 			= 0.02
	res@tmYLLabelFontHeightF 			= 0.02
	res@gsnLeftStringFontHeightF 		= 0.02
	res@gsnRightStringFontHeightF 		= 0.02
	res@tiMainFontHeightF 				= 0.025
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	res@tiMainOffsetYF 					= 0.01
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
	rescn@cnLevelSelectionMode 			= "ExplicitLevels"
	rescn@cnLevels 						= (/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	
	rescn@cnLineColor 					= "black"
	rescn@cnLineThicknessF 				= 2.5
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnInfoLabelOn         = False
    rescn@cnLineLabelsOn        = False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	
	resvc@vcGlyphStyle 			= "LineArrow"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	resvc@vcMinDistanceF        = 0.02
	resvc@vcLineArrowThicknessF = 1.5

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefLengthF           = 0.1
	resvc@vcRefMagnitudeF        = 25
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	
	resvc@vcRefAnnoParallelPosF   = 1.0
	resvc@vcRefAnnoOrthogonalPosF = 0.01
	resvc@vcRefAnnoString1        = "$VMG$ m/s"
	resvc@vcRefAnnoString2On      = False
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	resvc@vcRefAnnoFontHeightF    = 0.015
	
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
	
	resxy@tiMainString = " (b) AAR at 35~S~o~N~N 127.5~S~o~N~E (20180615-20180715)"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString 				= "time"
	resxy@tiYAxisString 				= "AAR"
	resxy@gsnLeftString 				= "mm"
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart 				= True
	
	; X 轴
	; tmXBMode 及相应属性

	resxy@tmXBMode 				= "Explicit"
	resxy@tmXBValues 			= (/0,5,10,15,20,25,30/)
	resxy@tmXBLabels 			= (/"6.14","6.19","6.24","6.29","7.04","7.09","7.14"/)
	time_xy 					= ispan(1,31,1)
	resxy@tmXBMinorValues 		= time_xy
	
	; (可选) XY线及柱状的颜色
	
	resxy@gsnXYBarChart 		= True
	resxy@gsnXYBarChartColors 	= "Blue"
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	resxy@tiMainFontHeightF 			= 0.025
	resxy@tiXAxisFontHeightF 			= 0.025
	resxy@tiYAxisFontHeightF 			= 0.025
	resxy@tmXBLabelFontHeightF 			= 0.02
	resxy@tmYLLabelFontHeightF 			= 0.02
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
	resMark 							= True
	resMark@gsMarkerIndex 				= 14
	resMark@gsMarkerThicknessF 			= 3.5
	resMark@gsMarkerColor 				= "Yellow"
	resMark@gsMarkerSizeF 				= 0.015
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图

	plot   = gsn_csm_contour_map_ce(wks, tp1, res)
	
	; 降水=30 mm

	plotcn = gsn_csm_contour(wks, tp1, rescn)
	
	; 700百帕风场

	plotvc = gsn_csm_vector(wks, u, v, resvc)
	
	; 添加标记
	
	plotM  = gsn_add_polymarker(wks, plot, 127.5, 35, resMark)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	
	plotxy = gsn_csm_xy(wks, time_xy, tp2, resxy)
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
	resPanel                            = True
    resPanel@gsnPanelYWhiteSpacePercent = 1.8
    resPanel@gsnPanelXF                 = 2.5
    resPanel@gsnPanelYF                 = 0.6

    resPanel@gsnPanelMainString         = "183048042 JINZI"
    resPanel@gsnPanelMainFontHeightF    = 0.02
    
    overlay(plot, plotvc)
	overlay(plot, plotcn)
	gsn_panel(wks,(/plot,plotxy/),(/2,1/),resPanel)
	
	; =========================================================================
	; 去掉空白空间
	
	delete(wks)
	system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
	system("rm " + wksName + ".eps -f") 

end