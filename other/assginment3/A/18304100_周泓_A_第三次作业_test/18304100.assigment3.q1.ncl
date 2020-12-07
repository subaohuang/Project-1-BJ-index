begin

	; 读取 .nc 文件
    path                                ="F:\NCLlesson\ncl\"
    file1                               ="tp.assignment3.nc"
    file2                               ="uwnd.assignment3.nc"
    file3                               ="vwnd.assignment3.nc"
	f1                                  =addfile(path+file1, "r")
    f2                                  =addfile(path+file2, "r")
    f3                                  =addfile(path+file3, "r")

	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	
	time1                               =f1->time
    time1A                              =cd_calendar(time1, -2)
    time1ind1                           =ind(time1A .eq. 20180627)

    time2                               =f2->time
    time2A                              =cd_calendar(time2, -2)
    time2ind1                           =ind(time2A .eq. 20180627)

	time3                               =f3->time
    time3A                              =cd_calendar(time3, -2)
    time3ind1                           =ind(time3A .eq. 20180627)

	; 2018年6月15日至7月15日
	
	time1ind2                           =ind(time1A .ge. 20180615 .and. time1A .le. 20180715)
    ;time1B1                             =f1->time(time1ind2)
    ;time1B2                             =cd_calendar(time1B1, -2)-20180000
    ;print(time1B2)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm

	tp1                                 =flt2dble(f1->tp(time1ind1,:,:))
    ; printVarSummary(tp1)
	tp11                                =tp1*1000.
	;printMinMax(tp11, 0)
    copy_VarMeta(tp1, tp11)
    tp11@units                          ="mm"

	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	lat1                                =f1->latitude
    lon1                                =f1->longitude
    lat1ind                             =ind(lat1 .eq. 35)
    ;print(lat1ind)
    lon1ind                             =ind(lon1 .eq. 127.5)
    ;print(lon1ind)	
    tp2                                 =f1->tp(time1ind2,lat1ind,lon1ind)
    ;printMinMax(tp2, 0)
    tp22                                =1000*tp2
    copy_VarMeta(tp2, tp22)
    tp22@units                          ="mm"
    ;printVarSummary(tp22)
	
	; 提取 2018年6月27日 700百帕的风场
	
	level1                              =f2->level
    level1ind1                          =ind(level1 .eq. 700)
    ;print(level1ind1)
    u1                                  =f2->uwnd(time2ind1,level1ind1,:,:)
    ;print(u1)

    level2                              =f3->level                          
	level2ind1                          =ind(level2 .eq. 700)
    ;print(level2ind1)
    v1                                  =f3->vwnd(time3ind1,level2ind1,:,:)
    ;print(v1)
	
	; ========================================================================
	; 打开工作空间
	
	wks                                 =gsn_open_wks("eps", "18304100.assignment3.q1")
	
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
	
	res@tiMainString                    ="Tp and Wind of 2018.6.27"
	res@tiMainFontHeightF               =0.015

	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
	res@mpMinLatF                       =15
    res@mpMaxLatF                       =45
    res@mpMinLonF                       =100
    res@mpMaxLonF                       =140
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	res@mpDataBaseVersion               ="HighRes"
	
	; 色板
	; cnFillPalette
	
	res@cnFillPalette                   ="precip4_11lev"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
    res@cnLevelSelectionMode            ="ManualLevels"
    res@cnMaxLevelValF                  =20
    res@cnMinLevelValF                  =0
    res@cnLevelSpacingF                 =4
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
    rescn@cnLevelSelectionMode          ="ManualLevels"
    rescn@cnMaxLevelValF                =30
    rescn@cnLevelSpacingF               =6
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	
	rescn@cnLineColor                   ="Pink"
    ;rescn@cnLineThicknessesF            =0.01
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnLineLabelsOn                =False
    rescn@cnInfoLabelOn                 =False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	
	resvc@vcGlyphStyle                  ="CurlyVector"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	resvc@vcLineArrowThicknessF         =1.5
    resvc@vcMinDistanceF                =0.01	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefLengthF                  =0.1
    resvc@vcRefMagnitudeF               =30
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	
	resvc@vcRefAnnoParallelPosF         =0.85
    resvc@vcRefAnnoOrthogonalPosF       =-0.015
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	resvc@vcRefAnnoFontHeightF          =0.015
	
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
	
	resxy@tiMainString                  ="TP from 20180615 to 20180715 on 35~S~o~N~N 127.5~S~o~N~E"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiXAxisString                 ="Date"
    resxy@tiYAxisString                 ="Total Precipation (mm)"
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart                 =True
	
	; X 轴
	; tmXBMode 及相应属性

    date2                               =ispan(1,31,1)
    resxy@tmXBMode                      ="Explicit"
    resxy@tmXBMinorOn                   =True
    resxy@tmXBMinorValues               =date2    
    resxy@tmXBValues                    =ispan(2,27,5)
    resxy@tmXBLabelsOn                  =True
    resxy@tmXBLabels                    =(/"06.16","06.21","06.26","07.01","07.06","07.11"/) 

	; (可选) XY线及柱状的颜色
	
	resxy@gsnXYBarChartColors2           ="BlueDarkRed18"
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	resxy@tiMainFontHeightF             =0.015
    resxy@tiXAxisFontHeightF            =0.015
    resxy@tiYAxisFontHeightF            =0.015
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF

	resmarker                           =True
	resmarker@gsMarkerIndex             =16
    resmarker@gsMarkerColor             ="Red"
    resmarker@gsMarkerSizeF             =0.01
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	
    plota1                             =gsn_csm_contour_map_ce(wks, tp11, res)
	
	; 降水=30 mm
	
    plottp                             =gsn_csm_contour(wks, tp11, rescn)
	overlay(plota1, plottp)
	; 700百帕风场
	
	plota2                             =gsn_csm_vector(wks, u1, v1, resvc)
    overlay(plota1, plota2)

	; 添加标记
	
	addplot                            =gsn_add_polymarker(wks, plota1, 127.5, 35, resmarker)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	
	plotb                              =gsn_csm_xy(wks, date2, tp22, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
    resPanel                            =True
    resPanel@gsnPanelMainString         ="ZhouHong (18304100)"
    resPanel@gsnPanelMainFontHeightF    =0.015
    resPanel@gsnPanelLabelBar           =False
    resPanel@gsnPanelYWhiteSpacePercent =1
    resPanel@gsnPanelYF                 =(/-0.1,0.65/)
	gsn_panel(wks, (/plota1,plotb/), (/1,2/), resPanel)
	
	; =========================================================================
	; 去掉空白空间
	
    delete(wks)

    system("convert -geometry 3000x3000 -density 600x600 -trim " + "18304100.assignment3.q1" + ".eps " + "18304100.assignment3.q1" + ".png")   
    system("rm " + "18304100.assignment3.q1" + ".eps -f") 

end