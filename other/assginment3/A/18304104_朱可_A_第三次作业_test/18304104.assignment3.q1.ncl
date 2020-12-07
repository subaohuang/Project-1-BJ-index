begin

	; 读取 .nc 文件

    f1 = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc", "r")
    f2 = addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc", "r")
    f3 = addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc", "r")
	
    ;print(f1)	查看数据类型 （time,lat,lon）
    ;print(f2)  查看数据类型 （time,level,lat,lon）

    tp      = f1 ->tp
    uwnd    = f2 ->uwnd 
    vwnd    = f3 ->vwnd


    ;printVarSummary(tp) lat(60:0) lon(90:160) 纬度反向，需要转换
    ;printVarSummary(uwnd) lat(0:60) lon(90:150)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)


    time_tp = cd_calendar(f1->time, 2)
    time_u  = cd_calendar(f2->time, 2)
    time_v  = cd_calendar(f3->time, 2)

	; 2018年6月27日
	
	time_tp_0627      = ind(time_tp.eq.20180627)
	;print(time_tp_0627) 验证代码是否正确
    time_u_0627       = ind(time_u.eq.20180627)
    time_v_0627       = ind(time_v.eq.20180627)

	; 2018年6月15日至7月15日
	
	time_tp_67   = ind(time_tp.ge.20180615.and.time_tp.le.20180715)
    time_u_67    = ind(time_u.ge.20180615.and.time_u.le.20180715)   
    time_v_67    = ind(time_v.ge.20180615.and.time_v.le.20180715)

    time_uv_67    = time_u(time_u_67)
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	tp_0627 = tp(time|time_tp_0627,latitude|::-1,longitude|:)
    printVarSummary(tp_0627)
	tp_0627   = tp_0627*10^3
    printVarSummary(tp_0627)
    exit()
    printVarSummary(tp_0627)


	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	tp_67   = tp(time_tp_67,{35},{127.5})
    tp_67     = tp_67*10^3 
    print(tp_67)
	
	; 提取 2018年6月27日 700百帕的风场
	
	u_0627    = uwnd(time_u_0627,{700},:,:)
    v_0627    = vwnd(time_v_0627,{700},:,:) 
	
	; ========================================================================
	; 打开工作空间
	
    wksName   = "18304104.assignment3.q1"
	wks       = gsn_open_wks("eps", wksName)
	
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
	
    res@vpHeightF                       = 0.4 			
	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
    ;rescn@trGridType                      = "TriangularMesh"
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	
	res@gsnCenterString                 = "(a) 20180627 East Asian precipitation and 700hPa wind field"
    res@gsnLeftString                   = "precipitation/(mm)"
    res@gsnRightString                  = "wind/(ms~S~-1~N~)"
    res@gsnLeftStringFontColor          = "gray40"
    res@gsnRightStringFontColor         = "gray40"   


	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
	res@mpMinLatF                  = 15
    res@mpMaxLatF                  = 45
    res@mpMinLonF                  = 100
    res@mpMaxLonF                  = 140

	; 高分辨率地图
	; mpDataBaseVersion
	
	res@mpDataBaseVersion          = "HighRes"
    res@mpGeophysicalLineColor     = "gray40"
	
	; 色板
	; cnFillPalette
	
	res@cnFillPalette              = "precip_11lev" 
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
	printMinMax(tp_0627, 1)    ;明确6月27日降水范围
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMaxLevelValF       = 120
    res@cnMinLevelValF       = 0
    res@cnLevelSpacingF      = 10
	


	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	

    res@gsnCenterStringFontHeightF       = 0.024
	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	res@gsnCenterStringOrthogonalPosF   = 0.1 
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	
    rescn@cnLevelSelectionMode    = "ExplicitLevels"
	rescn@cnLevels                = 30

	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	
	rescn@cnLineColor             = "black"
    rescn@cnLineThicknessF        = 2.
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnLineLabelsOn         = False
    rescn@cnInfoLabelOn          = False
	

    ;rescn@trGridType                      = "TriangularMesh"
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	
	resvc@vcGlyphStyle           = "CurlyVector"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	resvc@vcLineArrowThicknessF  = 1
	resvc@vcMinDistanceF         = 0.02

    resvc@vcLineArrowColor       = "black"

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
    ; wind = wind_speed(u_0627, v_0627)
	; printMinMax(wind, 1)   明确6月27日风速的范围
    resvc@vcRefMagnitudeF    = 20
    resvc@vcRefLengthF       = 0.1
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	
	resvc@vcRefAnnoParallelPosF      = 1.1
    resvc@vcRefAnnoOrthogonalPosF    = -0.03
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.55
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题

    resxy@tiMainString                  = "(b) 20180615-20180715 precipitation at 35~S~o~N~N 127.5~S~o~N~E "
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	resxy@tiYAxisString                 = "precipitation / (mm)"
    resxy@tiXAxisString                 = "time / (day)" 
	
	; gsnXYBarChart
    
    resxy@gsnXYBarChart                 = True	
	

    x_time_67                          = ispan(-15, 15, 1)
    x_time_67_major                    = (/-14,-9,-4,1,6,11/)

	resxy@tmXBMode                     = "Explicit"
    resxy@tmXBValues                   = x_time_67_major
    resxy@tmXBLabels                   = (/"0616","0621","0626","0701","0706","0711"/)
	resxy@tmXBMinorOn                  = True
    resxy@tmXBMinorValues              = x_time_67


    ;X轴主刻度显示 6月16, 21, 26日及7月1, 6及11日，副刻度显示所有日子 

	; (可选) XY线及柱状的颜色
	
    resxy@xyLineColor                 = "black"
    resxy@gsnYRefLine                 = 0.0
    resxy@gsnAboveYRefLineColor       = "blue"
    resxy@gsnBelowYRefLineColor       = "red"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	res@tiMainFontHeightF             = 0.024
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
	resmark                     = True
    resmark@gsMarkerIndex       = 5
    resmark@gsMarkerColor       = "red"
    resmark@gsMarkerSizeF       = 0.005
    resmark@gsMarkerThicknessF  = 3.
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot_a  = gsn_csm_contour_map_ce(wks, tp_0627, res)

	; 降水=30 mm
	plot_a_cn = gsn_csm_contour(wks, tp_0627, rescn)
	
	; 700百帕风场
	
    plot_a_vc = gsn_csm_vector(wks, u_0627, v_0627, resvc)
	
	; 添加标记
	
	overlay(plot_a, plot_a_cn)
    overlay(plot_a, plot_a_vc)

    plot_a_marker = gsn_add_polymarker(wks, plot_a, 127.5, 35, resmark)



	; --------------------------------------------------------------------------
	; 图 (b)
	

	plot_b = gsn_csm_xy(wks, x_time_67, tp_67, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
    resPanel                            = True
    resPanel@gsnPanelBottom             = 0.15
    resPanel@gsnPanelYWhiteSpacePercent = 3.
    
    ; 标题
    resPanel@gsnPanelMainString         = "18304104 ZK"
    resPanel@gsnPanelMainFontHeightF    = 0.03
    
    ; 输出展板
    gsn_panel(wks,(/plot_a,plot_b/),(/2,1/),resPanel)
        
    ; 转换图像类型

    delete(wks)
    
    system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
    system("/bin/rm -f " + wksName + ".eps")  

	; =========================================================================
	; 去掉空白空间

    
    

end