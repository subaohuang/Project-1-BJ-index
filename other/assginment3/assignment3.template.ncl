begin

	; 读取 .nc 文件
	f1          = addfile("tp.assignment3.nc", "r")
	f2          = addfile("uwnd.assignment3.nc", "r")
    f3          = addfile("vwnd.assignment3.nc", "r")
    tp          = f1->tp
    uwnd        = f2->uwnd
    vwnd        = f3->vwnd

	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time1       = cd_calendar(f1->time, 2)
    time1Ind1   = ind(time1.eq.20180627)

    time2       = cd_calendar(f2->time, 2)
	time2Ind    = ind(time2.eq.20180627)

	; 2018年6月15日至7月15日
	time1Ind2   = ind(time1.ge.20180615.and.time1.le.20180715)
	timevalue   = time1(time1Ind2)
    timelabel   = timevalue(ispan(1,26,5))
    ; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	data1        =  tp(time1Ind1,:,:)

    data_1       = data1*1000
    copy_VarMeta(data1, data_1)
    data_1@units ="mm"


	; ; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; ; 单位需要由 m 转为 mm
	data2       = tp(time1Ind2,{35.},{127.5})
    data_2      = data2*1000
    copy_VarMeta(data2, data_2)
    data_2@units ="mm"

	; ; 提取 2018年6月27日 700百帕的风场
	data3       = uwnd(time2Ind,{700},:,:)
	data4       = vwnd(time2Ind,{700},:,:)


	; ; ========================================================================
	; ; 打开工作空间
	wks         =gsn_open_wks("eps", "18304003.assignment3")
	
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	res@gsnLeftString 				    = " "
	res@gsnRightString 				    = " "
	
	resxy	 							= res
	
	res@gsnAddCyclic					= False
	rescn								= res
	resvc 								= res
				
	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	res@tmXTOn                          = False
    res@tmYROn                          = False
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; ; res: 图 (a) 底图
	
	; 标题
	res@tiMainString                    = "2018.6.27_precipitation_700wind"
	
	
	; ; 地图范围
	res@mpMinLatF                       = 15.
    res@mpMaxLatF                       = 45.
    res@mpMinLonF                       = 100.
    res@mpMaxLonF                       = 140.
	
	; ; 高分辨率地图
	res@mpDataBaseVersion               = "HighRes"
	res@mpDataResolution               = "FinestResolution"
    ; res@mpDataSetName                   = "Earth..4"
    res@mpOutlineBoundarySets           = "National"
	
	; ; 色板
	res@cnFillPalette                   = "BlueWhiteOrangeRed"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode            = "ExplicitLevels"
    res@cnLevels                        = ispan(0,80,8) + "mm"
	
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
	rescn@cnLevelSelectionMode             ="ExplicitLevels" 
    rescn@cnLevels                          = (/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor                      = "Green4"
    rescn@cnLineThicknessF                 = 4.0
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn                   = False
	rescn@cnInfoLabelOn                    = False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	resvc@vcGlyphStyle                     = "CurlyVector"
	resvc@vcMinDistanceF                   = 0.02
    
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	resvc@vcRefLengthF                     = 0.1
    resvc@vcRefMagnitudeF                  = 15
	
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoFontHeightF             = 0.01
    resvc@vcRefAnnoOrthogonalPosF          = -0.2
	resvc@vcRefAnnoParallelPosF            = 1.2
	      
	; (可选) 参考矢量信息的字体大小
	; vcRefAnnoFontHeightF
	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	resxy@vpWidthF 						= 0.6
	resxy@vpHeightF 					= 0.5
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	resxy@trYMinF                       = 0.
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString                  = "2018.6.15-2018.7.15_precipitation"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString                 = "time"
	resxy@tiYAxisString                 = "precipitation"
	resxy@gsnXYBarChart                 = True
	
	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode                      = "Explicit"
    resxy@tmXBLabelsOn                  = True
    ; resxy@tmXBValues                    = timelabel
    resxy@tmXBValues      = ispan(1, 30, 5)
  
    resxy@tmXBLabels                    = (/"0616","0621","0626","0701","0706","0711"/)
    resxy@tmXBMinorOn                   = True
	; resxy@tmXBMinorValues               = f1->time(time1Ind2)
    resxy@tmXBMinorValues = ispan(1, 31, 1)
    

	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors           = "skyblue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    resMark                               = True 
	resMark@gsMarkerIndex                 = 5
	resMark@gsMarkerColor                 = "red"
	resMark@gsMarkerSizeF                 = 0.025
    resMark@gsMarkerThicknessF            = 3.
	; =========================================================================
	; 画图
	plot               = new(2,graphic)
	; 图 (a)
	
	; 降水地图
	plot(0)            = gsn_csm_contour_map_ce(wks, data_1, res)    
	
	; 降水=30 mm
	plotOverlay_1      = gsn_csm_contour(wks, data_1, rescn)
    overlay(plot(0), plotOverlay_1)

	; 700百帕风场
    plotOverlay_2      = gsn_csm_vector(wks, data3, data4, resvc)
    overlay(plot(0), plotOverlay_2)

	; 添加标记
	plotM              = gsn_add_polymarker(wks,plotOverlay_2,127.5,35.,resMark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	; x                  = f1->time(time1Ind2)
    x = ispan(1, 31, 1)
    y                  = data_2

    ; print(timelabel)
    ; print(x)

    ; exit()
	plot(1)            = gsn_csm_xy(wks, x, y, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
    resPanel                            = True
    resPanel@gsnPanelMainString         = "18304003_CaiYaoting"
    resPanel@gsnPanelXWhiteSpacePercent = 2.   
    ; respanel@tiMainFontHeightF          = 0.03
    ; respanel@gsnPanelXF                 = 0.1
    ; respanel@gsnPanelYF                 = 0.02
	gsn_panel(wks, plot, (/1,2/), resPanel)
	
	
	
	; =========================================================================
	; 去掉空白空间
	delete (wks)            ; 必须先关闭 workstation
    system("convert -geometry 3000x3000 -density 600x600 -trim " + "18304003.assignment3" + ".eps " + "18304003.assignment3" + ".png")   
    system("rm " + "18304003.assignment3" + ".eps -f")

end