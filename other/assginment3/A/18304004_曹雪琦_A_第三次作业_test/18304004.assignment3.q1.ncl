begin

	; 读取 .nc 文件
    f = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc","r")

	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
    time       = f->time
    timeARR = cd_calendar(time,-2)
	; 2018年6月27日
    timeInd1    = ind(timeARR.eq.20180627)
	; 2018年6月15日至7月15日
	timeInd2   = ind(timeARR.ge.20180615.and.timeARR.le.20180715)

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp        = f->tp
    tph       = 1000*tp
    tph@units ="mm"
    latitude  = f->latitude
    longitude = f->longitude
    tph       = 1000*tp
    tph@units ="mm"
    tph!0     ="time"
    tph!1     ="latitude"
    tph!2     ="longitude"
    tph&longitude =longitude
    tph&time      =time
    tph&latitude  =latitude
    tp1       =tph(time|timeInd1,latitude|:,longitude|:) 
    
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2       =tph(time|timeInd2,{latitude|35},{longitude|127.5})
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	ufile     = "/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc"
    vfile     = "/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc"
    uf         = addfile(ufile,"r")
    vf         = addfile(vfile,"r")
    time1 = uf->time
    timeARR1 = cd_calendar(time1,-2)
    timeInd3    = ind(timeARR1.eq.20180627)
    u    = uf->uwnd(time|timeInd3,{level|700},lat|:,lon|:)
    v    = vf->vwnd(time|timeInd3,{level|700},lat|:,lon|:)
	
	
	; ========================================================================
	; 打开工作空间
	wks=gsn_open_wks("png", "18304004.assignment3.q1")
	
	
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
	res@tiMainString                    = "20180627 Precipitation and 700hPa wind field in East Asia"
	
	
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
	res@cnFillPalette                   ="BlueWhiteOrangeRed" 
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode            ="ExplicitLevels"
    res@cnLevels                        = ispan(0,100,10) + "mm"
	
	
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
    rescn@cnLevelSelectionMode          ="ExplicitLevels"
	rescn@cnLevels                      = 30
    rescn@cnLineThicknessF              = 3.5
    rescn@cnLineColor                   = "red"
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn        = False
    rescn@cnInfoLabelOn         = False

	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle             = "FillArrow"
    
    
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcMinDistanceF        = 0.05
    resvc@vcLineArrowThicknessF = 0.2
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF          = 0.05
    resvc@vcRefMagnitudeF       = 10
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoFontHeightF             = 0.01
    resvc@vcRefAnnoOrthogonalPosF          = -0.2
    resvc@vcRefAnnoParallelPosF            = 1.2
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    ; resvc@vcRefAnnoOrthogonalPosF=0.2
	
	
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
	resxy@tiMainString                 = "2018.6.15-2018.7.15 precipitation"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
    resxy@tiXAxisString                ="time"
	resxy@tiYAxisString                ="precipitation"
	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart                =True
	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode                      = "Explicit"
    resxy@tmXBLabelsOn                  = True
    resxy@tmXBValues                    = ispan(1, 30, 5)
    resxy@tmXBLabels                    = (/"0616","0621","0626","0701","0706","0711"/)
    resxy@tmXBMinorOn                   = True
    resxy@tmXBMinorValues               = ispan(1, 31, 1)
	
	
	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors           = "blue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    resMark                    =True
	resMark@gsMarkerIndex      = 6
    resMark@gsMarkerSizeF      = 0.02
    resMark@gsMarkerColor      = "red"
    resMark@gsMarkerThicknessF =5.
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)l  
	plot              = new(2,graphic)
    plot(0)           = gsn_csm_contour_map_ce(wks, tp1, res)
    plotOverlay1      = gsn_csm_contour(wks, tp1, rescn)
    overlay(plot(0), plotOverlay1)
    plotOverlay2      =gsn_csm_vector(wks, u, v, resvc)
    overlay(plot(0), plotOverlay2)
    plotM             =gsn_add_polymarker(wks, plotOverlay2, 127, 35, resMark)

    ; --------------------------------------------------------------
    ; 9) 输出展版

    

	; 降水地图
	
	
	; 降水=30 mm
	
	
	; 700百帕风场
	
	
	; 添加标记
	
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	x                  = ispan(1,31,1)
    y                  = tp2
    plot(1)            = gsn_csm_xy(wks, x, y, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel                            = True
    resPanel@gsnPanelMainString         = "18304004caoxueqi"
    resPanel@gsnPanelXWhiteSpacePercent = 2  

    gsn_panel(wks, plot, (/1,2/), resPanel)
	
	
	
	; =========================================================================
	; ; 去掉空白空间
	; delete(wks)
 ;    system("convert -geometry 3000x3000 -density 600x600 -trim " + "18304004.assignment3.q1" + ".eps " + "18304004.assignment3.q1" + ".png")   
 ;    system("rm " + "18304004.assignment3.q1" + ".eps -f")


end