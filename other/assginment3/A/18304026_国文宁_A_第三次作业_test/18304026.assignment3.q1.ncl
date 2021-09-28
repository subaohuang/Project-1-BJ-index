begin

	; 读取 .nc 文件
    f1  = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc", "r")
	f2  = addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc", "r")
    f3  = addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc", "r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time =f1->time
    cd_date = cd_calendar(time, -2) 
    timelnd1 =ind(cd_date .eq. 20180627)
	data1 =cd_date(timelnd1)
	
	; 2018年6月15日至7月15日
	timelnd2 = ind (cd_date .ge. 20180615 .and. cd_date .le. 20180715)
	data2 =cd_date(timelnd2)
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	p1 = f1->tp
    p11 =p1(timelnd1,:,:)
    p111 = 1000*p11
 
	
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	p2 = f1->tp
    p21 = p2(timelnd2,{35},{127.5})
    p211 = 1000*p21

	
	; 提取 2018年6月27日 700百帕的风场
	u=f2->uwnd(timelnd1,{700},:,:)
    v=f3->vwnd(timelnd1,{700},:,:)
	
	
	; ========================================================================
	; 打开工作空间
	wks =gsn_open_wks("eps", "18304026.assignment3.q1.test")
	
	
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
	res@tiMainString=  "18304026GuoWenning"
	
    res@gsnCenterString="Precipitation in the east of Asia on June 27 2018"
    
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF =15
	res@mpMaxLatF =45
    res@mpMinLonF =100
	 res@mpMaxLonF =140
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion = "MediumRes"

	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette = "NCV_jaisnd"

	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode    = "ManualLevels"
    res@cnMaxLevelValF          = 30.
    res@cnMinLevelValF          = 0.
    res@cnLevelSpacingF         = 3.
    res@cnLineThicknessF        = 0.5
    
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@gsnCenterStringFontHeightF     = 0.02
    res@tiMainFontHeightF           = 0.03
    res@tmXBLabelFontHeightF        = 0.015
    res@tmYLLabelFontHeightF        = 0.015
	

	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF              = -0.035
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@cnFillOn   =False
    rescn@cnLinesOn  =True
	rescn@cnLabelScalingMode ="ExplicitLevels"
	rescn@cnLevels = (/30/)
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor = "Pink"
    rescn@cnLineThicknessF = 0.023

	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnLineLabelsOn =False
    rescn@cnInfoLabelOn = False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle       = "LineArrow"   ;调整箭头类型
    resvc@vcRefMagnitudeF     = 20.                   ;调整箭头长度
    resvc@vcRefLengthF        = 0.1       
    resvc@vcMinDistanceF  = 0.025
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF = 0.015
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcMinAnnoParallelPosF = -0.37
    resvc@vcRefAnnoOrthogonalPosF = 0.45
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	resvc@vcRefAnnoFontHeightF = 0.012
	
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
	resxy@tiMainString="(Precipitation in 35~S~o~N~N 127.5~S~o~N~E) Jun 15 - Jul 15  "
	resxy@tiXAxisString = "date"
	resxy@tiYAxisString = "Precipitation"
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	
	
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart = True
	resxy@gsnYRefLine                = 0.
    resxy@gsnAboveYRefLineColor   = "Red"
    resxy@gsnBelowYRefLineColor   = "Blue"
    

	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBLabelsOn = True
    resxy@tmXBMode = "Manual"
	resxy@tmXBTickEndF     = 20180711
    resxy@tmXBTickStartF   = 20180616
    resxy@tmXBTickSpacingF = 5
	
	; (可选) XY线及柱状的颜色
	resxy@xyLineColor            = "Yellow"
    
    resxy@xyMonoDashPattern = True

	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
    resxy@tiMainFontHeightF           = 0.03
    resxy@tmXBLabelFontHeightF        = 0.015
    resxy@tmYLLabelFontHeightF        = 0.015
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    
	resMark                = True
    resMark@gsMarkerIndex       = 5
    resMark@gsMarkerThicknessF  = 3.
    resMark@gsMarkerColor       = "Green4"
    resMark@gsMarkerSizeF       = 0.024

	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	
	
	; 降水=30 mm
    plot       = gsn_csm_contour_map_ce(wks, p111, res)
    plotcn=gsn_csm_contour(wks, p111, rescn)
    overlay(plot, plotcn)

	 
	; 700百帕风场
	plotvc = gsn_csm_vector(wks, u, v, res)
	overlay(plot, plotvc)
	; 添加标记
	
	plotM         = gsn_add_polymarker(wks,plot,127.5,35,resMark)
	; --------------------------------------------------------------------------
	; 图 (b)
	
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
	res@gsnPanelXF = -0.034
    res@gsnPanelYF = 0.041
    
	plotxy=gsn_csm_xy(wks,data2, p211, resxy)
    
    ; resource of panel plot
    
    resPanel                                = True
    resPanel@gsnPanelLabelBar               = False  ; set a common label bar 
    resPanel@lbLabelFontHeightF             = 0.015 ; font size of the common label bar
    resPanel@gsnPanelXWhiteSpacePercent     = 5.    ; white space on x-axis between two plots in the panel plot
    resPanel@gsnPanelYWhiteSpacePercent     = 5.    ; white space on y-axis between two plots in the panel plot
    
    ; Use gsn_panel to draw a panel plot 
    
    plotpanel   = new(2,graphic)
    plotpanel(0)=plot
    plotpanel(1)=plotxy
    gsn_panel(wks,plotpanel,(/1,2/), resPanel)
   
	 draw(plotpanel)
     frame(wks)
    
	; =========================================================================
	; 去掉空白空间
    delete (wks)           
   system("convert -geometry 3000x3000 -density 600x600 -trim " +"18304026.assignment3.q1.test" + ".eps " + "18304026.assignment3.q1" + ".png")   
   system("rm " +"18304026.assignment3.q1.test" + ".eps -f") 
	

end