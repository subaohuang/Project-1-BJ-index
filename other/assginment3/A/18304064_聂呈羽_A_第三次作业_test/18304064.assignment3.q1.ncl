begin

	; 读取 .nc 文件
    f1=addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc", "r")
    f2=addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc", "r")
	f3=addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc", "r")
    print(f2)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	time1=f1->time
    time1:=cd_calendar(time1,-5)
    time2=f2->time
    time2:=cd_calendar(time2,-2)
    time3=f3->time
    time3:=cd_calendar(time3,-2)


	; 2018年6月27日
	y1=time1(:,0)
    m1=time1(:,1)
	d1=time1(:,2)	
    t1=ind(y1.eq.2018.and.m1.eq.06.and.d1.eq.27)
	; 2018年6月15日至7月15日
	t2=ind(y1.eq.2018.and.((m1.eq.06.and.d1.ge.15).or.(m1.eq.07.and.d1.le.15)))

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	rf1=f1->tp(time|t1,latitude|:,longitude|:)
	rf1=rf1*1000
    rf1@units="mm"
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	rf2=f1->tp(time|t2,{latitude|35},{longitude|127.5})
    rf2=rf2*1000
	rf2@units="mm"
	; 提取 2018年6月27日 700百帕的风场
    time4 = cd_calendar(f3->time,-5)
    y2=time4(:,0)
    m2=time4(:,1)
    d2=time4(:,2)   
    t2:=ind(y2.eq.2018.and.m2.eq.06.and.d2.eq.27)
    print(t1)
    print(t2)
	exit()
	u=f2->uwnd(time|t1,{level|700},lat|:,lon|:)
	v=f3->vwnd(time|t1,{level|700},lat|:,lon|:)
	
	; ========================================================================
	; 打开工作空间
	
	wks=gsn_open_wks("eps","18304064.assignment3.q1")
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	res@gsnLeftString 					= ""
	res@gsnRightString 					= ""
	
	resxy	 							= res
	res@gsnAddCyclic                   = False
	rescn								= res
	resvc 								= res
				
	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	
	res@gsnCenterString                    ="(a)2018 June 28th Rainfall and 700hPa U-V "
	
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
	
	res@cnFillPalette="WhiteBlue"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode="ManualLevels"
    res@cnMinLevelValF=0
    res@cnMaxLevelValF=120
    res@cnLevelSpacingF=10 
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF=0.02
	res@tmYLLabelFontHeightF=0.02
    res@tiMainFontHeightF=0.025
    res@gsnRightString="mm"
    res@gsnRightStringFontHeightF         = 0.02
    res@gsnRightStringOrthogonalPosF      = 0.03
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	res@gsnCenterStringOrthogonalPosF=0.085
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode="ExplicitLevels"
    rescn@cnLevels=30
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineThicknessF=1.5
    rescn@cnLineColor="navy"
	
	
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
    resvc@vcMinDistanceF     = 0.025
	
	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefMagnitudeF    = 5.
    resvc@vcRefLengthF       = 0.01
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF    =1.
    resvc@vcRefAnnoOrthogonalPosF      =-0.47
	
	
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
	resxy@tiMainString="(b)2018 Daily Rainfall June 15th - July 15th"
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiYAxisString="Rainfall (mm)"
	resxy@tiXAxisString="Day Time"
	
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart=True
	
	; X 轴
    x1=ispan(1, 30, 5)
    x=ispan(0,30, 1)
	; tmXBMode 及相应属性
    resxy@tmXBMode="Explicit"
    resxy@tmXBValues=x1
    resxy@tmXBLabels=(/"June 16th", "June 21st", "June 26th","July 01st", "July 06th","July 11th"/)
    resxy@tmXBMajorLengthF=0.02
    resxy@tmXBMinorOn=True
    resxy@tmXBMinorValues=x
    resxy@tmXBMinorLengthF=0.01
	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors="steelblue2"
    resxy@tmXBLabelFontHeightF=0.015
    resxy@tmYLLabelFontHeightF=0.015
    resxy@gsnXYBarChartBarWidth=0.75
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
    resxy@tiMainFontHeightF=0.025
	resxy@tiXAxisFontHeightF=0.02
    resxy@tiYAxisFontHeightF=0.02
	resxy@gsnRightString="mm"
    resxy@gsnRightStringFontHeightF         = 0.02
    resxy@gsnRightStringOrthogonalPosF      = 0.03
    resxy@tiMainOffsetYF=0
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark=True
    resMark@gsMarkerIndex=7
    resMark@gsMarkerThicknessF=2.5
    resMark@gsMarkerColor="red"
    resMark@gsMarkerSizeF=0.015
	; =========================================================================
	; 画图
	plot                           = new(2,graphic)
	; 图 (a)
	
	; 降水地图
	plot(0)=gsn_csm_contour_map_ce(wks, rf1, res)
	
	; 降水=30 mm
	plotcn=gsn_csm_contour(wks, rf1, rescn)
	overlay(plot(0), plotcn)

	; 700百帕风场
	plotvc=gsn_csm_vector(wks, u, v, resvc)
    overlay(plot(0), plotvc)

	; 添加标记
	plotM=gsn_add_polymarker(wks, plot(0),127.5,35, resMark)
    marktext                               = (/"35~S~o~N~N,127.5~S~o~N~E"/)
    rest=True
    rest@txFontHeightF                  = 0.016
    rest@txFontColor                    = "red"
    polytext                               = gsn_add_text(wks, plot(0),marktext,128,37, rest)

	; --------------------------------------------------------------------------
	; 图 (b)
	
	plot(1)=gsn_csm_xy(wks,x,rf2,resxy)
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel=True
    resPanel@gsnPanelMainString="18304064 niechengyu"
    resPanel@gsnPanelMainFontHeightF=0.025
    resPanel@gsnPanelXF=0.2
    resPanel@gsnPanelYF=0.7
	gsn_panel(wks, plot, (/2,1/), resPanel)

	; =========================================================================
	; 去掉空白空间
	delete(wks)
    system("convert -geometry 3000x3000 -density 600x600 -trim 18304064.assignment3.q1.eps 18304064.assignment3.q1.png")
    system("rm 18304051.assignment3.q1.eps -f")

end