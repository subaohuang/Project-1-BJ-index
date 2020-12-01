begin

	; 读取 .nc 文件
    f     =addfile("tp.assignment3.nc", "r")
    u     =addfile("uwnd.assignment3.nc", "r")
    v     =addfile("vwnd.assignment3.nc", "r")
    ; print(f)   ;这个地方只有用print才可以显示具体信息，用printvarsummary是不可以的
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	time  = f->time
    timeu = v->time
    timev = v->time
    ; print(time)

	; 2018年6月27日
	timereadid=cd_calendar(timeu, -5)
    year      =timereadid(:,0)
    month     =timereadid(:,1)
    day       =timereadid(:,2)
	timeid0   =ind(year.eq.2018.and.month.eq.6.and.day.eq.27)
    ; print(timeid0)     ;此时输出的是514
    ; print(timereadid(1734,:))    ;进行检验，输出的是2018-6-27-11:30,所以第1734个是对应的时间点
    timeid1   =ind(year.eq.2018.and.month.eq.6.and.day.eq.27)
    ; print(timeid1)
    ; print(timereadid(1734,:))  ;进行检验，输出的是2018-6-27-11:30,所以第1734个是对应的时间点
	; 2018年6月15日至7月15日
	timeid2   =ind(year.eq.2018.and.month.eq.6.and.(day.ge.15.and.day.le.30))
    ; print(timeid2)  ;输出的是502-517
    ;进行检验
    ; print(timereadid(502,:))
    ; print(timereadid(517,:))

    timeid3   =ind(year.eq.2018.and.month.eq.7.and.(day.ge.1.and.day.le.15))
	; print(timeid3)  ;输出的是518-532
    ;进行检验
    ; print(timereadid(518,:))
    ; print(timereadid(532,:))
    ;综上，此时的时间序列是502-532

	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	pp         =f->tp
    pp         =pp*1000
    printVarSummary(pp)
    ppp0       =pp(514,{15:45},{100:140})
	printMinMax(ppp0, 0)    ;其范围大小是0-8
    printVarSummary(ppp0)
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	ppp1       =pp(time|502:532,{latitude|35},{longitude|127.5});注意这个地方第一维那里不用括号
    ; printVarSummary(ppp)   ;这个时候的单位是m，时间序列是31
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	; print(u)
 ;    print(v)
    u0       =u->uwnd
    uu       =u0(1734,{700},:,:)
    v0       =v->vwnd
    vv       =v0(1734,{700},:,:)
    ; printVarSummary(uu)  ;这个地方lat0-60；lon90-150；二维坐标均是25的大小
    ; printVarSummary(vv)
	printMinMax(uu, 0)
	
	
	; ========================================================================
	; 打开工作空间
	wks     =gsn_open_wks("pdf", "output")
	
	
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
	res@tiMainString                    = "The percipitation and wind of 700hPa in East Asia on Jun.27~S~th~N~.2018"
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
    res@gsnAddCyclic                    =False
	res@mpMinLatF                       =15
    res@mpMaxLatF                       =45
    res@mpMinLonF                       =100
    res@mpMaxLonF                       =140
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette                ="precip3_16lev"
	

	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnLevelSelectionMode         ="ManualLevels"
    res@cnMinLevelValF               =0
    res@cnMaxLevelValF               =165
    res@cnLevelSpacingF              =10
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF         =0.02
    res@tmYLLabelFontHeightF         =0.02
	res@tiMainFontHeightF            =0.018
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF               =0.015
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@gsnAddCyclic              = False
    rescn@gsnDraw                   =False
    rescn@gsnFrame                  =False
    rescn@gsnLeftString             =""
    rescn@gsnRightString            =""
	rescn@cnLevelSelectionMode      ="ExplicitLevels"
    rescn@cnLevels                  =30
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor               ="White"
	rescn@cnLineThicknessF          =3
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn            =False
    rescn@cnInfoLabelOn             =False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
    resvc@gsnAddCyclic              = False 
	resvc@gsnFrame                  =False
    resvc@gsnDraw                   =False
    resvc@gsnLeftString             =""
    resvc@gsnRightString            =""
	resvc@vcGlyphStyle              ="FillArrow"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF     =0.02

	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF              =0.01
	resvc@vcRefMagnitudeF           =5
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	; resvc@vcRefAnnoOn               =True
    ; resvc@vcRefAnnoParallelPosF     =0.5
    resvc@vcRefAnnoOrthogonalPosF   = -0.2
	
	resvc@vcRefAnnoString1On        =True
    resvc@vcRefAnnoString2On  =False
    resvc@vcRefAnnoString1          ="2m/s"
    resvc@vcRefAnnoBackgroundColor  ="White"
    resvc@vcRefAnnoPerimOn          =False

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
	resxy@tiMainString                  ="The percipitation from Jun.15~S~th~N~ to Jul.15~S~th~N~ in 2018 at 35~S~0~N~N & 127.5~S~0~N~E"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	
	
	
	; gsnXYBarChart
	
	
	
	; X 轴
	; tmXBMode 及相应属性

	
	
	; (可选) XY线及柱状的颜色
	
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
    plot                             =gsn_csm_contour_map_ce(wks, ppp0, res)

	
	; 降水=30 mm
	plotcn                           =gsn_csm_contour(wks, ppp0, rescn)
	overlay(plot, plotcn)
	; 700百帕风场
	plotvc                           =gsn_csm_vector(wks, uu, vv, resvc)
	overlay(plot, plotvc)
	; 添加标记
	
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
    draw(plot)
    frame(wks)	
	
	
	; =========================================================================
	; 去掉空白空间
	

end