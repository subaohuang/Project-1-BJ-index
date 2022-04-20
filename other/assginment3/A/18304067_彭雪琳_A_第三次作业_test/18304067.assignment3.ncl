begin

	; 读取 .nc 文件
   ftp=addfile("tp.assignment3.nc","r")
   u=addfile("uwnd.assignment3.nc", "r")
   v=addfile("vwnd.assignment3.nc", "r")
	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	; vwnd.assignment3.nc
	; 2018年6月27日
	tp=ftp->tp
	time=ftp->time
	time1=cd_calendar(time, 2)
	timeInd1=ind(time1.eq.20180627)
	
	
	; 2018年6月15日至7月15日
	timeInd2=ind(time1.ge.20180615.and.time1.le.20180715)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1=tp(timeInd1,:,:)
	tp1=(/tp1*1000/)
	tp1@units="mm"
	
	printMinMax(tp1, 0)

	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2=tp(time|timeInd2,{latitude|35},{longitude|127.5})
	tp2=(/tp2*1000/)
	tp2@units="mm"
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	uwnd=u->uwnd
	vwnd=v->vwnd
	timeu=u->time
	timev=v->time
	time3=cd_calendar(timeu,2)
	time4=cd_calendar(timev,2)
	timeInd3=ind(time3.eq.20180627)
	timeInd4=ind(time4.eq.20180627)
	u1=uwnd(time|timeInd3,{level|700},lat|:,lon|:)
	v1=vwnd(time|timeInd4,{level|700},lat|:,lon|:)
	
	
	; ========================================================================
	; 打开工作空间
	wks=gsn_open_wks("eps", "18304067.assignment.q1")
	
	
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
	res@gsnCenterString     ="20180627_East Asia_Precipitation and wind_speed"
	res@gsnCenterStringFontHeightF   =0.04
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF=15
	res@mpMaxLatF=45
	res@mpMinLonF=100
    res@mpMaxLonF=140
    res@mpCenterLonF=120
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette="BlueWhiteOrangeRed"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode="ManualLevels"
	res@cnMinLevelValF=-10
	res@cnMaxLevelValF=110
	res@cnLevelSpacingF=10
	res@cnLineThicknessF =2.5
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF  =  0.02
	res@tmYLLabelFontHeightF  =  0.02
    res@tiMainFontHeightF	  =  0.02
    res@gsnCenterStringFontHeightF   =  0.02
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF        =  0.01
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	rescn@cnFillOn             =False
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode ="ExplicitLevels"
	rescn@cnLevels             = 30

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor          ="Red"
	rescn@cnLineThicknessF     = 4
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn       =False
	rescn@cnInfoLabelOn        =False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle         ="CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
    resvc@vcLineArrowThicknessF    = 2
    resvc@vcMinDistanceF           = 0.02
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF      =    0.1
	resvc@vcRefMagnitudeF   =    25
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF   =1.0
	resvc@vcRefAnnoOrthogonalPosF =0.02
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.522
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString                  ="20180615-20180715 Precipitation of 35~S~o~N~N,127.5~S~o~N~E"
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString  =  "date"
	resxy@tiYAxisString  =  "total precipitation(mm)"
	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart=True
	resxy@gsnXYBarChartBarWidth         = 0.8
    resxy@gsnXYBarChartColors           = "skyblue"
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode     =  "Explicit"
    resxy@tmXBValues   = (/1,6,11,16,21,26/)
    resxy@tmXBLabels   = (/"6.16","6.21","6.26","7.01","7.06","7.11"/)
    resxy@tmXBMinorValues  = ispan(0,30,1)
    
	
	
	; (可选) XY线及柱状的颜色
    resxy@xyMarkLineMode                = "MarkLines"
    resxy@xyMarker                      = 1
    resxy@xyMarkerSizeF                 = 0.01	
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF            = 0.02
    resxy@tiXAxisFontHeightF           = 0.02
    resxy@tiYAxisFontHeightF           = 0.02
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF  
	resMark                            = True
	resMark@gsMarkerIndex              =  16
	resMark@gsMarkerColor              =  "Black"
	resMark@gsMarkerSizeF              = 0.008
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plotbase =   gsn_csm_contour_map_ce(wks, tp1, res)
	
	; 降水=30 mm
	plotcn   = gsn_csm_contour(wks, tp1, rescn)
	
	; 700百帕风场
	plotOverlay   = gsn_csm_vector(wks, u1, v1, resvc)
	
	; 添加标记
	lon0 = (/127.5/)
	lat0 = (/35.0/)
	a = gsn_add_polymarker(wks, plotbase, lon0, lat0, resMark)

	overlay(plotbase, plotcn)
	overlay(plotbase, plotOverlay)


	; --------------------------------------------------------------------------
	; 图 (b)
	x=ispan(0,30,1)
	plotxy = gsn_csm_xy(wks, x, tp2, resxy)

	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)

	resPanel  =  True
    resPanel@gsnPanelMainString          = "18304067(Peng Xuelin)"
	plot = new(2,graphic)
    plot(0)=plotbase
    plot(1)=plotxy
	gsn_panel(wks,plot,(/1,2/),resPanel)
	; =========================================================================
	; 去掉空白空间

    delete (wks)            ; 必须先关闭 workstation
    system("convert -geometry 3000x3000 -density 600x600 -trim " + "18304067.assignment.q1" + ".eps " + "18304067.assignment.q1" + ".png")   
    system("rm " + "18304067.assignment.q1" + ".eps -f")

end