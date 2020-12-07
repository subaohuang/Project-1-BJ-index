begin

	; 读取 .nc 文件
	f1=addfile("tp.assignment3.nc","r")
	f2=addfile("uwnd.assignment3.nc","r")
	f3=addfile("vwnd.assignment3.nc","r")

	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	time1=f1->time
	time2=f2->time
	time3=f3->time
	time1convert=cd_calendar(time1, 2)
	time2convert=cd_calendar(time2, 2)
	time3convert=cd_calendar(time3, 2)

	; 2018年6月27日
	time1ind=ind(time1convert.eq.20180627)
	time2ind=ind(time2convert.eq.20180627)
	time3ind=ind(time3convert.eq.20180627)
	
	; 2018年6月15日至7月15日
	time1ind_=ind(time1convert.ge.20180615.and.time1convert.le.20180715)
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp=f1->tp
	tp=(/tp*1000/)
	tp@units :="mm"
	
	tp1=tp(time1ind,:,:)
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2=tp(time1ind_,{35},{127.5})
	
	; 提取 2018年6月27日 700百帕的风场
	uwnd=f2->uwnd(time2ind,{700},:,:)
	vwnd=f3->vwnd(time3ind,{700},:,:)
	
	; ========================================================================
	; 打开工作空间
	pltName="18304091.assignment3.q1"
	wks=gsn_open_wks("eps",pltName)
	
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
	res@gsnCenterString="(a)Total precipitation(mm) and 700-hPa wind(ms~S~-1~N~) in East Asia(20180627)"
	
	
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
	res@cnFillPalette="MPL_GnBu"
	
	
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
	res@tmXBLabelFontHeightF=0.024
	res@tmYLLabelFontHeightF=0.024
	res@gsnCenterStringFontHeightF=0.024
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	; 调节色板
	res@lbOrientation="Vertical"
	res@lbBoxEndCapStyle="TriangleBothEnds"
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode="ExplicitLevels"
	rescn@cnLevels=30

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor="black"
	rescn@cnLineThicknessF=3.
	
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
	resvc@vcLineArrowThicknessF=1.5
	resvc@vcMinDistanceF=0.03
	
	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF=0.1
	resvc@vcRefMagnitudeF=18
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF=1.27
	resvc@vcRefAnnoOrthogonalPosF=0.03
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	resvc@vcRefAnnoFontHeightF=0.01
	
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
	resxy@tiMainString="(b)Total precipitation(mm) on 20180615-20180715 (35~S~o~N~N 127.5~S~o~N~E)"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="Time"
	resxy@tiYAxisString="Total precipitation(mm)"
	
	; gsnXYBarChart
	resxy@gsnXYBarChart=True
	
	; X 轴
	; tmXBMode 及相应属性                         
	t=cd_calendar(tp2&time, 2)
	twant=(/20180616,20180621,20180626,20180701,20180706,20180711/)
	twantind=get1Dindex(t, twant)

	resxy@tmXBMode="Explicit"
	resxy@tmXBValues=tp2&time(twantind)
	resxy@tmXBLabels=(/"0616","0621","0626","0701","0706","0711"/)
	resxy@tmXBMinorValues=tp2&time(:)
	
	; (可选) XY线及柱状的颜色
	resxy@gsnYRefLine           = 0.0 
	resxy@gsnAboveYRefLineColor = "blue"             
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF=0.024
	resxy@tiXAxisFontHeightF=0.024
	resxy@tiYAxisFontHeightF=0.024
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot=gsn_csm_contour_map(wks, tp1, res)
	
	; 降水=30 mm
	plotcn 	= gsn_csm_contour(wks,tp1,rescn)
	overlay(plot,plotcn)

	; 700百帕风场
	plotvc=gsn_csm_vector(wks, uwnd,vwnd,resvc)
	overlay(plot,plotvc)
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark=True
	resMark@gsMarkerIndex=0
	resMark@gsMarkerColor="red"
	resMark@gsMarkerSizeF=0.024
	resMark@gsMarkerThicknessF=3.
	plotM=gsn_add_polymarker(wks, plot,127.5,35,resMark)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plotxy=gsn_csm_xy(wks,tp2&time,tp2,resxy)

	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel=True
	resPanel@gsnPanelMainString="Yang Tingting 18304091"
	resPanel@gsnPanelMainFontHeightF=0.025
	; resPanel@gsnPanelXF=0.1
	; resPanel@gsnPanelYF=0.9
	resPanel@gsnPanelYWhiteSpacePercent = 5
	gsn_panel(wks,(/plot,plotxy/),(/2,1/),resPanel)
	
	; =========================================================================
	; 去掉空白空间
	delete (wks) 			
	system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
	system("rm " + pltName + ".eps -f") 

end