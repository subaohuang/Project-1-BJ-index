begin
;load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_inv_string.ncl"

	; 读取 .nc 文件
	; precip=addfile("D:\assignment3\tp.assignment3.nc", "r")
	; uwind=addfile("D:\assignment3\uwnd.assignment3.nc", "r")
	; vwind=addfile("D:\assignment3\vwnd.assignment3.nc", "r")

	precip				= addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc","r")
	uwind			= addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc","r")
	vwind 				= addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc","r")



	time_precip=precip->time
	time_u=uwind->time
	time_v=vwind->time
	fi_precip=cd_calendar(time_precip, -2)
	fi_u=cd_calendar(time_u, -2)
	fi_v=cd_calendar(time_v, -2)
	


	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	fi_tp_select1=ind(fi_precip .eq. 20180627)
	fi_u_select1=ind(fi_u .eq. 20180627)
	fi_v_select1=ind(fi_v .eq. 20180627)
	
	
	; 2018年6月15日至7月15日
	
	fi_tp_select2=ind(fi_precip .ge. 20180615 .and. fi_precip .le. 20180715)
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1=precip->tp(fi_tp_select1,:,:)
	tp1_use=tp1*1000.0		
	copy_VarCoords(tp1,tp1_use)
    printVarSummary(tp1_use)
    tp1_use@units="mm"
    printMinMax(tp1_use,0)

	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp_lon=precip->longitude
	tp_lat=precip->latitude
	tp_lon_sel=ind(tp_lat.eq.35)
	tp_lat_sel=ind(tp_lon.eq.127.5)
	tp2=precip->tp(fi_tp_select2,tp_lon_sel,tp_lat_sel)
	time_tp=precip->time(fi_tp_select2)

	;time_use=cd_calendar(time_tp, -2)
	time_six=cd_calendar(time_tp, 0)
	month=tointeger(time_six(:,1))
	day=tointeger(time_six(:,2))
	time_use=sprinti("%0.2i",month)+sprinti("%0.2i",day)
	;print(time_use)

	;print(time_use(1:30:5))
	tp2_use=tp2*1000.0
	copy_VarCoords(tp2, tp2_use)
	tp2_use@units="mm"

	printMinMax(tp2_use,0)
	
	
	; 提取 2018年6月27日 700百帕的风场
	
	uw=uwind->uwnd(fi_u_select1,1,:,:)
    vw=vwind->vwnd(fi_v_select1,1,:,:)
	printMinMax(uw, 0)
	printMinMax(vw, 0)
	
	; ========================================================================
	; 打开工作空间
	pltName="18304007.assignment3.q1.test"
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
	res@tiMainString 			= " The Figure of Precipitation(contour) the wind (vector) "	 
	res@gsnCenterString 		= " East Asia in 2018.6.27"
	
	; 地图范围

	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF=15
	res@mpMaxLatF=45
	res@mpMinLonF=100
	res@mpMaxLonF=145
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"

	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette="precip3_16lev"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode 	="ExplicitLevels"
	res@cnLevels 				=ispan(0,200,10)

	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF = 0.02
	res@tmYLLabelFontHeightF = 0.02
	res@tiMainFontHeightF 	 = 0.02
	res@gsnCenterStringFontHeightF =0.018

	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	res@tiMainOffsetYF 		 = 0.012
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode="ExplicitLevels"
	rescn@cnLevels=(/30/)

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor   	= "Red"
	rescn@cnLineThicknessF  = 5.0
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn 	= False
	rescn@cnInfoLabelOn 	= False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle

	resvc@vcGlyphStyle 		="CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF 	= 1.5
	resvc@vcMinDistanceF 			= 0.03
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF 	 		= 0.02
	resvc@vcRefMagnitudeF       = 8
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF = 1.15
	resvc@vcRefAnnoOrthogonalPosF = -0.25
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	resvc@vcRefAnnoFontHeightF =0.01
	
	
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
	resxy@tiMainString   	= "The Time Series of Precipitation at 35~S~o~N~N,127.5~S~o~N~E "
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString 	= "Time(unit=day)"
	resxy@tiYAxisString 	= "Precipitation(unit=mm)"

	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart 	=	True
	
	
	; X 轴
	;设定x轴:tmXBMode
	loc=ispan(1,31,1)
	resxy@tmXBMode 		="Explicit"
	resxy@tmXBValues 	=loc(1::5)
	resxy@tmXBLabels 	=time_use(1::5)
	resxy@tmXBMinorOn 	=True
	resxy@tmXBMinorValues = loc 

	
	
	
	; (可选) XY线及柱状的颜色
	
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF =0.02

	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resmark					=	True
	resmark@gsMarkerIndex 	=	1
	resmark@gsMarkerColor 	=	"Yellow"
	resmark@gsMarkerSizeF 	=  	0.04
	resText				= True
	resText@txFontHeightF 		= 0.02
	resText@txFontColor 		= "Yellow"	
	

	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plotmap=gsn_csm_contour_map_ce(wks, tp1_use, res)

	; 降水=30 mm
	plotcn=gsn_csm_contour(wks, tp1_use, rescn)
	overlay(plotmap,plotcn)
	
	
	; 700百帕风场
	plotvc=gsn_csm_vector(wks, uw, vw, resvc)
	overlay(plotcn,plotvc)

	; 添加标记
	
	plotM=gsn_add_polymarker(wks, plotmap, 127.5, 35.0, resmark)
	plotText		= gsn_add_text(wks,plotmap,"35~S~o~N~N,127.5~S~o~N~E",133.5,35.0,resText)
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plotxy=gsn_csm_xy(wks,ispan(1,31,1),tp2_use, resxy)

	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel = True
	resPanel@gsnPanelMainString="18304007 ChenAoxuan"
	;resPanel@gsnPanelXF = (/0,0/)
	;resPanel@gsnPanelYF = (/0.3,0.8/)
	plot=(/plotmap,plotxy/)
	gsn_panel(wks,plot,(/2,1/), resPanel)

	; =========================================================================
	; 去掉空白空间
	delete (wks) 			; 必须先关闭 workstation
	system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
	system("rm " + pltName + ".eps -f") 

end