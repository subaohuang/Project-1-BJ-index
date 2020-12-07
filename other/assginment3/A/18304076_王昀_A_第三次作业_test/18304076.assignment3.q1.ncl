begin

	; 读取 .nc 文件
    path         = "/cygdrive/d/18304076/"
	tpf 		 = addfile("tp.assignment3.nc","r")
	uf 		     = addfile("uwnd.assignment3.nc","r")
	vf 		     = addfile("vwnd.assignment3.nc","r")

	
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
    time            = tpf->time
    tYYYYMMDD       = cd_calendar(time, -2)
    time_ind1       = ind(tYYYYMMDD.eq.20180627) 

	
	
	; 2018年6月15日至7月15日
    time_ind2       = ind(tYYYYMMDD.ge.20180615.and.tYYYYMMDD.le.20180715) 
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1       = tpf->tp(time|time_ind1,level|:,lat|:,{lon|:})
	tp1_sub   = tp1*1000
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2       = tpf->tp(time|time_ind2,level|:,lat|35,{lon|127.5})
	tp2_sub   = tp2*1000
	
	
	; 提取 2018年6月27日 700百帕的风场
	uwnd 	= uf->uwnd(time|time_ind1,{level|700},lat|:,lon|:)
	vwnd 	= vf->vwnd(time|time_ind1,{level|700},lat|:,lon|:)
	
	print(time_ind1)
	print(time_ind2)
	print(tp1_sub)
	print(tp2_sub)
	print(uwnd)
	print(vwnd)
	
	; ========================================================================
	; 打开工作空间
	wks  = gsn_open_wks("png","18304076.assignment3.q1")
	
	
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
	res@tiMainString            = "2018.6.27东亚降水(单位mm)及700百帕风场(单位ms~S~-1~N~)"            
                        

	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF               = 15.
    res@mpMaxLatF               = 45.
    res@mpMinLonF               = 100.
    res@mpMaxLonF               = 140.
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion       = "HighRes"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette           = "MPL_BrBG"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode     = "ManualLevels"   ;??????
    res@cnMaxLevelValF           = 60.  
    res@cnMinLevelValF           = 0.
    res@cnLevelSpacingF          = 3.
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF    = 0.025                    ;设置X轴坐标字体高度，默认值0.02
    res@tmYLLabelFontHeightF    = 0.025                    ;设置Y轴坐标字体高度，默认值0.02
	
	res@tiMainFontHeightF       = 0.034 
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode 		= "ExplicitLevels"
    rescn@cnLevels 			        = 30.

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor			= "Black"
    rescn@cnLineThicknessF		= 3.5
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn                 = False
	rescn@cnInfoLabelOn                  = False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle        = "CurlyVector" 
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	res@vcLineArrowThicknessF           = 0.9
	res@vcMinDistanceF                  = 0.025
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefMagnitudeF               = 20.                          
    resvc@vcRefLengthF                  = 0.1 
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF         = 1.0  ;默认值：1.0
	resvc@vcRefAnnoOrthogonalPosF       = 0.02 ;默认值：0.02
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	res@vcRefAnnoFontHeightF            =0.02

	
	
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
	resxy@gsnXYBarChart           = True      ; 画柱状图
	resxy@tiMainString            = "2018.6.15-2018.7.15北纬35度，东经127.5度降水的时间序列"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString           = "时间"                                      ;x轴标题
    resxy@tiYAxisString           = "降水量(mm)"                                       ;y轴标题
    
	
	
	
	; gsnXYBarChart
	res@gsnXYBarChart             = True
	
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode                = "Manual"
    resxy@tmXBTickStartF          = 20180616
    resxy@tmXBTickEndF            = 20180711
    resxy@tmXBTickSpacingF        = 5

	
	
	; (可选) XY线及柱状的颜色
	
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark			        	= True
    resMark@gsMarkerIndex 		= 5
    resMark@gsMarkerThicknessF 	= 3.
    resMark@gsMarkerColor 		= "Blue"
    resMark@gsMarkerSizeF 		= 0.024

	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot        = gsn_csm_contour_map_ce(wks,tp1_sub,res)
	
	; 降水=30 mm
	plotcn   	= gsn_csm_contour(wks,tp1_sub,rescn)
	
	; 700百帕风场
	plotvc      = gsn_csm_vector_map_ce(wks,uwnd,vwnd,resvc)

	
	; 添加标记
	plotM		= gsn_add_polymarker(wks,plot,127.5,30.,resMark)

	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plotxy      = gsn_csm_xy(wks,time_ind2,tp2_sub,resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel 							= True
	resPanel@gsnPanelBottom				= 0.15
	resPanel@gsnPanelYWhiteSpacePercent = 5.
	
	
	
	; =========================================================================
	; 去掉空白空间
	

end