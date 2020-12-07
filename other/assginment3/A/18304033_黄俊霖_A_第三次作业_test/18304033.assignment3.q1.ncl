load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin

	; 读取 .nc 文件
    f_tp       = addfile("tp.assignment3.nc","r")
    f_uwnd     = addfile("uwnd.assignment3.nc","r")
    f_vwnd     = addfile("vwnd.assignment3.nc","r")
    tp     = f_tp->tp
	uwnd     = f_uwnd->uwnd
    vwnd     = f_vwnd->vwnd
    time_1=f_tp->time
    time_2=f_uwnd->time
    time_3=f_vwnd->time
	; printVarSummary(tp)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	time_1_cd=cd_calendar(time_1, -2)
	time_2_cd=cd_calendar(time_2, -2)
	time_3_cd=cd_calendar(time_3, -2)
    time1Ind=ind(time_1_cd.eq.20180627)
    time2Ind=ind(time_2_cd.eq.20180627)
    time3Ind=ind(time_3_cd.eq.20180627)
	; 2018年6月15日至7月15日
	time1Ind_1=ind(time_1_cd.ge.20180615.and.time_1_cd.le.20180715)
	time2Ind_1=ind(time_2_cd.ge.20180615.and.time_2_cd.le.20180715)
	time3Ind_1=ind(time_3_cd.ge.20180615.and.time_3_cd.le.20180715)
   ;========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp_rain=tp(time1Ind,:,:)
    tp_rain=(/tp_rain*1000./)
	;printVarSummary(tp_rain)
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
	tp_rain_1=tp(time1Ind_1,{35},{127.5})
	tp_rain_1=1000*tp_rain_1
    time_tp=time_1_cd(time1Ind_1)
    print(time_tp)
    
	;printVarSummary(tp_rain_1)
    
    ;print(time_tp)
	; 提取 2018年6月27日 700百帕的风场
	uwnd_700=uwnd(time1Ind,{700},:,:)
    vwnd_700=vwnd(time1Ind,{700},:,:)
	
	
	
	; ========================================================================
	; 打开工作空间
	wksName="18304033.assignment3.q1"
	wks=gsn_open_wks("eps", "18304033.assignment3.q1")
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	res@gsnLeftString 					= ""
	res@gsnRightString 					= ""
	res@gsnAddCyclic    =False
	resxy	 							= res
	
	res@gsnAddCyclic					= False
	rescn								= res
	resvc 								= res
				
	res@vpWidthF 						= 0.75
	res@vpHeightF			            = 0.5
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	
	res@tiMainString="(a) 20180627 precipitation and wind(700hPa)"
	
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
	res@cnFillPalette="MPL_Blues"
	res@lbBoxEndCapStyle="TriangleBothEnds"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnLevelSelectionMode    = "ExplicitLevels"
    res@cnLevels=(/0,5,10,30,50,80,120/)
    ; res@cnMinLevelValF          = 0.
    ; res@cnMaxLevelValF          = 200.
    ; res@cnLevelSpacingF         = 20.
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
	res@tiMainFontHeightF=0.024
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	rescn@cnFillOn=False
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@cnLevelSelectionMode         = "ExplicitLevels"
    rescn@cnLevels                     = 30
   
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor="Black"
	rescn@cnLineThicknessF=3.
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	
	rescn@cnLineLabelsOn=False
    rescn@cnInfoLabelOn=False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle      = "LineArrow"
    resvc@vcRefLengthF      = 0.08
    resvc@vcRefMagnitudeF   = 20
    resvc@vcMinDistanceF    = 0.04
    resvc@vcLineArrowThicknessF =1.5
    resvc@vcRefAnnoOn=True
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF=0.15
    resvc@vcRefAnnoSide="Right"
	resvc@vcRefAnnoOrthogonalPosF=0.015
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF

	
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.75
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	
	resxy@tiMainString="(b)2018.6.15-7.15 35~S~o~N~N,127.5~S~o~N~E precipitation"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="time(date)"
    resxy@tiYAxisString="precipitation(mm)"
	
    resxy@tiMainFontHeightF=0.024
    
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart=True
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode="Explicit"

    resxy@tmXBValues = ispan(1, sizeof(time_tp)/4 - 5, 5)

    time_labels = time_tp((/1, 6, 11, 16, 21, 26/))
	resxy@tmXBLabelFontHeightF=0.018
    resxy@tmXBLabelAngleF=60.
    ; resxy@tmXBLabelStride=1
    resxy@tmXBLabelDeltaF=2.
    resxy@tmXBMinorValues = ispan(0, 30, 1)
	; (可选) XY线及柱状的颜色
	resxy@gsnYRefLine           = 0.
    resxy@gsnXYBarChart  = True
    resxy@gsnAboveYRefLineColor =  "blue"      
    resxy@tmXBLabels = time_labels
	
    lab_x = new(sizeof(time_tp)/4, "integer")
 
    lab_x = ispan(0, sizeof(time_tp)/4 - 1, 1)
 ;    print(lab_x)
	; print(tp_rain_1)
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    resMark=True
	resMark@gsMarkerIndex=3
    resMark@gsMarkerColor="red"
    resMark@gsMarkerSizeF=0.02
	resMark@gsMarkerThicknessF=3.5
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot=gsn_csm_contour_map_ce(wks, tp_rain, res)
	
	; 降水=30 mm
	plotline=gsn_csm_contour(wks, tp_rain, rescn)
	overlay(plot, plotline)
	; 700百帕风场
	plotwind=gsn_csm_vector(wks, uwnd_700, vwnd_700, resvc)
	overlay(plot, plotwind)
	; 添加标记
    plotM=gsn_add_polymarker(wks, plot, 127.5, 35., resMark)

	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	; resTime                    = True
 ;    resTime@ttmFormat       = "%c %d"

 ;    ; resTime@ttmAxis         = "XB"
 ;    resTime@ttmMajorStride  =1
 ;    ; resTime@ttmMinorStride  =0
    
 ;    time_axis_labels(time_tp,resxy,resTime)
	plotB=gsn_csm_xy(wks,lab_x,tp_rain_1, resxy)
	

	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	plot_1   = new(2,graphic)
    plot_1(0)=plot
    plot_1(1)=plotB
    ; 展板设定
    
    resPanel                            = True


    ; 标题
    resPanel@gsnPanelMainString         = "18304033-hjl"
    resPanel@gsnPanelMainFontHeightF    = 0.024

    gsn_panel(wks,plot_1,(/2,1/), resPanel)
	draw(plot_1)
    frame(wks)
	
	; =========================================================================
	; 去掉空白空间
	delete(wks)
    
    system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
    system("/bin/rm -f " + wksName + ".eps")

end