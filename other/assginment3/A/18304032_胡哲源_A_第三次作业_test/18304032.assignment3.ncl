begin

	; 读取 .nc 文件
    tp      = addfile("/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc","r")
	uwnd    = addfile("/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc","r")
    vwnd    = addfile("/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc","r")
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
    timetp = tp->time
    timeu  = uwnd->time
    timev  = vwnd->time
    time_tp=cd_calendar(timetp, -2)
    time_u=cd_calendar(timeu, -2)
    time_v=cd_calendar(timev, -2)
	; 2018年6月27日
	timetp1    =ind(time_tp.eq.20180627)
    timeu1     =ind(time_u.eq.20180627)
	timev1     =ind(time_v.eq.20180627)
	; 2018年6月15日至7月15日
	timetp2    =ind(time_tp.ge.20180615.and.time_tp.le.20180715)
	timeu2     =ind(time_u.ge.20180615.and.time_u.le.20180715)
	timev2     =ind(time_v.ge.20180615.and.time_v.le.20180715)
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
    rain1m  = tp->tp
    rain1mm =rain1m(timetp1,:,:)
    rain1mm =1000*rain1mm
	
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
    rain2m  = tp->tp
    rain2mm=rain2m(timetp2,{35},{127.5})
    rain2mm =1000*rain2mm
	
	
	; 提取 2018年6月27日 700百帕的风场
	u1=uwnd->uwnd(timeu1,{700},:,:)
    v1=vwnd->vwnd(timev1,{700},:,:)
	
	
	
	; ========================================================================
	; 打开工作空间
	wkstype="eps"
    wksname="18304032.assignment3.q1.test"
	wks=gsn_open_wks(wkstype, wksname)
	
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
	res@tiMainString    = "(a) 06/27 Precipitation & 750hPa Wind"

	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF              = 15.
    res@mpMaxLatF              = 45.
    res@mpMinLonF              = 100.
    res@mpMaxLonF              = 140.
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion      ="HighRes"
    ; res@mpDataResolution                = "FinestResolution"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette          ="CBR_wet"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性

    res@cnLevelSelectionMode    = "ManualLevels"
    res@cnMaxLevelValF          = 100.
    res@cnMinLevelValF          = -5.
    res@cnLevelSpacingF         = 15.
    res@cnLineThicknessF        = 1.5
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	

    res@tiMainFontHeightF           = 0.02
    res@tmXBLabelFontHeightF        = 0.015
    res@tmYLLabelFontHeightF        = 0.015
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	res@pmLabelBarParallelPosF=0.35
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnFillOn   =False
    rescn@cnLinesOn  =True
    rescn@cnLevelSelectionMode    ="ExplicitLevels"
    rescn@cnLevels=(/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnMonoLineColor      = False 
    rescn@cnLineColors          ="Black"
	rescn@cnLineThicknessF     =3.
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnInfoLabelOn         = False
    rescn@cnLineLabelsOn        = False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle          = "LineArrow"
    

   
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF = 1
    resvc@vcMinDistanceF        = 0.025
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF          = 0.1
    resvc@vcRefMagnitudeF       = 20.
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF=1.1
    resvc@vcRefAnnoOrthogonalPosF=-0.225
	
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
	
	resxy@tiMainString="(b) Marked point(35~S~o~N~N 127.5~S~o~N~E) Jun 15 - Jul 15  Precipitation"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString              ="time"
    resxy@tiYAxisString              ="Precipitation(mm)"
	
	
	
	; gsnXYBarChart
	
	resxy@gsnXYBarChart        =True
	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode         ="Explicit"
    resxy@tmXBValues = ispan(0, 30, 1)
    resxy@tmXBLabels = (/ "","Jun 16","","","","","Jun 21","","","","","Jun 26",\
                        "","","","","Jul 1","","","","","Jul 6","","","","","Jul 11","","","",""/)
   
	
	
	; (可选) XY线及柱状的颜色
	resxy@gsnAboveYRefLineColor =  "blue" 
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
    resxy@tmXBLabelFontHeightF            = 0.015
    resxy@tmYLLabelFontHeightF            = 0.015
    resxy@tiMainFontHeightF               = 0.02
    resxy@gsnLeftStringFontHeightF        = 0.015
    resxy@gsnRightStringFontHeightF       = 0.015
	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	resMark                     = True
    resMark@gsMarkerIndex       = 5
    resMark@gsMarkerThicknessF  = 3.
    resMark@gsMarkerColor       = "Red"
    resMark@gsMarkerSizeF       = 0.024

    
	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot=gsn_csm_contour_map_ce(wks, rain1mm, res)
	
	; 降水=30 mm
	plotline=gsn_csm_contour(wks, rain1mm, rescn)
    overlay(plot, plotline)
	
	; 700百帕风场
	plotwind=gsn_csm_vector(wks, u1, v1, resvc)
    overlay(plot, plotwind)
	
	; 添加标记
	plotmarker=gsn_add_polymarker(wks, plot, 127.5, 35, resMark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	
	 plotxy=gsn_csm_xy(wks,ispan(0, 30, 1),rain2mm, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	resPanel           =True
    plotpanel   = new(2,graphic)
    plotpanel(0)=plot
    plotpanel(1)=plotxy
	resPanel@gsnPanelMainString         = "18304032hzy"
    gsn_panel(wks,plotpanel,(/2,1/), resPanel)
    draw(plotpanel)
    frame(wks)
	; =========================================================================
	; 去掉空白空间
	delete (wks)            ; 必须先关闭 workstation
    system("convert -geometry 3000x3000 -density 600x600 -trim " + wksname + ".eps " + wksname + ".png")   
    system("rm " + wksname + ".eps -f") 

end