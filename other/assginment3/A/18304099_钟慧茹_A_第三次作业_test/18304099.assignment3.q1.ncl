begin

	; 读取 .nc 文件

	tpfile="/home/ys17-19/lsh/Project/other/assginment3/tp.assignment3.nc"
    ufile ="/home/ys17-19/lsh/Project/other/assginment3/uwnd.assignment3.nc"
    vfile ="/home/ys17-19/lsh/Project/other/assginment3/vwnd.assignment3.nc"
    
    ftp   = addfile(tpfile,"r")
    fu    = addfile(ufile,"r")
    fv    = addfile(vfile,"r")
   
	tptime=ftp->time ;提取时间维
    utime =fu->time
    vtime =fv->time

    tptime_=cd_calendar(tptime, -2) 
    utime_ =cd_calendar(utime, -2)
    vtime_ =cd_calendar(vtime, -2)
    ; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	
	; 2018年6月27日
	tptime1   = ind(tptime_.eq.20180627)
    utime1    = ind(utime_.eq.20180627)
    vtime1    = ind(vtime_.eq.20180627)
	
	
	; 2018年6月15日至7月15日
	
	tptime2   =ind(tptime_.ge.20180615.and.tptime_.le.20180715)
	tp2time   = cd_calendar(ftp->time(tptime2),-2)
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	
	tp1= ftp->tp(tptime1,:,:)   
    tp1    = (/tp1*1000/)
    tp1@units             = "mm"
    tp1&longitude@units   = "degrees_east"
    tp1&latitude@units    = "degrees_north"
	lon = ftp->longitude
    
    ; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	
    tp2    = ftp->tp(tptime2,{35},{127.5})
    tp2    = (/tp2*1000/)
    tp2@units            = "mm"
	
	
	
	; 提取 2018年6月27日 700百帕的风场
	
	u               = fu->uwnd(time|utime1,{level|700},lat|:,lon|:)
    v               = fv->vwnd(time|vtime1,{level|700},lat|:,lon|:)
    u&lon@units     = "degrees_east"
    u&lat@units     = "degrees_north"
    v&lon@units     = "degrees_east"
    v&lat@units     = "degrees_north"

	
	
	; ========================================================================
	; 打开工作空间
	
    wksName ="18304099.assignment3.q1"
    wks     = gsn_open_wks("eps", wksName)
	
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
	res@gsnCenterString                 = "precipitation in Jun 27 and wind in 700hPa"
    res@gsnCenterStringFontHeightF      = 0.028
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	
	res@mpMinLatF          = 15.
    res@mpMaxLatF          = 45.
    res@mpMinLonF          = 100.
    res@mpMaxLonF          = 140.
    res@mpCenterLonF       = 120.
	
	; 高分辨率地图
	; mpDataBaseVersion
	
	res@mpDataBaseVersion   = "HighRes"
    
	
	; 色板
	; cnFillPalette
	
	res@cnFillPalette       = "precip3_16lev"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode    = "ExplicitLevels"
    res@cnLevels                = ispan(0, 200, 10)
    res@lbBoxEndCapStyle        = "TriangleBothEnds"
    res@lbOrientation           = "Vertical"
	
	
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
	rescn@cnLevelSelectionMode   = "ExplicitLevels"
    rescn@cnLevels               = 30.
    

	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor            = "darkorchid2"
    rescn@cnLineDashPattern      = 0
    rescn@cnLineThicknessF       = 2.5
    
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn         = False
    rescn@cnInfoLabelOn          = False
    
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	
	resvc@vcGlyphStyle          = "CurlyVector"
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
    resvc@vcLineArrowThicknessF  =  2.
	
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	
	resvc@vcRefLengthF           = 0.08
    resvc@vcRefMagnitudeF        = 25.
    
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF         = 1.2
    resvc@vcRefAnnoOrthogonalPosF       = 0.1
    
    resvc@vcRefAnnoPerimOn              = False
    resvc@vcRefAnnoBackgroundColor      = -1
    resvc@vcRefAnnoString2              = "m/s"
   
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF          =0.016
	
	
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
	
	resxy@gsnCenterString               = "total precipitation in Jun 15-Jul 15"
    resxy@gsnCenterStringFontHeightF    = 0.028
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	
	
	
	resxy@tiXAxisString                 = "time(day)"
    resxy@tiYAxisString                 = "precipitation (mm)"
    
	; gsnXYBarChart
	
	resxy@gsnXYBarChart                 = True
	
	; X 轴
	; tmXBMode 及相应属性

	resxy@tmXBOn                        = True
    resxy@tmXBMode                      = "Explicit"
    resxy@tmXBValues                    = (/1,6,11,16,21,26/)
	resxy@tmXBLabels                    = (/"6.16","6.21","6.26","7.1","7.6","7.11"/)
    resxy@tmXBMinorValues               = fspan(0,30,31)

	; (可选) XY线及柱状的颜色
	resxy@gsnXYBarChartColors  =  "cornflowerblue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	resxy@tiXAxisFontHeightF   =0.02
	resxy@tiYAxisFontHeightF   =0.02
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
	
	
	resmk                         = True
    resmk@gsMarkerIndex           = 7
    resmk@gsMarkerColor           = "red"
    resmk@gsMarkerSizeF           = 0.02
    resmk@gsMarkerThicknessF      = 3.0
    x               = 127.5 
    y               = 35
    
	; =========================================================================
	; 画图
	
	; 图 (a)
	plot               = new(2,graphic)
	; 降水地图
	
    plot(0)       = gsn_csm_contour_map(wks, tp1, res)
    
    ; 降水=30 mm
    plotline     = gsn_csm_contour(wks, tp1, rescn)
    
    ; 700百帕风场
    plotvector     = gsn_csm_vector(wks, u, v, resvc)
    
    ;  添加标记
    
    plotmark     = gsn_add_polymarker(wks, plot(0), x, y, resmk)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)

    x2 = fspan(0,30,31)
	plot(1)        = gsn_csm_xy(wks, x2, tp2, resxy)
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	
	
	overlay(plot(0), plotvector)
    overlay(plot(0), plotline)
   

    resPanel                           = True
    resPanel@gsnPanelMainFontHeightF   = 0.02
    ;resPanel@gsnPanelYF                = 0.7
    ;resPanel@gsnPanelTop               = 0.9
    ;resPanel@gsnPanelDebug             = True  
    resPanel@gsnPanelMainString        = "18304099.zhonghuiru" 

    gsn_panel(wks, plot,(/2,1/), resPanel)  ;输出展板
	; =========================================================================
	; 去掉空白空间

    delete(wks)
    system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
    system("/bin/rm -f " + wksName + ".eps")
  

end