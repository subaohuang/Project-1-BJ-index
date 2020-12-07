begin

    ; 读取 .nc 文件
    path="/home/ys17-19/lsh/Project/other/assginment3/"
    fil1="uwnd.assignment3.nc"
    fil2="vwnd.assignment3.nc"
    fil3="tp.assignment3.nc"
    u=addfile(path+fil1, "r")
    v=addfile(path+fil2, "r")
    t=addfile(path+fil3, "r")
    uwnd=u->uwnd
    vwnd=v->vwnd
    tp=t->tp
    time1=u->time
    time2=t->time

    
    
    ; ========================================================================
    ; 提取 .nc 文件指定的时间位置
    ; (建议使用 cd_calendar 及 ind)
    ctime1=cd_calendar(time1, 2)
    ctime2=cd_calendar(time2, 2)

    ; 2018年6月27日
    time1_index=ind(ctime1.eq.20180627)
    time2_index=ind(ctime2.eq.20180627)

    
    
    
    ; 2018年6月15日至7月15日
    time3_index=ind(ctime2.ge.20180615.and.ctime2.le.20180715)
    time_final=ctime2(time3_index)-2018

	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1=1000*tp(time2_index,:,:)
	copy_VarMeta(tp(time2_index,:,:), tp1)

	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm

	tp2=1000*tp(time3_index,{35},{127.5})
	
	; 提取 2018年6月27日 700百帕的风场
	
  	us=uwnd(time1_index,{700},:,:)

  	vs=vwnd(time1_index,{700},:,:)


	; ========================================================================
	; 打开工作空间
	pltName="18304098.assignment3.q1.ncl "
    wks = gsn_open_wks("png", pltName)
	
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	
	rescn                               =res
	resvc                               =res
	resxy                               =res
	

	
	
	res@vpWidthF 						= 0.7
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
  
	
	; 标题
    res@tiMainString="20180627 precipitation and 700hPa wind"
    res@gsnLeftString="precipitation (mm)"
    res@gsnRightString="speed (m/s)"
	
	
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF=15.
	res@mpMaxLatF=45.
	res@mpMinLonF=100.
	res@mpMaxLonF=140.
    
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"

	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette="MPL_BrBG"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode="ManualLevels"
	res@cnMaxLevelValF=205.
	res@cnMinLevelValF=0.
	res@cnLevelSpacingF=10.
	
	
	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF=0.02
	res@tmYLLabelFontHeightF=0.02
	res@tiMainFontHeightF=0.034
	res@gsnCenterStringFontHeightF=0.02

	
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode="ExplicitLevels"
	rescn@cnLevels=30.


	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor="Red"
	rescn@cnLineThicknessF=3.5
	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn=False
	rescn@cnInfoLabelOn=False
	
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle="CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF=0.025
	resvc@vcMinDistanceF=0.025

	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF=0.1
	resvc@vcRefMagnitudeF=20
	
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	resvc@vcRefAnnoParallelPosF=0.1
	resvc@vcRefAnnoOrthogonalPosF=0.2
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
	resvc@vcRefAnnoFontHeightF=0.02
	
	
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
	
	
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="20180615-20180715"
	resxy@tiYAxisString="precipitation in 35~S~o~N~N,127.5~S~o~N~E"
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart=True
	
	
	; X 轴
	; tmXBMode 及相应属性
	resxy@tmXBMode="Explicit"
	resxy@tmXBValues=(/616,621,626,701,706,711/)
	resxy@tmXBMinorValues=time_final
	
	; (可选) XY线及柱状的颜色
	color=read_colormap_file("circular_1")
    ; resxy@xyLineColors=color(0:29)

	

	
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF=0.034
	resxy@tiXAxisFontHeightF=0.02
	resxy@tiYAxisFontHeightF=0.02

	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	;可选) gsMarkerThicknessF
	lat=35.
	lon=127.5
	resmark=True
	resmark@gsMarkerIndex=8
	resmark@gsMarkerColor="Blue"
	resmark@gsMarkerSizeF=0.024
	resmark@gsMarkerThicknessF=3.

	respanel=True
	respanel@gsnPanelMainString="18304098ZLX"
	

	
	
	; =========================================================================
	; 画图
	plot=new(2,graphic)
	; 图 (a)
	
	
	; 降水地图
	
	; plot(0)=gsn_contour_map(wks,tp1,res)
	; plot(0) = gsn_contour_map(wks,tp1,res)
	; ; 降水=30 mm
	; plot2=gsn_contour(wks,tp1,rescn)
	
	; ; 700百帕风场
	; plot3=gsn_vector(wks,us,vs,resvc)

    plot(0)=gsn_csm_contour_map_ce(wks,tp1,res)
	
	; 降水=30 mm
	plot2=gsn_csm_contour(wks,tp1,rescn)
	
	; 700百帕风场
	plot3=gsn_csm_vector(wks,us,vs,resvc)



	overlay(plot(0),plot2)
	overlay(plot(0),plot3)

	
	; 添加标记
	plotmark=gsn_add_polymarker(wks, plot(0),lat,lon,resmark)
	
	
	; --------------------------------------------------------------------------
	; 图 (b)
	plot(1)=gsn_csm_xy(wks,time_final,tp2,resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	gsn_panel(wks, plot, (/1,2/),respanel)
	
	
	
	; =========================================================================
	; 去掉空白空间
	; delete (wks) 			
	; system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
	; system("rm " + pltName + ".eps -f") 



end