begin
; 1.  利用 tp.assignment3.nc, uwnd.assignment3.nc和vwnd.assignment3.nc完成assignment3.template.ncl 
; 并且绘画展板显示图 (a) 及图 (b):


; b)  图 (b) 绘画2018年6月15日至7月15日北纬35度，东经127.5度降水的时间序列
; i.  使用柱状图
; ii. X轴主刻度显示 6月16, 21, 26日及7月1, 6及11日，副刻度显示所有日子 
; iii.    X，Y轴加上合适的标题
; c)  展板的标题为学号及拼音，图 (a) 及图 (b) 加上合适标题

	; 读取 .nc 文件 time*lecel*lat*lon
filename1="uwnd.assignment3.nc"
filename2="vwnd.assignment3.nc"
filename3="tp.assignment3.nc"

    f1=addfile(filename1, "r")
    f2=addfile(filename2, "r")
    f3=addfile(filename3, "r")
    t1=f1->time
    t2=f3->time

    tp=f3->tp
    uwnd=f1->uwnd
    vwnd=f2->vwnd
    ; printVarSummary(tp)
    ; printVarSummary(uwnd)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	t1:=cd_calendar(t1, -2)
    t2:=cd_calendar(t2, -2)

    ; print(d1)
    ; print(d2)

	; 2018年6月27日
	d1=ind(t1.eq.20180627) ;风场
    d2=ind(t2.eq.20180627) ;降水
	
	
	; 2018年6月15日至7月15日降水
	d3=ind(t2.ge.20180615.and.t2.le.20180715)
    ; print(d3)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	tp1=tp(d2,{15:45},{100:140})
    p=tp1*1000
    p!0="lat" 
    p&lat=tp1&latitude
    p!1="lon"
    p&lon=tp1&longitude
    ; printVarSummary(p)
    p@units="mm"
    p@long_name="tp(mm)"
    tp1:=p
    ; printVarSummary(tp1)
    tt=ispan(1, 31, 1)


	
 

	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	tp2=tp(d3,{35},{127.5})
    printVarSummary(tp2)
    p2=tp2*1000
    p2!0="time" 
    p2&time=tp2&time


    p2@units="mm"
    p2@long_name="tp(mm)"
    tp2:=p2
    printVarSummary(tp2)

 


	; 提取 2018年6月27日 700百帕的风场
	u=uwnd(d1,{700},{15:45},{100:140})
    v=vwnd(d1,{700},{15:45},{100:140})
	
	; ========================================================================
	; 打开工作空间
	wks=gsn_open_wks("eps", "18304080.assignment3.q1")
	
	
	; =========================================================================
	; 定义绘图属性
	
	res 								= True
	res@gsnDraw 						= False
	res@gsnFrame 						= False
	res@gsnLeftString 					= ""
	res@gsnRightString 					= ""
	res@gsnCenterString=" "
	resxy	 							= res
	
	res@gsnAddCyclic					= False
	rescn								= res
	resvc 								= res
				
	 res@vpWidthF 						= 0.7

				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	res@tmXTOn=False
    res@tmYROn=False
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@gsnCenterString="(a) Precipitation and wind field of 700hPa of east Asia in 20180627"
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMaxLatF=45
    res@mpMinLatF=15
    res@mpMaxLonF=140
    res@mpMinLonF=100
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion="HighRes"

	; ; 色板
	; cnFillPalette
	res@cnFillPalette="BlueWhiteOrangeRed"
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    res@cnLevelSelectionMode="ExplicitLevels"
	res@cnLevels = ispan(0,160,15)

	
	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@gsnCenterStringFontHeightF=0.025
	
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
    

	; 等值线间距
	; cnLevelSelectionMode 及相应属性
    rescn@cnLevelSelectionMode="ExplicitLevels"
    rescn@cnLevels          = 30



    rescn@cnLabelMasking=True
    rescn@cnLineLabelBackgroundColor=-1
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
    rescn@cnLineColor           = "Green"	
    rescn@cnLineThicknessF      = 3.5	
	
	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
    rescn@cnInfoLabelOn         = False
	rescn@cnLineLabelsOn=False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量

	; 矢量样式
	; vcGlyphStyle
    resvc@vcGlyphStyle="CurlyVector"	
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF=2
    resvc@vcMinDistanceF=0.025
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF

	resvc@vcRefLengthF = 0.1
    resvc@vcRefMagnitudeF = 20
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF 
    resvc@vcRefAnnoOrthogonalPosF = -0.45
    ; resvc@vcRefAnnoOn=True

	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoFontHeightF = 0.01
    resvc@vcRefAnnoString2On = False
	resvc@vcRefAnnoBackgroundColor=-1
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0.
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题2018年6月15日至7月15日北纬35度，东经127.5度
    resxy@gsnCenterString = "(b) Precipitation on Jun15 to Jul15, 2018 (35~S~o~N~N,127.5~S~o~N~E)"
    resxy@gsnCenterStringOrthogonalPosF = 0.04

	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="time(Day)"
	resxy@tiYAxisString="Precipitation(mm)"
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart=True
    resxy@gsnXYBarChartBarWidth=1
  
	
	
	; X 轴
    ; b)  图 (b) 绘画2018年6月15日至7月15日北纬35度，东经127.5度降水的时间序列
; i.  使用柱状图
; ii. X轴主刻度显示 6月16, 21, 26日及7月1, 6及11日，副刻度显示所有日子 
; iii.    X，Y轴加上合适的标题
; c)  展板的标题为学号及拼音，图 (a) 及图 (b) 加上合适标题
	; tmXBMode 及相应属性
    resxy@tmXBMode="Explicit"
	resxy@tmXBValues=(/2,7,12,17,22,27/)
	resxy@tmXBLabels=(/"6-16","6-21","6-26","7-01","7-06","7-11"/)
    resxy@tmXBMinorOn=True
    resxy@tmXBMinorValues=ispan(1, 31, 1)
	; (可选) XY线及柱状的颜色
    
	  resxy@gsnXYBarChartColors="Blue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	resxy@tiMainFontHeightF=0.02
    resxy@gsnCenterStringFontHeightF = 0.025
    resxy@tiXAxisFontHeightF = 0.025
    resxy@tiYAxisFontHeightF = 0.025
    resxy@tmXBLabelFontHeightF = 0.02
    resxy@tmYLLabelFontHeightF = 0.02

	
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF

    resMark             = True
; resMark@gsMarkerIndex       = 5  
resMark@gsMarkerThicknessF  = 3 
  
resMark@gsMarkerIndex       = 15
resMark@gsMarkerSizeF       =20
resMark@gsMarkerColor       = "Red"
x=127.5
y=35

	
	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
	plot1=gsn_csm_contour_map_ce(wks, tp1, res)
	; 降水=30 mm
	plotcn1=gsn_csm_contour(wks, tp1, rescn)
    overlay(plot1, plotcn1) 
	; ; 700百帕风场
	plotvc=gsn_csm_vector(wks, u, v, resvc)

    overlay(plot1, plotvc)

	; 添加标记
plotM1        = gsn_add_polymarker(wks,plot1,x,y,resMark)




; draw(plot1)
	; frame(wks)
	; --------------------------------------------------------------------------
	; 图 (b)
	plot2=gsn_csm_xy(wks, tt,tp2, resxy)
; 	draw(plot2)
; frame(wks)
	
	; c)  展板的标题为学号及拼音，图 (a) 及图 (b) 加上合适标题
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
	respanel=True
    ; respanel@gsnPanelXF=(/0.4,0.4/)
    ; respanel@gsnPanelYF=(/0.9,0.5/)
	; respanel@gsnPanelMainString="18304080+wulei"
 ;    respanel@gsnPanelMainFontHeightF=0.02
 ; respanel@gsnPanelXWhiteSpacePercent=5
	;    gsn_panel(wks,(/plot1,plot2/),(/2,1/),respanel)
            ; respanel = True
    respanel@gsnPanelMainString = "18304080(wu lei)"
    respanel@gsnPanelMainFontHeightF = 0.02
    respanel@gsnPanelYWhiteSpacePercent = 2.

  gsn_panel(wks, (/plot1,plot2/), (/2,1/), respanel)
    
    
        ; draw(plot2)
    ;     ; frame(wks)
    ; ; 转换图像类型
    delete(wks)

    system("convert -geometry 1600x1600 -density 600x600 -trim " + "18304080.assignment3.q1" + ".eps " + "18304080.assignment3.q1" + ".png")
    system("/bin/rm -f " + "18304080.assignment3.q1" + ".eps")

	
	

end