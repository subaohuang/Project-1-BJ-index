begin

	; 读取 .nc 文件
    dir="/home/ys17-19/lsh/Project/other/assginment3/"
    fu=addfile(dir+"uwnd.assignment3.nc", "r")
    fv=addfile(dir+"vwnd.assignment3.nc", "r")
	ftp=addfile(dir+"tp.assignment3.nc", "r")
	 ;print(fu)
     ;print(fv)
 ;    print(ftp)
	; ========================================================================
	; 提取 .nc 文件指定的时间位置
	; (建议使用 cd_calendar 及 ind)
	t=fu->time
    time=cd_calendar(t, -2)

    ttp=ftp->time
    time_tp=cd_calendar(ttp, -3)
    ;print(time_tp)
	; 2018年6月27日
    time1_tpind=ind(time_tp.eq.2018062711)
    time1_uvind=ind(time.eq.20180627)   ;
	;print(time1_ind)
	
	; 2018年6月15日至7月15日  
    time2_tpind=ind(time_tp.le.2018071711.and.time_tp.ge.2018061511)

    ;print(time2_ind)
	
	
	; ========================================================================
	; 提取 2018年6月27日的降水数据
	; 单位需要由 m 转为 mm
	time1_tp_m=ftp->tp(time1_tpind,{15:45},{100:140})

    time1_tp_mm=time1_tp_m
    time1_tp_mm=(/time1_tp_m(:,:)*1000/)
    time1_tp_mm@unit="mm"
    printVarSummary(time1_tp_mm)
    
	
	; 提取 2018年6月15日至7月15日北纬35度，东经127.5的降水数据
	; 单位需要由 m 转为 mm
	time2_tp_m=ftp->tp(time2_tpind,{35},{127.5})
	time2=ftp->time(time2_tpind)
    t2=cd_calendar(time2, -2)
	time2_tp_mm=time2_tp_m
    time2_tp_mm=(/time2_tp_m(:)*1000/)
    time2_tp_mm@unit="mm"
	
	; 提取 2018年6月27日 700百帕的风场
	u=fu->uwnd ( time1_uvind, {700}, {15:45},{100:140})
    v=fv->vwnd ( time1_uvind, {700}, {15:45},{100:140})
	
	
	
	; ========================================================================
	; 打开工作空间
    pltName="18304018.assignment3.q1"
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
				
	res@vpWidthF 						= 0.9
				
	res@cnFillOn 						= True		
	res@cnLinesOn						= False
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; res: 图 (a) 底图
	
	; 标题
	res@tiMainString                     =" precipitation and 700hPa wind of June 27th(2018)"
	res@tiMainFontThicknessF          =1.0
	
	; 地图范围
	; mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF
	res@mpMinLatF                        =15
    res@mpMaxLatF                        =45
    res@mpMinLonF                        =100
    res@mpMaxLonF                        =140
	
	
	; 高分辨率地图
	; mpDataBaseVersion
	res@mpDataBaseVersion                 ="HighRes"
	
	
	; 色板
	; cnFillPalette
	res@cnFillPalette                     ="cmp_b2r"
	
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	res@cnLevelSelectionMode              ="ManualLevels"
    printMinMax(time1_tp_mm, 0)
	res@cnMaxLevelValF                    =120
    res@cnMinLevelValF                    =0
	res@cnLevelSpacingF                   =5

	; (可选) 字体大小
	; XY轴：tmXBLabelFontHeightF, tmYLLabelFontHeightF
	; 主标题: tiMainFontHeightF
	; 副标题: gsnCenterStringFontHeightF
	res@tiMainFontHeightF                  =0.02
   
	; (可选) 主标题Y方向移动位置
	; tiMainOffsetYF
	
	
	; ---------------------------------------------------------------------------
	; rescn: 图 (a) 等值线：降水=30mm
	
	; 等值线间距
	; cnLevelSelectionMode 及相应属性
	rescn@cnLevelSelectionMode     ="ExplicitLevels"
    rescn@cnLevels                 =(/30/)
	
	; 等值线颜色及粗细度
	; cnLineColor, cnLineThicknessF
	rescn@cnLineColor              ="springgreen4"
    rescn@cnLineThicknessF         =3

	; 不显示等值线标签及信息
	; cnLineLabelsOn, cnInfoLabelOn
	rescn@cnLineLabelsOn             =False
    rescn@cnInfoLabelOn              =False
	
	; ---------------------------------------------------------------------------
	; resvc: 图 (a) 矢量
	
	; 矢量样式
	; vcGlyphStyle
	resvc@vcGlyphStyle                 ="CurlyVector"
	
	
	; (可选) 矢量的粗细度及最短距离(NDC单位)
	; vcLineArrowThicknessF, vcMinDistanceF
	resvc@vcLineArrowThicknessF        =1.
    resvc@vcMinDistanceF               =0.01
	
	

	; 参考矢量的长度(NDC单位)及其表示的数值
	; vcRefLengthF, vcRefMagnitudeF
	resvc@vcRefLengthF     =0.03
    resvc@vcRefMagnitudeF  =5
	
	; 参考矢量信息的X,Y方向的位置
	; vcRefAnnoParallelPosF, vcRefAnnoOrthogonalPosF
	;resvc@vcRefAnnoParallelPosF         =0.5
    resvc@vcRefAnnoOrthogonalPosF        =-1.6
	
	
	; (可选) 参考矢量信息的字体大小
	;vcRefAnnoFontHeightF
    resvc@vcRefAnnoString2On            =False
	resvc@vcRefAnnoString1On             =True
    resvc@vcRefAnnoString1              ="5 m/s"
	
	; ---------------------------------------------------------------------------
	; resxy: 图 (b) XY 图
	
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@trYMinF 						= 0
	
	; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	; 以柱状图表示
	; 标题
	resxy@tiMainString                  ="precipitation from June 15th to July 15th (2018)"
	
	; X,Y轴标题
	; tiXAxisString, tiYAxisString
	resxy@tiXAxisString="Day"
    resxy@tiYAxisString="precipitation (mm)"
	
	
	
	; gsnXYBarChart
	resxy@gsnXYBarChart                 =True

	
	; X 轴
	; tmXBMode 及相应属性
    resxy@tmXBMode                      ="Explicit"
    resxy@tmXBValues                    =ispan(5,30,5)
    resxy@tmXBLabels                    =(/"Jun 16","Jun 21" ,"Jun 26","Jul 1","Jul 6","Jul 11"/)


    resxy@tmXTOn                        =True
    ;res@tmXTBorderOn                    =True
    resxy@tmXTLabelsOn                  =True
    resxy@tmXUseBottom                  =False
    resxy@tmXTMode                      ="Explicit"
    resxy@tmXTValues                    =ispan(1,29,1)
    resxy@tmXTLabels                    =(/"Jun16","17","18","19","20","21","22","23","24","25","26","27","28","29","30",\
                                          "Jul1","2","3","4","5","6","7","8","9","10","11","12","13","14"/)
	resxy@tmXTLabelFontHeightF          =0.01
    ; (可选) XY线及柱状的颜色
	resxy@gsnAboveYRefLineColor           ="lightskyblue"
	
	
	; (可选) 字体大小
	; 主标题: tiMainFontHeightF
	; X,Y轴标题: tiXAxisFontHeightF, tiYAxisFontHeightF
	
	resxy@tiMainFontHeightF                =0.02
	
	; ---------------------------------------------------------------------------
	; 标记 (北纬35度，东经127.5度)
	; gsMarkerIndex, gsMarkerColor, gsMarkerSizeF
	; (可选) gsMarkerThicknessF
    resMark                     = True
    resMark@gsMarkerIndex       = 12
    resMark@gsMarkerThicknessF  =4.
    resMark@gsMarkerColor       = "Blue"
    resMark@gsMarkerSizeF       = 0.024

	
	; =========================================================================
	; 画图
	
	; 图 (a)
	
	; 降水地图
    plot=new(2,graphic)
    plot(0)=gsn_csm_contour_map(wks, time1_tp_mm, res)
	
	;降水=30 mm
    overlay_plot1=gsn_csm_contour(wks,time1_tp_mm , rescn)
	overlay(plot(0),overlay_plot1)
	
	; 700百帕风场
    ;printVarSummary(u)
    ;printVarSummary(v)
	overlay_plot2=gsn_csm_vector(wks, u, v, resvc)
	overlay(plot(0),overlay_plot2)
	; 添加标记
	plotM      = gsn_add_polymarker(wks,plot(0),127.5,35.,resMark)

	
	; --------------------------------------------------------------------------
	; 图 (b)

	plot(1)=gsn_csm_y(wks,  time2_tp_mm, resxy)
	
	
	
	; =========================================================================
	; 绘画展板
	; gsnPanelXF 设定图(a)及(b)在X方向的位置(NDC坐标)
	; gsnPanelYF 设定图(a)及(b)在Y方向的位置(NDC坐标)
    respanel=True
    respanel@gsnPanelXF=0.2
    respanel@gsnPanelYF=0.2
    respanel@gsnPanelMainString="18304018 Chong Mei"
	gsn_panel(wks, plot, (/2,1/), respanel)

    delete (wks)           ; 必须先关闭 workstation
    system("convert -geometry 3000x3000 -density 600x600 -trim " + pltName + ".eps " + pltName + ".png")   
    system("rm " + pltName + ".eps -f") 
	; =========================================================================
	; 去掉空白空间


    ;运行前先设置环境变量：export NCARG_RANGS="/cygdrive/C/users/阿欣/ncl/rangs/"


end