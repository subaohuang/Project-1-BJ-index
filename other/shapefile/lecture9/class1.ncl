begin
	
	; 读取文件
	
	f  		= addfile("wind10m.msl.20180916.06UTC.nc","r")

	; printVarSummary(f)
	; 提取数据
	
	u10 	= short2flt(f->u10(0,:,:))
	v10 	= short2flt(f->v10(0,:,:))
	printVarSummary(u10)
	; 运算风速
	
	wind10	= wind_speed(u10,v10)
	wind10  = (/wind10*3.6/)
    ; print(wind10)
	printMinMax(wind10, 0)
	; ============================================
	; 打开工作空间
	
	wksName = "16 Sep 2018 06UTC"
	wks 	= gsn_open_wks("X11",wksName)
    ; display
	; ============================================
	; 定义绘图属性

	res 				= True
	res@gsnDraw 		= False
	res@gsnFrame 		= False
	res@gsnLeftString 	= ""
	res@gsnRightString 	= ""
	
	rescn 				= True
	
	res@gsnAddCyclic				= False
	res@tmXTOn 						= False
	res@tmYROn 						= False
	
	res@mpGeophysicalLineThicknessF = 0.8
	res@pmTickMarkDisplayMode 		= "Always"
	
	res@lbLabelAutoStride			= False
	res@lbBoxEndCapStyle			= "TriangleBothEnds"
	
	; ！！！！！！！！！！！！！！！！！！！！！
	; 完成以下部分
	
	
	; 加上标题 (tiMainString, gsnCenterString)
	res@tiMainString                 ="10-m wind speed (km/h)"
	res@gsnCenterString              ="16 Sep 2018 06UTC"
	
	
	; 地图范围 (mpMinLatF, mpMaxLatF, mpMinLonF, mpMaxLonF)
	res@mpMinLatF                    =10
	res@mpMaxLatF                    =30
	res@mpMinLonF                    =100
	res@mpMaxLonF                    =135
	
	
	
	
	; 以填色表示，不显示等值线 (cnFillOn, cnLinesOn)
	res@cnFillOn                     =True
	res@cnLinesOn                    =False
	

	; 设定等值线的方式 (cnLevelSelectionMode)
	; ManualLevels -> cnMinLevelValF, cnMaxLevelValF, cnLevelSpacingF
	; ExplicitLevels -> cnLevels
	; res@cnLevelSelectionMode         ="ManualLevels"
	; res@cnMinLevelValF               =12
	; res@cnMaxLevelValF               =140
	; res@cnLevelSpacingF              =1


	; 颜色表 (cnFillPalette)
	;用Wind_17lev
	res@cnFillPalette                ="wind_17lev"

	; 字体大小
	
	res@gsnCenterStringFontHeightF 	= 0.02
	res@tiMainFontHeightF 	 		= 0.015
	res@tmXBLabelFontHeightF 		= 0.015
	res@tmYLLabelFontHeightF 		= 0.015
	
	
	; (optional) 色板在 Y 轴移动的位置 (pmLabelBarOrthogonalPosF)
	res@pmLabelBarOrthogonalPosF    =0.15
	
	; (optional) 主标题在 Y 轴移动的位置 (tiMainOffsetYF)
	res@tiMainOffsetYF              =-0.015
	
	; ================================================
	; 等值线 31 km/h 和 63 km/h 的绘图属性
	
	; 设定等值线的方式 (cnLevelSelectionMode)
	; ManualLevels -> cnMinLevelValF, cnMaxLevelValF, cnLevelSpacingF
	; ExplicitLevels -> cnLevels
	res@cnLevelSelectionMode        ="Explicit"
	res@cnLevels                    =(/31,63/)
	


	; 等值线的粗细度 (cnLineThicknessF 或 cnLineThicknesses)
	res@cnLineThicknesses           =(/2,2/)
	
	; 等值线的颜色 (cnLineColors)
	
	rescn@cnMonoLineColor		= False ; 必须用 cnMonoLineColor 才能用不同颜色
	res@cnLineColors            =(/"Black","Blue"/)
	
	; 不显示等值线的信息
	
	rescn@cnInfoLabelOn 		= False
	rescn@cnLineLabelsOn 		= False
	
	; ================================================
	; 画图
	plot                        =gsn_csm_contour_map(wks, wind10, res)
	draw(plot)
	frame(wks)
	end 