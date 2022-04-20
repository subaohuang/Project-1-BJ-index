begin

	f1 		= addfile("tp.assignment3.nc","r")
	f2 		= addfile("uwnd.assignment3.nc","r")
	f3 		= addfile("vwnd.assignment3.nc","r")
	
	; ========================================================================
	
	time1 	= f1->time

	time1Arr = cd_calendar(time1,-2)
	time1Arrb = cd_calendar(time1,-5)  ;;;YYYY,MM,DD,HR,MIN,SEC
	dayArrb  = time1Arrb(:,2)          ;;;DD
	
	time1Ind = ind(time1Arr.eq.20180627)	 ;;;;index 2018.06.27的降水数据 

	time1Indb = ind(time1Arr.ge.20180615 .and. time1Arr.le.20180715) ;;;index 2018.06.15-27的降水数据 
	dayPlot   = dayArrb(time1Indb)          ;;;;day : 15-30,1-15
	
	
	dayIndAll = ind(dayPlot%1.eq.0)    ;;提取index


	dayInd 	  = ind((dayPlot-1)%5.eq.0)  ;;;提取5天为间隔的index

	
	time2 	= f2->time
	time2Arr = cd_calendar(time2,-2)
	time2Ind = ind(time2Arr.eq.20180627)
	
	time3 	= f3->time
	time3Arr = cd_calendar(time3,-2)
	time3Ind = ind(time3Arr.eq.20180627)
	
	; ========================================================================
	
	tp	 	= f1->tp(time1Ind,:,:)    ;;;2018.6.27的降水数据
	tp 		= (/tp*1000./)            ;;; mm convert to m  或者使用 copy_VarCoords
	
	u 		= f2->uwnd(time2Ind,{700},:,:)        ;;2018.6.27的u风数据
	v 		= f3->vwnd(time3Ind,{700},:,:)        ;;2018.6.27的v风数据
	
	tpTime 	= f1->tp(time1Indb,{35},{127.5})
	tpTime	= (/tpTime*1000./)
	
	cnLevels = (/1,5,10,20,30,40,50,60,80,100/)
	
	wksName = "precipitation_20180627"
	wks 	= gsn_open_wks("eps",wksName)
	
	; =========================================================================
	
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
	
	res@mpMinLatF 						= 15.
	res@mpMaxLatF 						= 45.
	res@mpMinLonF 						= 100.
	res@mpMaxLonF 						= 140.

	res@tmXBLabelFontHeightF 			= 0.028
	res@tmYLLabelFontHeightF 			= 0.028
	
	res@lbBoxEndCapStyle 				= "TriangleBothEnds"
	res@lbLabelAutoStride				= False
	res@pmTickMarkDisplayMode 			= "Always"
	
	res@pmLabelBarParallelPosF 			= 0.4
	res@pmLabelBarOrthogonalPosF 		= 0.125
	
	; res@mpDataBaseVersion				= "RANGS_GSHHS"
    res@mpDataBaseVersion               = "HighRes"
	res@mpDataResolution				= "FinestResolution"
	
	res@cnFillPalette 					= "precip_11lev"
	res@cnLevelSelectionMode 			= "ExplicitLevels"
	res@cnLevels 						= cnLevels
	
	res@mpGeophysicalLineColor			= "Gray40"
	res@mpGeophysicalLineThicknessF		= 0.8
	
	res@tiMainFontHeightF 				= 0.028
	res@gsnCenterStringFontHeightF		= 0.028
	res@tiMainOffsetYF					= -0.015
	
	; ---------------------------------------------------------------------------
	
	rescn@cnLevelSelectionMode			= "ExplicitLevels"
	rescn@cnLevels 						= 30.
	rescn@cnLineLabelsOn 				= False
	rescn@cnInfoLabelOn 				= False
	rescn@cnLineThicknessF 				= 2.
	rescn@cnLineColor 					= "Red"
	
	; ---------------------------------------------------------------------------
	
	resvc@vcGlyphStyle 					= "CurlyVector"
	resvc@vcRefMagnitudeF 				= 8.
	resvc@vcRefLengthF	  				= 0.02
	resvc@vcMinDistanceF				= 0.005
	resvc@vcLineArrowThicknessF 		= 2.
	
	resvc@vcRefAnnoString2On			= False
	resvc@vcRefAnnoParallelPosF			= 0.95
	resvc@vcRefAnnoOrthogonalPosF		= -0.25
	resvc@vcRefAnnoFontHeightF			= 0.02
	
	; ---------------------------------------------------------------------------
	
	resxy@gsnXYBarChart					= True
	resxy@vpWidthF 						= 0.7
	resxy@vpHeightF 					= 0.4
	
	resxy@tmXTOn 						= False
	resxy@tmYROn 						= False
	
	resxy@tmXBLabelFontHeightF 			= 0.028
	resxy@tmYLLabelFontHeightF 			= 0.028
	
	resxy@trYMinF 						= 0.
	
	resxy@tmXBMode 						= "Explicit"
	resxy@tmXBValues 					= dayInd
	resxy@tmXBLabels 					= dayPlot(dayInd)
	resxy@tmXBMinorValues 				= dayIndAll
	
	resxy@xyLineColor 					= "Black"
	resxy@gsnXYBarChartColors 			= "Blue"
	
	resxy@tiMainFontHeightF 			= 0.028
	resxy@tiXAxisFontHeightF 			= 0.028
	resxy@tiYAxisFontHeightF 			= 0.028
	
	; ---------------------------------------------------------------------------
	
	resMark								= True
	resMark@gsMarkerIndex 				= 5
	resMark@gsMarkerColor				= "White"
	resMark@gsMarkerSizeF 				= 0.008
	resMark@gsMarkerThicknessF 			= 3.
	
	; ---------------------------------------------------------------------------
	
	resShp								= True
	resShp@gsLineThicknessF				= 0.8
	resShp@gsLineColor					= "Gray40"
	
	; =========================================================================
	
	res@tiMainString 					= "~F21~(a) 27 Jun 2018"
	res@gsnCenterString 				= "precipitation (shading, mm/day) and 700-hPa wind (vector, m/s)"
	
	plot 								= gsn_csm_contour_map_ce(wks,tp,res)
	plotcn	 							= gsn_csm_contour(wks,tp,rescn)
	plotvc 								= gsn_csm_vector(wks,u,v,resvc)
	
	overlay(plot,plotcn)
	overlay(plot,plotvc)
	
	
	resxy@tiMainString 					= "(b) precipitation at 35~S~o~N~N,127.5~S~o~N~E (15 Jun-15 Jul 2018)"
	resxy@tiYAxisString 				= "precipitation (mm/day)"
	resxy@tiXAxisString 				= "date"
	
	plot2 								= gsn_csm_y(wks,tpTime,resxy)
	
	; =========================================================================
	
	plotMark							= gsn_add_polymarker(wks,plot,127.5,35,resMark)
	
	; =========================================================================
	
	resPanel	 						= True
	resPanel@gsnPanelMainString			= "Student Name (ID)"
	resPanel@gsnPanelMainFontHeightF 	= 0.014
	resPanel@gsnPanelYF					= (/0.86,0.43/)
	
	gsn_panel(wks,(/plot,plot2/),(/2,1/),resPanel)
	
	delete(wks)
	
	system("convert -geometry 1600x1600 -density 600x600 -trim " + wksName + ".eps " + wksName + ".png")
	system("/bin/rm -f " + wksName + ".eps")
	

end