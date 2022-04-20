begin

	pathData 	= 	"/home/ys17-19/lsh/Project/other/shapefile/lecture9/"
	f 		= addfile(pathData + "slp.20050310-20050311.nc","r")
	slp 	= short2flt(f->slp)
	slp 	= (/slp/100./)
	
	dslp 	= slp(1,:,:)
	dslp 	= (/slp(1,:,:)-slp(0,:,:)/)
	
	; ================================================
	; add shapefiles
	
	pathShape 	= 	"/home/ys17-19/lsh/Project/other/shapefile/shapefile_China/"
	
	china 		= pathShape + "china.shp"
	scs 		= pathShape + "south_china_sea.shp"
	river 		= pathShape + "river.shp"
	province 	= pathShape + "province.shp"
	tp = "/home/ys17-19/lsh/Project/other/shapefile/shapefile_Tibet/DBATP_Polygon.shp"
	
	; ==========================================================================
	; open workstation
	
	wksType 	= "pdf"
	wksName 	= "slp_china"
	
	wks 		= gsn_open_wks(wksType,wksName)
	
	; ==========================================================================
	; set resources for drawing SLP
	
	res 						= True
	res@gsnDraw 				= False
	res@gsnFrame 				= False
	res@gsnMaximize 			= True
	
	res@cnFillOn 				= True

	res@cnLevelSelectionMode 	= "ManualLevels"
	res@cnMinLevelValF 			= -20
	res@cnMaxLevelValF 			= 20
	res@cnLevelSpacingF 		= 5
	res@cnFillPalette 			= "BlueWhiteOrangeRed"
	
	res@cnLinesOn 				= False
	
	res@mpMinLatF 				= 0.
	res@mpMaxLatF 				= 60.
	res@mpMinLonF 				= 70.
	res@mpMaxLonF 				= 140.	
	
	res@tmXBLabelFontHeightF 	= 0.02
	res@tmYLLabelFontHeightF 	= 0.02
	
	res@pmTickMarkDisplayMode 	= "Always"

	plot = gsn_csm_contour_map(wks, dslp, res)
	
	; set resources for shapefiles

	; =============================================================================
	; plot SLP
	


    resShape 					= True
	; add boundaries of China, Chinese provinces, and South China Sea

	resShape@gsLineColor= "black"
	
	plot_add = gsn_add_shapefile_polylines(wks, plot,china, resShape)
    resShape@gsLineColor= "black"
    plot_add2 = gsn_add_shapefile_polylines(wks, plot,scs, resShape)
    resShape@gsLineColor= "blue"
    plot_ad3 = gsn_add_shapefile_polylines(wks, plot,river, resShape)
    resShape@gsLineColor= "black"
    plot_add4 = gsn_add_shapefile_polylines(wks, plot,province, resShape)
     resShape@gsLineColor= "black"
     resShape@gsFillColor  = "gray"
     plot_add5 = gsn_add_shapefile_polygons(wks, plot,tp, resShape)



	draw(plot)
	frame(wks)
	
	
end