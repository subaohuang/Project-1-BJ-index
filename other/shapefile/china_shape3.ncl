load "shapefile_utils.ncl"

begin
	
	f 		= addfile("MIX_Asia_mask.generic.025x025.nc","r")
	fShape  = addfile("/home/ys17-19/lsh/Project/other/shapefile/shapefile_China/china.shp","r")

	time 			= f->time
	lat 			= f->lat
	lon 			= f->lon
	MASK3D 			= f->MASK(:,:,:)
	MASK1 			= f->MASK(0,:,:)

	opt 			= True
	opt@shape_var	= "NAME"				
	opt@shape_names	= (/"China"/)
	opt@return_mask = True
	
	fShapeName 		= "shapefile_China/china.shp"
	China	 		= shapefile_mask_data(MASK1,fShapeName,opt)
	MASK2D			= mask(MASK1,China,1)

	; MASK = where(MASK.eq.MASK@_FillValue, 0, MASK)
	MASK2D 			= where(ismissing(MASK2D), 0, MASK2D)
	copy_VarCoords(MASK1,MASK2D)

	MASK 			= conform_dims(dimsizes(MASK3D), MASK2D, (/1,2/))
	copy_VarCoords(MASK3D, MASK)
	printVarSummary(MASK)

	system("date")

    ; 画图验证
 	; res                  = True
	; res@mpMinLatF 		 = 10.
	; res@mpMaxLatF 		 = 60.
	; res@mpMinLonF 	     = 60.
	; res@mpMaxLonF 	     = 150.
	
 	; res@cnFillOn         = True
 	; res@cnFillPalette    = "BlueWhiteOrangeRed"
	
	; res@mpDataBaseVersion = "RANGS_GSHHS"
	; res@mpDataResolution	= "FinestResolution"

	; wks                  = gsn_open_wks("png","map2")

 	; plot                 = gsn_csm_contour_map_ce(wks,MASK(0,:,:),res)

 	;写入新的mask文件
 	; file_name = "MIX_Asia_mask.generic.025x025_China.nc"
 	; system("rm -f" + file_name)
 	; file = addfile(file_name, "c")

 	; ntim = 1
 	; nlat = 720
 	; nlon = 1440

 	; dim_names = (/"time","lat","lon"/)
 	; dim_sizes = (/ntim, nlat, nlon/)
 	; dim_unlim = (/False,False,False/)

 	; filedimdef(file, dim_names, dim_sizes, dim_unlim)

 	; filevardef(file, "time", typeof(time), getvardims(time))
 	; filevardef(file, "lat", typeof(lat), getvardims(lat))
 	; filevardef(file, "lon", typeof(lon), getvardims(lon))
 	; filevardef(file, "MASK", typeof(MASK), getvardims(MASK))

 	; filevarattdef(file, "time", time)
 	; filevarattdef(file, "lat", lat)
 	; filevarattdef(file, "lon", lon)
 	; filevarattdef(file, "MASK", MASK)

 	; file->time = (/time/)
 	; file->lat = (/lat/)
 	; file->lon = (/lon/)
 	; file->MASK = (/MASK/)

 	filename = "./MIX_Asia_mask.generic.025x025.nc"
 	f2 = addfile(filename, "w")
 	

 	f2->MASK = MASK
print("done")
end