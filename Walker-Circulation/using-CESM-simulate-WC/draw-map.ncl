  begin 
    res = True
    res@gsnFrame = False
    res@tmXBOn  = False
    res@tmYLOn  = False
    res@tmYROn  = False
    res@tmXTOn  = False
   
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 30
    res@mpMaxLonF             = 390

    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
    res@mpFillOn = False
    res@gsnMaximize = True
  
  ; wks  = gsn_open_wks("png","maponly")    ; send graphics to PNG file
       ; draw global map

    pltType = "png"
    pltName = "map"
    wks = gsn_open_wks(pltType, pltName)

    plot = gsn_csm_map(wks,res)   
    frame(wks)

    system("convert -geometry 1800x1800 -density 600x600 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
  

    end 