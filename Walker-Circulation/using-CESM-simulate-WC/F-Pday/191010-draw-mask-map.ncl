;test
begin
  ;;;----------------------------------------------------------读入数据
    fils1 = "~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc"
    f1 = addfile(fils1, "r")
    time = f1->time
    lon = f1->lon
    lat = f1->lat
    
    ; SST_cpl = f1->SST_cpl(:,:,:)
    SST_cpl = lonFlip(f1->SST_cpl(:,:,:))
    
  ;;;------去掉大西洋的海温增温趋势，海盆的大小定义为30N-30S ,80w(280E)-20E---------

     ;; 将0-360E的数据先换为-180-180的数据  
     ; SST_cpl_atlantic = SST_cpl(:,{-30:30},{-80:20})
     SST_cpl_atlantic = SST_cpl(:,:,:)
     

      printVarSummary(SST_cpl_atlantic)

     ;;-------cesm中的海温在陆地上也有数据，所以这里不能将陆地的温度也去掉--------
       fils2 = "~/lsh/data/landsea.nc"
       f2 = addfile(fils2, "r")
       ; lsdata = lonFlip(f2->LSMASK)

       lsdata = f2->LSMASK
       printVarSummary(lsdata)
       
       lsm = landsea_mask(lsdata,SST_cpl_atlantic&lat,SST_cpl_atlantic&lon) 
   
       printVarSummary(lsm)
      

      wks = gsn_open_wks("pdf", "mask_map")
      
  cmap = read_colormap_file("wh-bl-gr-ye-re") ; read color data
  
    res = True
    res@mpFillOn         = False                ; do not color-fill the map
    res@gsnTickMarksOn   = False                ; turn off all tick marks
    res@mpPerimOn        = True                 ; turn the map perimeter on
    res@mpPerimDrawOrder = "PostDraw"           ; draw the map perimeter last
    res@gsnDraw          = False                ; do not draw the plot
    res@gsnFrame         = False                ; do not advance the frame
    res@cnLinesOn        = False                ; turn off the contour lines
    res@cnLineLabelsOn   = False                ; turn off the contour line labels
    res@cnLevelSelectionMode = "ExplicitLevels" ; explicitly set the levels via cnLevels
    res@cnLevels         = (/1.,2.,3.,4./)      ; set the levels
    res@cnFillOn         = True                 ; turn on color fill
    res@cnFillPalette    = cmap( (/60,100,20,140,5/), :)     ; use subsetted color map
    ; res@cnFillMode       = "RasterFill"         ; use raster fill         
    ;  res@cnFillColors     = (/60,100,20,140,5/)  ; set the colors that will be used to color fill
    res@lbLabelStrings   = ispan(0,4,1)         ; labels for the labelbar boxes
    res@lbLabelAlignment = "BoxCenters"         ; put the labels in the center of the label bar boxes
    res@lbTitleString    = "0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf"  ; labelbar title
    res@lbTitleFontHeightF = 0.0125             ; labelbar title font height


  ; res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = min(lsm&lat)
  ; res@mpMaxLatF             = max(lsm&lat)
  ; res@mpMinLonF             = min(lsm&lon)
  ; res@mpMaxLonF             = max(lsm&lon)
  ; res@mpCenterLonF          =avg(lsm&lon)



       plot = gsn_contour(wks, lsm, res)

       draw(plot)
       frame(wks)



end 

   