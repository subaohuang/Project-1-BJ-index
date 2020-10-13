begin
;;;Φ»╗σàÑhadisst
title=(/"industry","power","residential","transportation"/)
z_total=new((/4,120,130/),float)
plot=new(4,graphic)


wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/hang")
  do i=0,3
    print("*******************title***************="+title(i))
  diri="/mnt/f/share/hang/6月份BC排放源/2016_06__"+title(i)+"__BC.nc"
  f1=addfile(diri, "r")
  delete(diri)
  
  z=f1->z
  x_range=f1->x_range
  y_range=f1->y_range
  spacing=f1->spacing
  delete(f1)

  x=toint((x_range(1)-x_range(0))*10.0)
  y=toint((y_range(1)-y_range(0))*10.0)

  lat=fspan(y_range(0), y_range(1),y)
  lon=fspan(x_range(0), x_range(1),x)
 
  printVarSummary(z)
  z_mean_sum=reshape(z, (/y,x/)) 
  z_mean_sum=(/z_mean_sum(::-1,:)/)
  delete(x)
  delete(y)

  z_total(i,:,:)=z_mean_sum

  
  
  z_mean_sum!0="lat"
  z_mean_sum!1="lon"
  z_mean_sum&lat=lat
  z_mean_sum&lon=lon
  z_mean_sum&lat@units="degree_north"
  z_mean_sum&lon@units="degree_east"

     printVarSummary(z_mean_sum)

    cmap=read_colormap_file("MPL_hot")
    cmap=(/cmap(::-1,0:3)/)

    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic=False
    res@cnLinesOn=False   
    res@cnFillOn=True
    res@cnLevelSelectionMode="ExplicitLevels"
    ; res@cnLevels=(/0.01,0.1,1,2.5,5,7.5,10,15,20,30,40/) 
    res@cnLevels=(/0.01,0.1,1,2.5,5,7.5,10,15,20,30,40/) 
  
    ; res@cnFillPalette = cmap
    res@cnFillPalette = "BlAqGrYeOrReVi200"
    ; res@gsnSpreadColorStart=10
    ; res@gsnSpreadColorEnd=255
    res@mpFillDrawOrder="PreDraw"
    res@cnLineLabelsOn=False

    res@mpLimitMode="LatLon"
    res@mpMinLatF             = min(lat)
    res@mpMaxLatF             = max(lat)
    res@mpMinLonF             = min(lon)
    res@mpMaxLonF             = max(lon)

    res@mpDataSetName = "Earth..4"
    res@mpDataBaseVersion = "MediumRes"
    res@mpGeophysicalLineThicknessF= 4. ; double the thickness of geophysical boundaries 
    res@mpNationalLineThicknessF= 1.
    res@mpOutlineBoundarySets ="NoBoundaries"
    res@mpOutlineSpecifiers = (/"China:states","Thailand","Laos"/)
    res@lbLabelBarOn=False
    res@gsnCenterString=title(i)

    res@gsnCenterStringOrthogonalPosF=0.02
    if i.eq.3
      res@gsnCenterStringOrthogonalPosF=0.01
    end if

    ; res@pmLabelBarHeightF=0.05
    ; res@pmLabelBarOrthogonalPosF=0.15
    ; res@pmLabelBarParallelPosF=0.5
    ; res@lbTitleOn        =  True                ; turn on title
    ; res@lbTitleString    = "Mg/month/grid"              ; title string
    ; res@lbTitlePosition  = "Right"              ; title position
    ; res@lbTitleFontHeightF= .01              ; make title smaller
    ; res@lbTitleDirection = "Across"             ; title direction
    
    res@tmYROn=False
    res@tmXTOn=False
    res@tmXBOn=False
    res@tmYLOn=False


   plot(i)=gsn_csm_contour_map(wks, z_mean_sum, res)
   delete(z_mean_sum)

   
  end do
  

  z_total_mean=dim_sum_n(z_total, 0)
  z_total_mean!0="lat"
  z_total_mean!1="lon"
  z_total_mean&lat=lat
  z_total_mean&lon=lon
  z_total_mean&lat@units="degree_north"
  z_total_mean&lon@units="degree_east"
  
  totalres=res
  totalres@gsnCenterString="total"
  totalres@gsnCenterStringFontHeightF=0.02
  totalres@tmXBOn=True
  totalres@tmYLOn=True
  totalres@tmYLMode="ExplicitLevels"
  totalres@tmYLValues=(/35,40/)
  totalres@tmYLLabels=(/"35N","40N"/)
  totalres@tmYLMinorValues=ispan(31,43,1)


  totalres@tmYLLabelFontHeightF=0.008
  totalres@tmXBLabelFontHeightF=0.008


  totalres@vpHeightF=0.45
  totalres@vpWidthF=0.45
  totalres@vpXF=0.045
  totalres@vpYF=0.748
  totalres@lbLabelBarOn=True
  totalres@pmLabelBarHeightF=0.08
  totalres@pmLabelBarOrthogonalPosF=0.06
  totalres@lbTitleOn        =  True                ; turn on title
  totalres@lbTitleString    = "Mg/month/grid"              ; title string
  totalres@lbTitlePosition  = "Bottom"              ; title position
  totalres@lbTitleFontHeightF= .01              ; make title smaller
  totalres@lbTitleDirection = "Across"             ; title direction

   plot_total=gsn_csm_contour_map(wks, z_total_mean, totalres)

     ;;;µ╖╗σèáσî║σƒƒ
   ;;indian
   pgx_indian=(/115.8,117,117,115.8,115.8/)
   pgy_indian=(/39.4,39.4,40.5,40.5,39.4/)
   ;;;

   str1=unique_string("dum1")
   pgres=True
   pgres@gsLineColor="red"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   plot@$str1$ = gsn_add_polyline(wks, plot_total, pgx_indian, pgy_indian, pgres)





   
 ;;;panel
  pres=True
  pres@gsnPanelLabelBar=False
  pres@gsnFrame=False
  ; pres@pmLabelBarHeightF=0.05
  ; pres@pmLabelBarOrthogonalPosF=0.15
  ; pres@pmLabelBarParallelPosF=0.5
  pres@lbTitleOn        =  True                ; turn on title
  pres@lbTitleString    = "Mg/month/grid"              ; title string
  pres@lbTitlePosition  = "Right"              ; title position
  pres@lbTitleFontHeightF= .01              ; make title smaller
  pres@lbTitleDirection = "Across"             ; title direction
  pres@gsnPanelLeft=0.5
  pres@gsnPanelDebug=True



  gsn_panel(wks, plot, (/2,2/), pres)
  
    ; drawNDCGrid(wks)
draw(plot_total)
frame(wks)



  end 