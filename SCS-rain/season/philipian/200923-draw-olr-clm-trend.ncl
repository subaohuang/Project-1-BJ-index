begin
   ;;;;
   startmon = 197901
   endmon = 201612


   ;;;region
   lat1 = -40
   lat2 = 40
   lon1 = 40
   lon2 = 280
   ;;;;


   f6 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
   date:=cd_calendar(f6->time, 1)
   timeind :=ind(date.ge.startmon.and.date.le.endmon)
   olr = short2flt(f6->olr(timeind,{lat1:lat2},{lon1:lon2}))
   ;;;

 

    wks = gsn_open_wks("pdf", "OLR-clm")



    season = (/"MAM","JJA","SON","DJF"/)
    plot  = new(4,graphic)
    iplot = new(4,graphic)


do i = 0,3 


  olr_annual = month_to_season(olr, season(i))
  
   olr_clm = dim_avg_n_Wrap(olr_annual, 0)
;==============================================================
; Mission 4. draw plot
; ==============================================================

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  ; res@cnFillPalette ="CBR_coldhot"
  ; res@cnFillPalette = "cmp_b2r"

  color = read_colormap_file("MPL_PuBu")
  res@cnFillPalette = color(::-1,:)



  res@cnFillOn=True
  res@cnLinesOn= False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -180.0
  ; res@cnMaxLevelValF       =   180.0
  ; res@cnLevelSpacingF      =  30.0
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   210.0
  res@cnMaxLevelValF       =   290.0
  res@cnLevelSpacingF      =  10.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@lbLabelStride = 2 
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             = 125
  ; res@mpMaxLonF             = 160
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2
  ;  ;   lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.5
  
  res@gsnScalarContour        = True


  res@gsnLeftString  = season(i)+"OLR CLM"
  plot(i) = gsn_csm_contour_map(wks, olr_clm, res)
end do 

pres = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)

end 
 




  


   