;;********************************************************
;;画散点图,得到关系
;;
;;lsh 2020-11-04
;**********************************************************
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end




begin

startmon = 198201
endmon = 201912

; startmon = 197901
; endmon = 201512



 ;==============================================================
; Mission 3. calculate OLR
; ==============================================================
    f6 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
    ;;;read data
    lat1 = -40
    lat2 =  40
    lon1 =  0
    lon2 = 360


 
    date:=cd_calendar(f6->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    olr = short2flt(f6->olr(timeind,{lat1:lat2},{lon1:lon2}))

    plot = new(4,graphic)
    iplot = new(4,graphic)
   
   yearnum = 2019 -1982+1
  season = (/"MAM","JJA","SON","DJF"/)
  wks = gsn_open_wks("pdf", "OLR-trend-fourseason")
 do i=0,3
    olr_annual = month_to_season(olr, season(i))
    olr_rc := regCoef_n(ispan(1, yearnum, 1), olr_annual, 0, 0)
    copy_VarCoords(olr_annual(1,:,:),olr_rc)
    olr_rc = (/olr_rc*yearnum/)
    prob_olr_rc = rc_test(olr_rc)


    ;;;; draw olr 
    res3=True
    res3@gsnFrame=False
    res3@gsnDraw=False
    res3@gsnAddCyclic = False

    res3@cnFillOn=True
    res3@cnLinesOn= False
    

    res3@cnLineLabelsOn=False
  
    res3@lbLabelBarOn=False

    res3@lbLabelStride=1
    res3@lbOrientation="horizontal"
    res3@pmLabelBarOrthogonalPosF = 0.17

   
    res3@mpLimitMode           = "LatLon"
    res3@mpMinLatF             = lat1
    res3@mpMaxLatF             = lat2
    res3@mpMinLonF             = lon1
    res3@mpMaxLonF             = lon2

    ; res3@tmXBValues = (/100,120,140,160/) 
    ; res3@tmXBLabels = (/100,120,140,160/) +""+"E" 

    res3@mpCenterLonF          = (res3@mpMinLonF+res3@mpMaxLonF)/2.0

    res3@tmXTOn=False
    res3@tmYROn=False
    res3@tmXBOn=True
    res3@tmYROn=False

    res3@tmXBLabelFontHeightF=0.025
    res3@tmYLLabelFontHeightF=0.025
    res3@tmXBMajorThicknessF=1.0
    res3@tmYLMajorThicknessF=1.0
    res3@tmBorderThicknessF=1.5

    res3@gsnScalarContour        = True
    res3@mpFillOn =False

    res3@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res3@cnMinLevelValF        := -12.0        ; set min contour level
    res3@cnMaxLevelValF        := 12.0          ; set max contour level
    res3@cnLevelSpacingF       := 2.0          ; set contour interval
    ; res3@cnFillPalette = "BlueWhiteOrangeRed"
    res3@cnFillPalette = "BlueRed"
    ; res3@cnLevelSelectionMode = "ExplicitLevels"

    ; gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    ; res3@cnLevels := (/-20.0,-15.0,-10.0,-5.0,-2.5,2.5,5.0,10.0,15.0,20.0/)
  


    ; res3@lbLabelStrings = (/-20.0,-10.0,10.0,20.0/)


    ; res3@cnFillColors := (/2,5,30,40,3,146,177,200,251/) -4

    
    res3@gsnLeftString=" "+season(i)+" OLR trend " 
    res3@gsnLeftStringFontHeightF = 0.03
    res3@lbAutoManage = False


    plot(i) = gsn_csm_contour_map(wks,olr_rc,res3)


  ;;;打点
    ;;;打点
    res_pro2=True
    res_pro2@gsnAddCyclic = False
    res_pro2@gsnDraw=False
    res_pro2@gsnFrame=False
    res_pro2@cnFillOn = True

    res_pro2@cnLevelSelectionMode = "ExplicitLevels"
    res_pro2@cnLevels = (/0,0.1/)
    res_pro2@lbLabelBarOn = False


    res_pro2@cnLinesOn=False


    res_pro2@cnMonoFillPattern = True
    res_pro2@cnFillPattern = 17
    ; res_pro2@cnFillDotSizeF = 0.0045
    res_pro2@cnFillDotSizeF = 0.002
    res_pro2@cnMonoFillColor = False

  
    res_pro2@cnLineLabelsOn=False
    res_pro2@tmYRMode = "Automatic"
    res_pro2@cnInfoLabelOn=False
    res_pro2@cnFillColors = (/"Transparent","White","Transparent"/)

    iplot(i)=gsn_csm_contour(wks,prob_olr_rc,res_pro2)

    overlay(plot(i), iplot(i))
  end do
  ;;;;

pres = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/2,2/), pres)

end 