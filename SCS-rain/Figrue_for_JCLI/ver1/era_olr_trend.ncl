
  
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
 ;==============================================================
; Mission 3. calculate OLR
; ==============================================================

    fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/toa_netthermal/toa_netthermal_monthly_*")
    f3 = addfiles(fils3, "r")

    startmon =197901
    endmon = 201512

      startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1

    date=cd_calendar(f3[:]->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude
    lon = f3[1]->longitude
 
    ;;;region 
    lat1 = -90
    lat2 = 90
    lon1 = 0
    lon2 = 360
    ;;;


    olr = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)

    do i = 0,dimsizes(timeind)/12-1

    olr(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->ttr(:,:,:))
    end do


    olr = (/olr/(24*3600)/) 

    season = (/"MAM","JJA","SON","DJF"/)

    plot = new(4,graphic)
    iplot = plot

      str = ispan(1,4,1)+""

    wks   = gsn_open_wks("pdf", "era-OLR-trend-79-15")


 do i=0,3
    olr_annual = month_to_season(olr, season(i))


    if season(i).eq."DJF" then 
    olr_annual(0:yearnum-2,:,:) = olr_annual(1:yearnum-1,:,:)
    olr_annual(yearnum-1,:,:) = olr(yearnum*12-1,:,:)
    end if

    olr_rc := regCoef_n(ispan(1, yearnum, 1), olr_annual, 0, 0)
    copy_VarCoords(olr_annual(1,:,:),olr_rc)
    olr_rc = (/-1*olr_rc*yearnum/)
    prob_olr_rc = rc_test(olr_rc)


    ;;;; draw olr 
    res3=True
    res3@gsnFrame=False
    res3@gsnDraw=False
    res3@gsnAddCyclic = False

    res3@cnFillOn=True
    res3@cnLinesOn= False
    

    res3@cnLineLabelsOn=False
    if i.eq.2 then 
    res3@lbLabelBarOn=True
    else 
    res3@lbLabelBarOn=False
    end if
    res3@lbLabelStride=1
    res3@lbOrientation="horizontal"
    res3@pmLabelBarOrthogonalPosF = 0.17

   
    res3@mpLimitMode           = "LatLon"
    res3@mpMinLatF             = -40
    res3@mpMaxLatF             = 40
    res3@mpMinLonF             = 30
    res3@mpMaxLonF             = 290

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
    res3@cnFillPalette = "BlueWhiteOrangeRed"
    ; res3@cnFillPalette = "BlueRed"
    ; res3@cnLevelSelectionMode = "ExplicitLevels"

    ; gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    ; res3@cnLevels := (/-20.0,-15.0,-10.0,-5.0,5.0,10.0,15.0,20.0/)
    ; res3@lbLabelStrings = res3@cnLevels
    res3@lbLabelFontHeightF = 0.02
    res3@lbAutoManage = False


    ; res3@lbLabelStrings = (/-20.0,-10.0,10.0,20.0/)


    ; res3@cnFillColors := (/2,5,30,40,3,146,177,200,251/) -4

    
    res3@gsnLeftString=season(i)+" OLR trend " 
    res3@gsnLeftStringFontHeightF = 0.03
    res3@lbAutoManage = False


    plot(i) = gsn_csm_contour_map(wks,olr_rc,res3)

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.5
    pgres@gsLineDashPattern=0
    ;;;画box
    lat1 = 0
    lat2 = 15
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

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


pres= True
gsn_panel(wks, plot, (/2,2/), pres)

end
  ;;;;