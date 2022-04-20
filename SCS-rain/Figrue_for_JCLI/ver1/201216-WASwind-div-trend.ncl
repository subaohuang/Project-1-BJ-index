  ;;;;
;==============================================================
; Mission 4. plot(3) 
; ==============================================================
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
    ;;;region 
      lat1 = -5
      lat2 = 45
      lon1 = 90
      lon2 = 170
    ;;;
    yearnum = 33
   
    f3 = addfile("/home/yangsong3/data-observation/linshh/data/WASwind/waswind_v1_0_1.monthly.nc", "r")

    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.201112)
    time = f3->time(timeind)
    lat = f3->lat
    lon = f3->lon


    u = f3->u(timeind,{lat1:lat2},{lon1:lon2})
    v = f3->v(timeind,{lat1:lat2},{lon1:lon2})

    div = uv2dv_cfd(u,v,u&lat,u&lon,2)
    
    copy_VarCoords(u, div)
  printVarSummary(div)
    div = (/div*10^6/)
    div = (/-1*div/)
   
     w = div
season = (/"MAM","JJA","SON","DJF"/)
str = ispan(1,4,1) +""
wks = gsn_open_wks("pdf", "waswind-div-trend")

plot = new(4,graphic)
iplot = new(4,graphic)
do i = 0,3

    w_annual = month_to_season(w, season(i))

    if season(i).eq."DJF" then 
        w_annual(0:yearnum-2,:,:) = w_annual(1:yearnum-1,:,:)
        w_annual(yearnum-1,:,:) = w(yearnum*12-1,:,:)
    end if
    ; yearnum = 37
    wgt_w_rc := regCoef_n(ispan(1,yearnum,1), w_annual, 0, 0)
    copy_VarCoords(w_annual(1,:,:), wgt_w_rc)

    prob_w = rc_test(wgt_w_rc)
    wgt_w_rc = (/wgt_w_rc*yearnum/) 

    res4=True
    res4@gsnFrame=False
    res4@gsnDraw=False
    res4@gsnAddCyclic = False
  
    res4@cnFillOn=True
    res4@cnLinesOn=False

    res4@cnLineLabelsOn=False
    
    if i.eq.2 then 
    res4@lbLabelBarOn=True
    else 
    res4@lbLabelBarOn=False
    end if
    res4@lbLabelStride=1
    res4@lbOrientation="horizontal"
    res4@pmLabelBarOrthogonalPosF = 0.17


    res4@mpLimitMode           = "LatLon"
    res4@mpMinLatF             = -5
    res4@mpMaxLatF             = 30
    res4@mpMinLonF             = 100
    res4@mpMaxLonF             = 170 
    res4@mpCenterLonF          = (res4@mpMinLonF+res4@mpMaxLonF)/2.0
    res4@tmXTOn=False
    res4@tmYROn=False
    res4@tmXBOn=True
    res4@tmYROn=False
    res4@tmXBMode="Explicit"
    res4@tmYLMode="Explicit"


    res4@tmXBLabelFontHeightF=0.025
    res4@tmYLLabelFontHeightF=0.025
    res4@tmXBMajorThicknessF=1.0
    res4@tmYLMajorThicknessF=1.0
    res4@tmBorderThicknessF=1.5
    res4@tmXBValues = (/100,120,140,160/) 
    res4@tmXBLabels = (/100,120,140,160/) +""+"E" 

    res4@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
    res4@cnFillPalette = "MPL_PuOr"
    ; gsn_define_colormap(wks2, "MPL_PuOr")
    
    res4@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
    res4@lbLabelStrings = res4@cnLevels
    res4@lbLabelFontHeightF = 0.018
    ; res4@cnFillColors = (/22,32,45,60,-1,80,97,100,110/) -4

    
    res4@gsnLeftString=" "+season(i)+" low-level wind convergence trend "
    res4@gsnLeftStringFontHeightF = 0.03
    plot(i) = gsn_csm_contour_map(wks, wgt_w_rc, res4)

    
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=0
    ;;;画box
    lat1 = 0
    lat2 = 15
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)


  ;;;打点 ;;;打点
    res_pro4=True
    res_pro4@gsnAddCyclic = False
    res_pro4@gsnDraw=False
    res_pro4@gsnFrame=False
    res_pro4@cnFillOn = True

    res_pro4@cnLevelSelectionMode = "ExplicitLevels"
    res_pro4@cnLevels = (/0.0,0.1/)
    res_pro4@lbLabelBarOn = False

    res_pro4@cnLinesOn=False
    ; res_pro4@cnMonoLineColor =False
    ; res_pro4@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro4@cnMonoFillPattern = True
    res_pro4@cnFillPattern = 17
    res_pro4@cnFillDotSizeF = 0.001
    res_pro4@cnMonoFillColor = False
    res_pro4@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro4@cnLevels = ispan(-6,6,1)
    res_pro4@cnLineLabelsOn=False
    res_pro4@tmYRMode = "Automatic"
    res_pro4@cnInfoLabelOn=False



    iplot(i)=gsn_csm_contour(wks,prob_w,res_pro4)

    overlay(plot(i), iplot(i))
end do


pres = True

gsn_panel(wks, plot, (/2,2/), pres)

end