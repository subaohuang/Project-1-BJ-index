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
    startmon = 197901
    endmon = 201512 

   

    fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
    f3 = addfiles(fils3, "r")

    date:=cd_calendar(f3[:]->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude
    lon = f3[1]->longitude


    div_r = new((/dimsizes(time),7, dimsizes(lat),dimsizes(lon)/),float)


    do i = 0,dimsizes(timeind)/12-1

       div_r(i*12:(i+1)*12-1,:,:,:) = short2flt(f3[i]->d(:,{1000:850},:,:))

    end do


    ;   div_r = new((/dimsizes(time), dimsizes(lat),dimsizes(lon)/),float)


    ; do i = 0,dimsizes(timeind)/12-1

    ;    div_r(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->d(:,{850},:,:))

    ; end do
    order = (/"(a)","(b)","(c)","(d)"/)
    season = (/"MAM","JJA","SON","DJF"/)
    str = season

    plot = new(4,graphic)
    iplot = plot

    div = dim_avg_n_Wrap(div_r(:,:,{-10:30},{100:180}),1)
    ; div = div_r
    div = (/div*10^6/)
    div = (/-1*div/)
    wks   = gsn_open_wks("pdf", "divergence-trend.pdf")
 

do i = 0,3  

    div_annual = month_to_season(div, season(i))

     if season(i).eq."DJF" then 
        div_annual(0:yearnum-2,:,:) = div_annual(1:yearnum-1,:,:)
        div_annual(yearnum-1,:,:) = div(yearnum*12-1,:,:)
    end if

    yearnum = dimsizes(div_annual&time)

    div_rc := regCoef_n(ispan(1, yearnum, 1), div_annual, 0, 0)
    copy_VarCoords(div_annual(1,:,:), div_rc)
    div_rc = (/div_rc*yearnum/)

    prob_div = rc_test(div_rc)
        ; div_rc = where(prob_div.le.0.1, div_rc, div_rc@_FillValue)


    res4=True
    res4@gsnFrame=False
    res4@gsnDraw=False
    res4@gsnAddCyclic = False
  
    res4@cnFillOn=True
    res4@cnLinesOn=False

    res4@cnLineLabelsOn=False
    
    res4@lbLabelBarOn=False

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


    res4@tmYLValues = (/0,10,20,30/) 
    res4@tmYLLabels = (/"0","10N","20N","30N"/) 
    res4@tmYLMinorValues = ispan(-5,30,5)

    res4@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
    res4@cnFillPalette = "MPL_PuOr"
    ; res4@cnFillPalette = "MPL_bwr"

    ; gsn_define_colormap(wks, "MPL_PuOr")
    
    res4@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
    res4@lbLabelStrings = res4@cnLevels
    res4@lbLabelFontHeightF = 0.018
    ; res4@cnFillColors = (/22,32,45,60,-1,80,97,100,110/) -4

    res4@gsnLeftString=order(i)+" "+season(i)+" low-level wind convergence trend "
    res4@gsnLeftStringFontHeightF = 0.03

    plot(i) = gsn_csm_contour_map(wks, div_rc, res4)

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
    iplot(i)=gsn_csm_contour(wks,prob_div,res_pro4)

    overlay(plot(i), iplot(i))


    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.5
    pgres@gsLineDashPattern=1
    ;;;画box
    lat1 = 0
    lat2 = 20
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)


end do

pres = True
pres@gsnPanelLabelBar = True
pres@lbLabelStride = 1
pres@lbLabelFontHeightF = 0.012

gsn_panel(wks, plot, (/2,2/), pres)


end