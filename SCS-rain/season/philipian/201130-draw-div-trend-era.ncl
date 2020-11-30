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

    season = (/"MAM","JJA","SON","DJF"/)
    plot = new(4,graphic)

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
    ;


    div = dim_avg_n_Wrap(div_r(:,:,{-10:30},{100:180}),1)

    div = (/div*10^5/)
    div = (/-1*div/)
    wks   = gsn_open_wks("pdf", "low-level-divergence-trend.pdf")

do i = 0,3  
        div_annual = month_to_season(div, season(i))
        yearnum = dimsizes(div_annual&time)

        div_rc := regCoef_n(ispan(1, yearnum, 1), div_annual, 0, 0)
        copy_VarCoords(div_annual(1,:,:), div_rc)
        div_rc = (/div_rc*yearnum/)

        prob_div = rc_test(div_rc)
        ; div_rc = where(prob_div.le.0.1, div_rc, div_rc@_FillValue)


    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic  = False

    res@cnFillPalette="NCV_blu_red"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=1
    res@pmLabelBarOrthogonalPosF = 0.17
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 30
    res@mpMinLonF             =100
    res@mpMaxLonF             = 170
    res@lbOrientation="horizontal"
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpFillOn = False
 

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025
    res@gsnLeftStringFontHeightF = 0.03
    res@tmYLMinorOn = True

    res@tmXBValues = (/100,120,140,160/) 
    res@tmXBLabels = (/100,120,140,160/) +""+"E" 
    ; res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levelW
  res@cnMinLevelValF       =   -0.1
  res@cnMaxLevelValF       =   .1
  res@cnLevelSpacingF      =  0.01
    res@cnFillPalette = "precip_diff_12lev"


    ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ; ; res@cnLevels = (/-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0,2.5/)
    ; ; res@cnFillColors= (/4,5,6,7,8,9,10,11,12,13,14/) - 2

    ; res@cnLevels = (/-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/5,6,7,8,10,11,12,13,14/) - 2

    res@lbLabelFontHeightF = 0.02
    res@lbAutoManage = False

    res@gsnLeftString=season(i) + " divergence trend " 
    res@gsnLeftStringOrthogonalPosF  = -0.12

    plot(i) = gsn_csm_contour_map(wks, div_rc, res)

end do

pres = True

gsn_panel(wks, plot, (/2,2/), pres)


end