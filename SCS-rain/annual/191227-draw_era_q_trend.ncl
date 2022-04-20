;;********************************************************
;;计算水汽通量算度，并且分解水汽通量散度
;; ERA-interim:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
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
  ;;;read era data
		path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
		var_name = (/"shum"/)
		filein = path_wind + var_name +".mon.mean.nc"
        
    ;;;region 
    lat1 = -10
    lat2 = 50
    lon1 = 60
    lon2 = 160
    ;;;

    ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;

    f3 = addfile(filein(0),"r")
    date :=cd_calendar(f3->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
    ;;;
   print("=====finish read data==========")

   ;;;计算趋势
    
    q_bar = month_to_annual(q, 1)
    
    
    pver = q_bar&level
    pver = pver*100
    pver@units="Pa"

    dp = dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm = conform_dims(dimsizes(q_bar),dp,1)          ; confirm dp to caculate sum(var*dp)

    q_wgt =  q_bar(:,1,:,:)   
    q_wgt = (/dim_sum_n(dp_cfm*q_bar,1)/)               ; sum(var*dp)
    q_wgt_clm = dim_avg_n_Wrap(q_wgt, 0)

    q_rc  = regCoef_n(ispan(1, yearnum, 1), q_wgt, 0, 0)
    copy_VarCoords(q_wgt(1,:,:), q_rc)
    prob1 = rc_test(q_rc)
    q_rc = (/q_rc*yearnum/)

    q_percent = q_rc
    q_percent = (/100*q_rc/q_wgt_clm/)

    wks   = gsn_open_wks("pdf", "q的趋势")
    plot=new(1,graphic)
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic = False

    res@cnFillPalette="MPL_RdBu"
    res@cnFillOn=True
    res@cnLinesOn=False
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -8.0
    res@cnMaxLevelValF       =   8.0
    res@cnLevelSpacingF      =  2
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=2
    res@lbOrientation="Vertical"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 45
    res@mpMinLonF             = 90
    res@mpMaxLonF             = 160
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False

    res@tmXBLabelFontHeightF=0.015
    res@tmYLLabelFontHeightF=0.015
    res@tmXBMajorThicknessF=3.0
    res@tmYLMajorThicknessF=3.0
    res@tmBorderThicknessF=4.0


    plot(0)=gsn_csm_contour_map(wks, q_percent, res)

    str1 = unique_string("dum1")  ; "unique_string" will return a unique    

    ; ;;添加区域
    ;;warm pool
    pgres = True
    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=0

    str = (/"1"/)

    ;;;画box
    lat1 = 0
    lat2 = 15
    lon1 = 110
    lon2 = 155

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    
    plot@$str(0)$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)
    

    ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False

    ; res_pro@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res_pro@cnMinLevelValF       =   0
    ; res_pro@cnMaxLevelValF       =   0.05
    ; res_pro@cnLevelSpacingF      =  0.05

    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(1,graphic)

    iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
    

    i=0
    overlay(plot(i), iplot(i))
    
   draw(plot)
   frame(wks)
   end
    



