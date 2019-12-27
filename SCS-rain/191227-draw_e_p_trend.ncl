;;********************************************************
;;计算蒸发的趋势以及查看P-E的大小
;;oaflux的蒸发数据，文件的单位为cm/day，短型
;;以及GCCP的降水数据，
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
;****************************************
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
   ;;read  eva data 
    fils = systemfunc("ls /home/yangsong3/data-observation/OAflux/monthly/evaporation/evapr_oaflux_"+"*"+".nc")

    f1 = addfiles(fils, "r")
    date = cd_calendar(f1[:]->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    evp_r =short2flt(f1[:]->evapr(timeind,:,:))
    evp_r@missing_value = -9.96921e+36
    evp_r@_FillValue =   -9.96921e+36
    evp_r&lat@units = "degrees_north"
    evp_r&lon@units = "degrees_east"
    
    evp_r = where(evp_r.eq.3276.6, -9.96921e+36, evp_r)
    evp =  evp_r*10 / 365.0   ;;convert t0 day
    copy_VarMeta(evp_r, evp)
    delete(evp_r)
   
    printVarSummary(evp)
    
    
    ;;;处理数据
        ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

        yo = fspan(-88.75, 88.75, 72)
        xo = fspan(1.25, 358.75, 144)

        evp_new = linint2_Wrap(evp&lon, evp&lat, evp, True, xo, yo, 0)

        delete(evp_new&LAT)
        delete(evp_new&LON)
        evp_new!1 = "lat"
        evp_new!2 = "lon"
        evp_new&lat = yo
        evp_new&lon = xo
        printVarSummary(evp_new)
        ;;;

        ; ;;116.25存在插值后缺侧的现象，现在插值下
        ; m =  linmsg(u_av_1({100:120}),0)
        ; copy_VarCoords(u_av_1({100:120}), m)
        ; u_av_1({116.25}) = m({116.25})
        ; delete(m)
        ; ;;;
;;;;

;;;read pre data

    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.201312)
    p=f2->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
    printVarSummary(p)
    
;;;

;;;计算 趋势

    evp_annual = month_to_annual(evp_new, 1)
    p_annual  = month_to_annual(p, 1)
    


    p_e = p_annual
    p_e = (/p_annual - evp_annual/)
    copy_VarCoords(p_annual, p_e)
  
    dim = dimsizes(p_annual)

    p_rc = regCoef_n(ispan(1, dim(0), 1), p_annual, 0, 0)
    copy_VarCoords(p_annual(1,:,:), p_rc)
    prob1=rc_test(p_rc)
    p_rc=(/p_rc*dim(0)/)


    evp_rc = regCoef_n(ispan(1, dim(0), 1), evp_annual, 0, 0)
    copy_VarCoords(p_annual(1,:,:), evp_rc)

    prob2=rc_test(evp_rc)
    evp_rc=(/evp_rc*dim(0)/)
  
    p_e_rc = regCoef_n(ispan(1, dim(0), 1), p_e, 0, 0)
    copy_VarCoords(p_e(1,:,:), p_e_rc)
    prob3=rc_test(p_e_rc)
    p_e_rc=(/p_e_rc*dim(0)/)
   
   copy_VarCoords(p_annual(1,:,:), prob2)

;;;

;;**************************************************************
;画图
;;***************************************************************
    wks   = gsn_open_wks("pdf", "蒸发和降水的趋势")
    plot=new(3,graphic)
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnLeftStringFontHeightF=20.0
    res@cnFillPalette="NCV_blu_red"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=2
    res@lbOrientation="Vertical"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 45
    res@mpMinLonF             =90
    res@mpMaxLonF             = 160
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"


    res@tmXBLabelFontHeightF=0.015
    res@tmYLLabelFontHeightF=0.015
    res@tmXBMajorThicknessF=3.0
    res@tmYLMajorThicknessF=3.0
    res@tmBorderThicknessF=4.0

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -2.5
    res@cnMaxLevelValF       =   2.5
    res@cnLevelSpacingF      =  0.25

    res@cnFillPalette="MPL_BrBG"
    res@gsnLeftString="(a)GPCP Precipitation"

    plot(0) = gsn_csm_contour_map(wks, p_rc, res)
     res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.2
      res@gsnLeftString="(b)OAflux evaporation"
    plot(1) = gsn_csm_contour_map(wks, evp_rc, res)
     res@cnMinLevelValF       =   -2.5
    res@cnMaxLevelValF       =   2.5
    res@cnLevelSpacingF      =  0.25
     res@gsnLeftString="(c)P-E "
    plot(2) = gsn_csm_contour_map(wks, p_e_rc, res)

   

  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   do i=0,2
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False

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
    iplot = new(3,graphic)

    iplot(0)=gsn_csm_contour(wks,prob1,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob2,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob3,res_pro)
  

    do i=0,2
    overlay(plot(i), iplot(i))
    end do

 pres = True
 gsn_panel(wks, plot, (/3,1/), pres)


end