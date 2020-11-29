;;********************************************************
;;计算水汽通量算度，并且分解水汽通量散度
;; ERA-interim:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
;**********************************************************
load  "/home/ys17-19/lsh/Project/other/guanshun/add_TP.ncl"
  
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

function cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
begin
  ;;;;Q component
    Q_rc  = regCoef_n(ispan(1, yearnum, 1), era_divQ_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), Q_rc)
    Q_rc = (/Q_rc*yearnum/)
    
    era_wgt_qu_delta_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_qu_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), era_wgt_qu_delta_rc)
    era_wgt_qu_delta_rc = (/era_wgt_qu_delta_rc*yearnum/)

    era_wgt_qv_delta_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_qv_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), era_wgt_qv_delta_rc)
    era_wgt_qv_delta_rc = (/era_wgt_qv_delta_rc*yearnum/)
 return [/Q_rc,era_wgt_qu_delta_rc,era_wgt_qv_delta_rc/]
 end

;==============================================================
; Mission 0. plot(0) 
; ==============================================================
begin
   ;;;choose season
   season = (/"JJA","SON","DJF"/)
   order =(/"(a)","(b)","(c)","(d)"\
          , "(e)","(f)","(g)","(h)","(i)"/)
   ;;;


    plot = new(9,graphic)
    iplot = plot

    str = ispan(1,9,1)+""

    wks   = gsn_open_wks("pdf", "fig4-2015-925hpa.pdf")

    ;;;time span
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;


    ; f1=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    f1=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201512)
    p=f1->precip(timeind,:,:)
   
    delete(p&time)
    p&time=date(timeind)
  
    p_pi = p

    delete(p)

 ;;;;;
;==============================================================
; Mission read  wind data
; ==============================================================
  ;;;read era data
  path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
  var_name = (/"uwnd","vwnd"/)
  filein = path_wind + var_name +".mon.mean.nc"

  ;;;region 
    lat1 = -10
    lat2 = 50
    lon1 = 60
    lon2 = 180
  ;;;


    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; timeind :=ind(date.ge.197901.and.date.le.201512)

    yearnum = 37

    u = short2flt(f1->u(timeind,{925},{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(1),"r")
    v = short2flt(f1->v(timeind,{925},{lat1:lat2},{lon1:lon2}))
 
    ;;;
print("=====finish read data==========")
;;;;


 
   do i=0,2

    p_annual = month_to_season(p_pi, season(i))

    if season(i).eq."DJF" then 
     p_annual(0:yearnum-2,:,:) = p_annual(1:yearnum-1,:,:)
     p_annual(yearnum-1,:,:) = p_pi(yearnum*12-1,:,:)
    end if

    
    rc_mam  := regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
    copy_VarCoords(p_annual(1,:,:), rc_mam)
    rc_mam = (/rc_mam*yearnum/)

    prob_mam=rc_test(rc_mam)
    copy_VarCoords(p_pi(1,:,:), prob_mam)
    ;;;
   
   ;;;calculate wind trend 
    u_annual = month_to_season(u(:,:,:), season(i))
    v_annual = month_to_season(v(:,:,:), season(i))

    ; if season(i).eq."DJF" then 
    ;  u_annual(0:yearnum-2,:,:) = u_annual(1:yearnum-1,:,:)
    ;  u_annual(yearnum-1,:,:) = u(yearnum*12-1,{850},:,:)
    ; end if

    ; if season(i).eq."DJF" then 
    ;  v_annual(0:yearnum-2,:,:) = v_annual(1:yearnum-1,:,:)
    ;  v_annual(yearnum-1,:,:) = v(yearnum*12-1,{850},:,:)
    ; end if

    u_rc  = regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
    copy_VarCoords(u_annual(1,:,:), u_rc)
    u_rc = (/u_rc*yearnum/)
    
    v_rc  = regCoef_n(ispan(1, yearnum, 1),v_annual, 0, 0)
    copy_VarCoords(u_annual(1,:,:), v_rc)
    v_rc = (/v_rc*yearnum/)

    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)
    

    ; prob_u
    if season(i).eq."JJA" then
    prob_u({10:20},{100:120}) = where(prob_u({10:20},{100:120}).le.0.4, 0.15, prob_u@_FillValue)   
    u_rc = where(prob_u.le.0.15.or.prob_v.le.0.15, u_rc, u_rc@_FillValue) 
    else
    u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)
    end if




    ; v_rc = where(prob_u.le.01.or.prob_v.le.0.1, v_rc, v_rc@_FillValue)
    u_rc_test = u_rc
    v_rc_test = v_rc


    u_rc_test = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

    printVarSummary(u_rc)

    delete(u_annual)
    delete(v_annual)   




    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic = False

    res@cnFillPalette="NCV_blu_red"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    

    if i.eq.2 then 
    res@lbLabelBarOn=True
    else 
    res@lbLabelBarOn=False
    end if

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
    res@tmYLMajorThicknessF = 0.01
    res@tmXBMajorThicknessF = 0.01
    ; res@tmYLMinorOn = True

    res@tmXBValues = (/100,120,140,160/) 
    res@tmXBLabels = (/100,120,140,160/) +""+"E" 
    res@tmBorderThicknessF=1.5
    res@cnFillPalette = "precip_diff_12lev"
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    ; res@cnLevels = (/-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/4,5,6,7,8,9,10,11,12,13,14/) - 2

    res@cnLevels = (/-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    res@cnFillColors= (/5,6,7,8,10,11,12,13,14/) - 2


    res@lbLabelFontHeightF = 0.02
    res@gsnLeftStringOrthogonalPosF  = -0.12

    res@lbAutoManage = False

    res@gsnLeftString=order(i*3)+" "+season(i)+" Precip trend " 


    res@gsnScalarContour        = True

    ;;;;;;;LineArrow
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcMinDistanceF           = 0.022          ; thin vectors
    
    ; res@vcRefAnnoJust = "TopLeft"
    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False
    res@vcRefAnnoOrthogonalPosF = -1.2
    res@vcRefMagnitudeF          = 1.6
    res@vcMinMagnitudeF          = 0.16
    res@vcRefLengthF             = 0.06     ; define length of v
    res@vcLineArrowThicknessF=1.0
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006
    res@vcLineArrowHeadMaxSizeF = 0.018           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.012
    


    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "2 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013

    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False

    plot(i*3) = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,rc_mam ,res)

    ; plot(i*3) = gsn_csm_contour_map(wks, rc_mam, res)



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
    plot@$str(i*3)$ = gsn_add_polyline(wks, plot(i*3), pgx_pool, pgy_pool, pgres)


    ; ;;;打点
    ; res_pro=True
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLinesOn=False
    ; res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    ; res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    ; res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    ; res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False

    ; iplot(i*3) = gsn_csm_contour(wks,prob_mam,res_pro)    

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; ; opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeMid     = 17         ; stipple
    ; ; opt@gsnShadeFillDotSizeF =  0.0025
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    ; iplot(i*3) = gsn_contour_shade(iplot(i*3)  , 0,0.1,opt)
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


    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False

  
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)

    iplot(i*3)=gsn_csm_contour(wks,prob_mam,res_pro)

    overlay(plot(i*3) , iplot(i*3) )
end do 


 ;==============================================================
; Mission 3. calculate OLR
; ==============================================================
    f6 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
    ;;;region
    lat1 = -5
    lat2 = 45
    lon1 = 90
    lon2 = 170

    date:=cd_calendar(f6->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    olr = short2flt(f6->olr(timeind,{lat1:lat2},{lon1:lon2}))

 do i=0,2
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
    if i.eq.2 then 
    res3@lbLabelBarOn=True
    else 
    res3@lbLabelBarOn=False
    end if
    res3@lbLabelStride=1
    res3@lbOrientation="horizontal"
    res3@pmLabelBarOrthogonalPosF = 0.17

   
    res3@mpLimitMode           = "LatLon"
    res3@mpMinLatF             = -5
    res3@mpMaxLatF             = 30
    res3@mpMinLonF             = 100
    res3@mpMaxLonF             = 170

    res3@tmXBValues = (/100,120,140,160/) 
    res3@tmXBLabels = (/100,120,140,160/) +""+"E" 

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
    res3@cnMinLevelValF        := -20.0        ; set min contour level
    res3@cnMaxLevelValF        := 20.0          ; set max contour level
    res3@cnLevelSpacingF       := 5.0          ; set contour interval
    res3@cnFillPalette = "BlueWhiteOrangeRed"
    ; res3@cnLevelSelectionMode = "ExplicitLevels"

    ; ; gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    ; res3@cnLevels := (/-20.0,-15.0,-10.0,-5.0,5.0,10.0,15.0,20.0/)
    res3@lbLabelStrings = res3@cnLevels
    res3@lbLabelFontHeightF = 0.02
    res3@lbAutoManage = False


    ; res3@lbLabelStrings = (/-20.0,-10.0,10.0,20.0/)


    ; res3@cnFillColors := (/2,5,30,40,3,146,177,200,251/) -4

    
    res3@gsnLeftString=order(i*3+1)+" "+season(i)+" OLR trend " 
    res3@gsnLeftStringFontHeightF = 0.03
    res3@lbAutoManage = False


    plot(i*3+1) = gsn_csm_contour_map(wks,olr_rc,res3)


    plot@$str(i*3+1)$ = gsn_add_polyline(wks, plot(i*3+1), pgx_pool, pgy_pool, pgres)

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

    iplot(i*3+1)=gsn_csm_contour(wks,prob_olr_rc,res_pro2)

    overlay(plot(i*3+1), iplot(i*3+1))
end do
  ;;;;
;==============================================================
; Mission 4. plot(3) 
; ==============================================================

    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"

    var_name := (/"wwnd"/)
    filein := path_wind + var_name +".mon.mean.nc"
        
    ;;;region 
      lat1 = -5
      lat2 = 45
      lon1 = 90
      lon2 = 170
    ;;;

    f1 = addfile(filein(0),"r")

    date := cd_calendar(f1->time, 1)

    timeind:=ind(date.ge.197901.and.date.le.201512)

    w = short2flt(f1->w(timeind,:,{lat1:lat2},{lon1:lon2}))
    w = (/w*-1/)


do i = 0,2

    w_annual = month_to_season(w, season(i))


    wgt_w = w_annual(:,{500},:,:)

    yearnum = 37
    wgt_w_rc := regCoef_n(ispan(1,yearnum,1), wgt_w, 0, 0)
    copy_VarCoords(wgt_w(1,:,:), wgt_w_rc)

    prob_w = rc_test(wgt_w_rc)
    wgt_w_rc = (/wgt_w_rc*yearnum/) *100  

       
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
    
    res4@cnLevels = (/-0.02,-0.015,-0.01,-0.005,0.005,0.01,0.015,0.02/)*100
    res4@lbLabelStrings = res4@cnLevels
    ; res4@cnFillColors = (/22,32,45,60,-1,80,97,100,110/) -4

    
    res4@gsnLeftString=order(i*3+2)+" "+season(i)+" 500 hPa Omega trend "
    res4@gsnLeftStringFontHeightF = 0.03
    plot(i*3+2) = gsn_csm_contour_map(wks, wgt_w_rc, res4)


    plot@$str(i*3+2)$ = gsn_add_polyline(wks, plot(i*3+2), pgx_pool, pgy_pool, pgres)


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



    iplot(i*3+2)=gsn_csm_contour(wks,prob_w,res_pro4)

    overlay(plot(i*3+2), iplot(i*3+2))

end do
    delete(w) 
  ;;;;

;==============================================================
; Mission  panel  
; ==============================================================
pres = True

gsn_panel(wks, plot, (/3,3/), pres)


end