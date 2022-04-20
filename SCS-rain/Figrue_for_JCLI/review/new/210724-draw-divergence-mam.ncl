;==============================================================
; Mission read  wind data
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

function vertical_avg(data)
begin
    psfc =1013.0
    plev = data&level/100.0
    
    ptop = min(plev)
    dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)

    dp_4d = conform_dims(dimsizes(data), dp, 1)
    
    q1dp = data(:,1,:,:)
    ; q1dp =  (/dim_sum_n_Wrap(data*dp_4d, 1)/)
    q1dp =  dim_sum_n_Wrap(data*dp_4d, 1)/dim_sum_n(dp_4d, 1)   ;;;积分
return q1dp
end 

begin
 startmon = 197901
 endmon = 201512

  ;;;read era data
  path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
  var_name = (/"uwnd","vwnd","shum/shum.total"/)
  filein = path_wind + var_name +".mon.mean.nc"

  
  ;;;region 
    lat1 = -10
    lat2 = 30
    lon1 = 90
    lon2 = 180
  ;;;


    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; timeind :=ind(date.ge.197901.and.date.le.201512)

    yearnum = 37

    u = short2flt(f1->u(timeind,:,{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(1),"r")
    v = short2flt(f1->v(timeind,:,{lat1:lat2},{lon1:lon2}))


    f3 = addfile(filein(2),"r")
    date :=cd_calendar(f3->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    timeind :=ind(date.ge.197901.and.date.le.201512)
    q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
    ;;;
print("=====finish read data==========")
;;;;


;;;;;
   u_mean = vertical_avg(u(:,{850:100},{-10:30},{100:180}))
   v_mean = vertical_avg(v(:,{850:100},{-10:30},{100:180}))
;;;



;;;计算散度 
    div_r = u
    div_r = uv2dv_cfd(u,v,u&latitude,u&longitude, 1)
    
    vor_r = u
    vor =  uv2vr_cfd(u,v,u&latitude,u&longitude, 1)
;;;
 



   div = vertical_avg(div_r(:,{850:100},{-10:30},{100:180}))
   ; div = vor_r(:,{850},{-10:30},{100:180})

    ; div = div_r

    div = (/div*10^6/)
    div = (/-1*div/)
   
    w = div

    w_annual = month_to_season(w,"MAM")

    ; yearnum = 37
    wgt_w_rc := regCoef_n(ispan(1,yearnum,1), w_annual, 0, 0)
    copy_VarCoords(w_annual(1,:,:), wgt_w_rc)

    prob_w = rc_test(wgt_w_rc)
    wgt_w_rc = (/wgt_w_rc*yearnum/) 
;;;


; ;;;风场异常变化导致的辐合
;     u_annual = month_to_season(u_mean,"MAM")
;     wgt_u_rc := regCoef_n(ispan(1,yearnum,1), u_annual, 0, 0)
;     copy_VarCoords(u_annual(1,:,:), wgt_u_rc)

;     v_annual = month_to_season(v_mean,"MAM")
;     wgt_v_rc := regCoef_n(ispan(1,yearnum,1), v_annual, 0, 0)
;     copy_VarCoords(v_annual(1,:,:), wgt_v_rc)

;     printVarSummary(wgt_u_rc)
;     printVarSummary(u&latitude)
;     div_2 = wgt_u_rc
;     div_2 = uv2dv_cfd(wgt_u_rc, wgt_v_rc, wgt_u_rc&latitude, wgt_v_rc&longitude, 1)

;     div_2 = (/div_2*10^6/)
;     div_2 = (/-1*div_2*yearnum/)
 

; ;;;


    



wks = gsn_open_wks("pdf", "divergence-850-100")
;;;;;
    res4=True
    res4@gsnFrame=False
    res4@gsnDraw=False
    res4@gsnAddCyclic = False
  
    res4@cnFillOn=True
    res4@cnLinesOn=False

    res4@cnLineLabelsOn=False
    res4@lbLabelBarOn=True
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
    res4@tmXBMajorThicknessF=2.0
    res4@tmYLMajorThicknessF=2.0
    ; res4@tmBorderThicknessF=1.0
    res4@tmXBValues = (/100,120,140,160/) 
    res4@tmXBLabels = (/100,120,140,160/) +""+"E" 
    res4@tmYLValues = (/0,10,20,30/) 
    res4@tmYLLabels = (/"0","10N","20N","30N"/) 

    res4@tmYLMinorValues = ispan(-5,30,5)

    res4@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
    res4@cnFillPalette = "MPL_PuOr"
    
    res4@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
    res4@lbLabelStrings = res4@cnLevels
    res4@lbLabelFontHeightF = 0.018

    
    res4@gsnLeftString = "(c) MAM low-level wind convergence trend "
    res4@gsnLeftStringFontHeightF = 0.03
    plot = gsn_csm_contour_map(wks, wgt_w_rc, res4)
    ; plot = gsn_csm_contour_map(wks, div_2, res4)




    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=1
    ;;;画box
    lat1 = 0
    lat2 = 20
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    ;;;
    str = "2"
    plot@$str$ = gsn_add_polyline(wks, plot, pgx_pool, pgy_pool, pgres)



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



    iplot = gsn_csm_contour(wks,prob_w,res_pro4)

    overlay(plot, iplot)
  ;;;;

draw(plot)
frame(wks)
end



