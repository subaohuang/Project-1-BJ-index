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
; function cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
;   begin
;   ;;;;Q component
 
;     Q_mk  = trend_manken(era_divQ_delta, False, 0)
;     Q_rc = Q_mk(1,:,:)
;     copy_VarCoords(era_divQ_delta(1,:,:), Q_rc)
;     Q_rc = (/Q_rc*yearnum/)
;     Q_rc@prob = Q_mk(0,:,:)
    
;     era_wgt_qu_delta_mk  = trend_manken(era_wgt_qu_delta, False, 0)
;     era_wgt_qu_delta_rc = era_wgt_qu_delta_mk(1,:,:)
;      copy_VarCoords(era_wgt_qu_delta(1,:,:), era_wgt_qu_delta_rc)

;     era_wgt_qu_delta_rc = (/era_wgt_qu_delta_rc*yearnum/)
;     era_wgt_qu_delta_rc@prob = era_wgt_qu_delta_mk(0,:,:)

;     era_wgt_qv_delta_mk  = trend_manken(era_wgt_qv_delta, False, 0)
;     era_wgt_qv_delta_rc = era_wgt_qv_delta_mk(1,:,:)
;     copy_VarCoords(era_wgt_qv_delta(1,:,:), era_wgt_qv_delta_rc)
;     era_wgt_qv_delta_rc = (/era_wgt_qv_delta_rc*yearnum/)
;     era_wgt_qv_delta_rc@prob = era_wgt_qv_delta_mk(0,:,:)

;  return [/Q_rc,era_wgt_qu_delta_rc,era_wgt_qv_delta_rc/]
;  end

begin


  ;;;read era data
        path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
        var_name = (/"uwnd","vwnd","shum/shum.total"/)
        filein = path_wind + var_name +".mon.mean.nc"
        
        ;;;region 
          lat1 = -10
          lat2 = 50
          lon1 = 60
          lon2 = 180
        ;;;

        ; ;;;global region 
        ; lat1 = -50
        ; lat2 = 50
        ; lon1 = 0
        ; lon2 = 360
        ; ;;;

  
        ;;;time span
        startmon =197901
        endmon = 201512

        startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
        endyear_r = str_split_by_length(endmon, 4 )

        startyear = stringtoint(startyear_r(0))
        endyear = stringtoint(endyear_r(0))

        yearnum = endyear - startyear +1
        ;;


        f1 = addfile(filein(0),"r")

        date=cd_calendar(f1->time, 1)
        timeind=ind(date.ge.startmon.and.date.le.endmon)

        u = short2flt(f1->u(timeind,:,{lat1:lat2},{lon1:lon2}))

        f2 = addfile(filein(1),"r")
        v = short2flt(f2->v(timeind,:,{lat1:lat2},{lon1:lon2}))

        
        f3 = addfile(filein(2),"r")
        date :=cd_calendar(f3->time, 1)
        timeind :=ind(date.ge.startmon.and.date.le.endmon)
        q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
   ;;;
   print("=====finish read data==========")



;==============================================================
; Mission 1. calculate Q
; ==============================================================
    era_qu  = q
    era_qv  = q
    era_qu  = u*q
    printVarSummary(era_qu)
    era_qv  = v*q
    era_wgt_qu = era_qu(:,0,:,:)
    era_wgt_qv = era_qv(:,0,:,:)


   ;;;------Cal DivQ-------
    ;-----------Follow is the vertical integral-----------
    pver=era_qu&level
    pver=pver*100
    pver@units="Pa"

    dp=dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(era_qu),dp,1)          ; confirm dp to caculate sum(var*dp)

    era_wgt_qu=(/dim_sum_n(dp_cfm*era_qu,1)/)               ; sum(var*dp)
    era_wgt_qv=(/dim_sum_n(dp_cfm*era_qv,1)/)               ; sum(var*dp)

    era_wgt_qu = (/era_wgt_qu/9.8/)
    era_wgt_qv = (/era_wgt_qv/9.8/)

    g_lon   =   era_wgt_qu&longitude
    g_lat   =   era_wgt_qu&latitude
    g_nlat  =   dimsizes(g_lat)

    DivQ = uv2dv_cfd(era_wgt_qu,era_wgt_qv,g_lat,g_lon,2)
    ; DivQ = uv2dvF_Wrap(era_wgt_qu(:,::-1,:),era_wgt_qv(:,::-1,:)) 
    DivQ   =   (/-DivQ*86400/) 
    copy_VarCoords(era_wgt_qu, DivQ)
    printVarSummary(DivQ)
    ; print(DivQ(0,:,:))
    print("=====finish DivQ calculate==========")
;;;;




 print("=====choose season==========")
;;;choose season
   season = (/"MAM","JJA","SON","DJF"/)
   plot  = new(4,graphic)
   iplot = new(4,graphic)
   plot_tp = plot
   add_a = (/"(a)","(b)","(c)","(d)"/)

   str = "" + ispan(1,4,1) + ""
    wks   = gsn_open_wks("pdf", "q_delta趋势以及箭头—四季-phipilian-curlwind")

do i = 0,3    
      era_divQ_bar = month_to_season(DivQ, season(i))
       if season(i).eq."DJF"
          era_divQ_bar(0:dimsizes(era_divQ_bar&time)-2,:,:) = era_divQ_bar(1:,:,:) 
          era_divQ_bar(dimsizes(era_divQ_bar&time)-1,:,:) = DivQ(dimsizes(DivQ&time)-1,:,:)
        end if 
      era_divQ_delta = dim_rmvmean_n_Wrap(era_divQ_bar, 0)

      printVarSummary(era_divQ_delta)

      ;;;得到divQ_bar的两个箭头
      era_wgt_qu_annual = month_to_season(era_wgt_qu, season(i))
        if season(i).eq."DJF"
          era_wgt_qu_annual(0:dimsizes(era_wgt_qu_annual&time)-2,:,:) = era_wgt_qu_annual(1:,:,:) 
          era_wgt_qu_annual(dimsizes(era_wgt_qu_annual&time)-1,:,:) = era_wgt_qu(dimsizes(era_wgt_qu&time)-1,:,:)
        end if 
      era_wgt_qu_delta =  dim_rmvmean_n_Wrap(era_wgt_qu_annual, 0)


      era_wgt_qv_annual = month_to_season(era_wgt_qv, season(i))
       if season(i).eq."DJF"
          era_wgt_qv_annual(0:dimsizes(era_wgt_qv_annual&time)-2,:,:) = era_wgt_qv_annual(1:,:,:) 
          era_wgt_qv_annual(dimsizes(era_wgt_qv_annual&time)-1,:,:) = era_wgt_qv(dimsizes(era_wgt_qv&time)-1,:,:)
        end if 
      era_wgt_qv_delta =  dim_rmvmean_n_Wrap(era_wgt_qv_annual, 0)

      era_wgt_qu_delta = (/era_wgt_qu_delta*9.8/)
      era_wgt_qv_delta = (/era_wgt_qv_delta*9.8/)
      ;;;

; ;==============================================================
; ; Mission 3.  trend
; ; ==============================================================
    ;;;;Q component
    Q_rc =  cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
    
    q_rc = Q_rc[0]
    qu_rc  =  Q_rc[1]
    qv_rc  =  Q_rc[2]

    prob_q = rc_test(q_rc)
    prob_qu = rc_test(qu_rc)
    prob_qv = rc_test(qv_rc)

    ; prob_q = 1- q_rc@prob
    ; prob_qu = 1- qu_rc@prob
    ; prob_qv = 1- qv_rc@prob

    copy_VarCoords(q_rc, prob_q)
    copy_VarCoords(q_rc, prob_qv)
    copy_VarCoords(q_rc, prob_qu)
    

    qu_rc_test = qu_rc
    qv_rc_test = qv_rc
    qu_rc_test = where(prob_qu.le.0.1.or.prob_qv.le.0.1, qu_rc,qu_rc@_FillValue)

;==============================================================
; Mission 6. draw plot
; ==============================================================
  colordata = read_colormap_file("/home/ys17-19/lsh/data/colordata/bluewhiteorangered_11color.rgb")

  ; gsn_define_colormap(wks, colordata)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette=colordata(::-1,:)

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -3.0
  res@cnMaxLevelValF       =   3.0
  res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2

  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             = 90
  ; res@mpMaxLonF             = 160

    ;;;global region
   res@mpMinLatF             = -5
   ; res@mpMaxLatF             = 50
   res@mpMaxLatF             = 30
   res@mpMinLonF             = 100
   res@mpMaxLonF             = 170
  ; ;;;

 
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.025
  res@tmYLLabelFontHeightF=0.025
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@gsnLeftStringFontHeightF = 0.03
  res@tmXBValues = (/100,120,140,160/) 
  res@tmXBLabels = (/100,120,140,160/) +""+"E" 



  res@mpFillOn = False
  res@gsnScalarContour        = True
    ;;;;;;;LineArrow
  ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  res@vcGlyphStyle="CurlyVector"
   res@vcMonoLineArrowColor = True
  res@vcLineArrowColor = "black"
  res@vcMinDistanceF           = 0.025           ; thin vectors
  
  res@vcRefLengthF             = 0.08     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcVectorDrawOrder="PostDraw"

  res@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  res@vcRefMagnitudeF          = 1000.0
  ; res@vcMinMagnitudeF          =100.0
  res@gsnLeftString = add_a(i)+" "+season(i) +" Integrated Moisture Flux Trend"
  
  plot(i)=gsn_csm_vector_scalar_map(wks,qu_rc,qv_rc, q_rc, res)


    resv                             = True 
    resv@gsnDraw                     = False
    resv@gsnFrame                    = False      ;don't advance frame
    ; resv@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcGlyphStyle="CurlyVector"

    resv@vcMinDistanceF           = 0.025           ; thin vectors

    resv@vcRefLengthF             = 0.08     ; define length of v
    resv@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ; resv@vcVectorDrawOrder="PostDraw"
    resv@vcMonoLineArrowColor = True
    resv@vcLineArrowColor = "blue"
    resv@vcLineArrowThicknessF=1.5
    ; resv@vcLineArrowHeadMaxSizeF=0.01
    ; resv@vcLineArrowHeadMinSizeF=0.01
    resv@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    resv@vcLineArrowHeadMinSizeF = 0.006

    resv@vcRefMagnitudeF          = 1000.0
    ; resv@vcMinMagnitudeF          =100.0

    resv@vcRefAnnoOn    = False
    resv@gsnLeftString         = " "
    resv@gsnRightString        = " "


    plot_vect := gsn_csm_vector(wks,qu_rc_test({-5:30},{100:170}),qv_rc_test({-5:30},{100:170}),resv)
    overlay (plot(i),plot_vect)
 
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
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
    ; res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot(i)=gsn_csm_contour(wks,prob_q,res_pro)    
    overlay(plot(i), iplot(i))

    
    plot_tp(i) = overlay_tp_horizontal(wks,-10,50,60,160,"gray")
    overlay(plot(i), plot_tp(i))

end do

 pres = True
 pres@gsnPanelLabelBar = True
 ; pres@lbOrientation = "vertical"
 pres@lbLabelFontHeightF = 0.02
 ; pres@gsnPanelFigureStrings = (/"a","b","c","d"/)
 ; pres@gsnPanelFigureStringsJust = "TopLeft"

 gsn_panel(wks, plot, (/2,2/), pres)
 

end 
    
   
   

  





    
   ;;;

   






