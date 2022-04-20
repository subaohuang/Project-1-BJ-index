;;********************************************************
;;计算ERA的散度趋势以及稳定度趋势
;;
;;lnshuheng 2020.03.17
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

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 
;==========================================================
;;;回归u，v,sst，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    var_name = (/"uwnd","vwnd"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
     ;;;region
     lat1 = -50
     lat2 = 50
     lon1 = 0
     lon2 = 360

    ;;;time span
    startmon =198001
    endmon = 202012

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    u =  f1->uwnd(timeind,{850},{lat1:lat2},{lon1:lon2})

    u = (/rmMonAnnCycTLL(u)/)

    f2 = addfile(filein(1),"r")
    v = f2->vwnd(timeind,{850},{lat1:lat2},{lon1:lon2})
    v =  (/rmMonAnnCycTLL(v)/)
  ;;;;;


;;;;
 
   ;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst  = lonFlip(f1->sst(timeind,:,:))
    ; sst  = short2flt(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
   ;;
;;;;;

;;;;;;split to 12 month  is a array size (12,yearnum)
  
     u_12 = split_12(u)
     v_12 = split_12(v)
     sst_12 = split_12(ssta)
;;;;;


;;;choose season
    plot  = new(1,graphic)
    iplot = new(1,graphic)

   ;;;处理成年数据 
     u_annual = month_to_season(u, "MAM")
     v_annual = month_to_season(v, "MAM")
     sst_annual = month_to_season(sst, "MAM") 

     u_annual = (/dim_avg_n(u_12(0:4,:,:,:),0)/)
     v_annual = (/dim_avg_n(v_12(0:4,:,:,:),0)/)
     sst_annual = (/dim_avg_n(sst_12(0:4,:,:,:),0)/)

   ;;;

      str = "" + ispan(1,4,1) + ""
      wks = gsn_open_wks("pdf", "u,v,SST趋势-NCEP-NCAR")

     ;;趋势
     u_rc := regCoef_n(ispan(1, yearnum, 1), u_annual, 0, 0)
     v_rc := regCoef_n(ispan(1, yearnum, 1), v_annual, 0, 0)
     sst_rc := regCoef_n(ispan(1, yearnum, 1), sst_annual, 0, 0)

     copy_VarCoords(u_annual(1,:,:), u_rc)
     copy_VarCoords(v_annual(1,:,:), v_rc)
     copy_VarCoords(sst_annual(1,:,:), sst_rc)

     sst_rc = (/sst_rc*yearnum/)
     u_rc = (/u_rc*yearnum/)
     v_rc = (/v_rc*yearnum/)
  


   ;;;显著性检验
    
    prob_u = rc_test(u_rc)
    prob_v = rc_test(v_rc)
    prob_sst = rc_test(sst_rc)
   
   ;;;
 
   ; d_rc = where(prob_d.le.0.1, d_rc, d_rc@_FillValue)
   ; cape_rc = where(prob_cape_rc.le.0.1, cape_rc, cape_rc@_FillValue)
   ; sta_s_rc = where(prob_sta_s.le.0.1, sta_s_rc, sta_s_rc@_FillValue)

   u_rc_test = u_rc
   v_rc_test = v_rc
   u_rc_test = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

;==============================================================
; Mission 4. draw plot
; ==============================================================

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
 
  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueRed"

  res@cnFillOn=True
  res@cnLinesOn= False
  
  res@cnLineThicknessF = 2.0
  res@gsnContourNegLineDashPattern = 14
  res@gsnContourPosLineDashPattern = 0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
  
  ; ; res@cnFillPalette = "CBR_drywet"
  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; res@cnLevels := (/-200,-150,-100,-50,50,100,150,200/)

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  ; res@lbOrientation="Vertical"
  res@pmLabelBarOrthogonalPosF = 0.2
  res@lbLabelStride = 2 
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             = 125
  ; res@mpMaxLonF             = 160
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  ;  ;   lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160
  res@mpFillOn = True
  res@mpFillColor = "gray"
  ; res@mpMonoFillColor = True
  

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.5
  
  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.02           ; thin vectors
  res@vcMinMagnitudeF          =0.6
  res@vcRefMagnitudeF          = 6.0
  res@vcRefLengthF             = 0.05     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcVectorDrawOrder="PostDraw"
  res@mpOutlineOn = True
  res@mpCountyLineColor = "black"
  ; res@mpFillOn =False
  ; res@mpFillColor = "Transparent"
  ; res@mpFillDrawOrder = "PreDraw"

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
  res@vcRefAnnoString1On = True
  res@vcRefAnnoString1 = "6.0 m/s"
  res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
  res@vcRefAnnoParallelPosF = 1.02
  res@vcRefAnnoFontHeightF = 0.01

  res@vcRefAnnoBackgroundColor = "Transparent"
  res@vcRefAnnoPerimOn = False

  res@gsnLeftStringOrthogonalPosF  = -0.12
  res@gsnLeftString  = "JFMAM SST&850hPa wind trend"

  plot  = gsn_csm_vector_scalar_map(wks,u_rc,v_rc,sst_rc ,res)

  ;;;;;;;LineArrow
  resv = True
  resv@gsnDraw = False
  resv@gsnFrame = False
  resv@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  resv@vcMinDistanceF           = 0.022           ; thin vectors
  res@vcMinMagnitudeF          =0.4
  resv@vcRefMagnitudeF          = 4.0
  resv@vcRefLengthF             = 0.1     ; define length of v
  resv@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  resv@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  resv@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  resv@vcLineArrowHeadMinSizeF = 0.006
  resv@vcVectorDrawOrder="PostDraw"

  resv@vcMonoLineArrowColor = True
  resv@vcLineArrowColor = "green"
  resv@vcLineArrowThicknessF=1.5

  resv@vcRefAnnoOn    = False
  resv@gsnLeftString         = " "
  resv@gsnRightString        = " "




  ; plot_vect := gsn_csm_vector(wks,u_rc_test({-5:30},{100:170}),v_rc_test({-5:30},{100:170}),resv)
  ; overlay (plot(i),plot_vect)

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -2.0
  ; res@cnMaxLevelValF       =   8.0
  ; res@cnLevelSpacingF      =  2.0
 
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 250
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   i = 0
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
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = gsn_csm_contour(wks,prob_sst,res_pro)
  
    overlay(plot, iplot)
draw(plot)
frame(wks)



end
    
  







    
   