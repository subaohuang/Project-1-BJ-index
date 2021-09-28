;;;;;
;==============================================================
; Mission read  wind data
; ==============================================================
begin
  ;;;read era data
  path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
  var_name = (/"uwnd","vwnd","slp"/)
  filein = path_wind + var_name +".mon.mean.nc"

  ;;;region 
    lat1 = -90
    lat2 = 90
    lon1 = 0
    lon2 = 360
  ;;;
    

    startmon = 197901
    endmon = 201512

    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; timeind :=ind(date.ge.197901.and.date.le.201512)

    ; yearnum = 37
     ; yearnum  =  endyear -  startyear

    u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(1),"r")
    v = short2flt(f1->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f1 = addfile(filein(2),"r")
    slp = short2flt(f1->msl(timeind,{lat1:lat2},{lon1:lon2}))
    slp = slp/100
    printMinMax(slp, 1)
 
    ;;;
print("=====finish read data==========")
;;;;
    


  season = (/"MAM","JJA","SON","DJF"/)
 

  plot = new(12,graphic)
 
    wks   = gsn_open_wks("pdf", "era-850wind-clm-indo-pacific")

    str  = ispan(1, 12, 1) +"month"
do i =0,11

    ; u_clm = dim_avg_n_Wrap(month_to_season(u, season(i)), 0)
    ; v_clm = dim_avg_n_Wrap(month_to_season(v, season(i)), 0)


    u_clm = dim_avg_n_Wrap(u(i::12,:,:), 0)
    v_clm = dim_avg_n_Wrap(v(i::12,:,:), 0)
    slp_clm = dim_avg_n_Wrap(slp(i::12,:,:), 0)


;;;;draw plot 

       
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic  = False


    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -20
    res@mpMaxLatF             = 20
    res@mpMinLonF             =120
    res@mpMaxLonF             = 290
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
    res@gsnLeftStringFontHeightF = 0.02
    res@tmYLMinorOn = True

    res@cnFillOn=True
    res@cnLinesOn=False
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   1004
    res@cnMaxLevelValF       =   1020
    res@cnLevelSpacingF      =  2
    ;   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; res@cnMinLevelValF       =   -1
    ; res@cnMaxLevelValF       =   1
    ; res@cnLevelSpacingF      =  .1
    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
   
     res@gsnScalarContour        = True

    ; res@tmXBValues = (/100,120,140,160/) 
    ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 


    ;;;;;;;LineArrow
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcMinDistanceF           = 0.03           ; thin vectors
    
    ; res@vcRefAnnoJust = "TopLeft"
    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False
    ; res@vcRefAnnoOrthogonalPosF = -1.2
    res@vcRefMagnitudeF          = 10.0
    ; res@vcMinMagnitudeF          =0.2
    res@vcRefLengthF             = 0.08     ; define length of v
    res@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006
    res@vcLineArrowHeadMaxSizeF = 0.018           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.012
    


    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "10 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.3              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013

    ; res@vcRefAnnoBackgroundColor = "Transparent"
    ; res@vcRefAnnoPerimOn = False

    res@tiMainString = str(i)
    res@tiMainFontHeightF = 0.02
    res@tiMainOffsetYF = -0.008
    res@gsnRightString = ""
    res@gsnLeftString = ""

    ; plot(i) = gsn_csm_vector_map(wks,u_clm,v_clm ,res)
   plot(i) = gsn_csm_vector_scalar_map(wks,u_clm,v_clm,slp_clm ,res)

  end do 

  pres = True
  pres@gsnPanelLabelBar = True
  pres@pmLabelBarOrthogonalPosF = -0.02
  pres@lbLabelFontHeightF = 0.01

  gsn_panel(wks, plot, (/4,3/), pres)

  end 