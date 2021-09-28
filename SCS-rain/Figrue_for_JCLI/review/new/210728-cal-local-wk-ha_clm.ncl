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
    diri="/home/yangsong3/data-observation/linshh/data/era/"
    fils1=systemfunc("ls "+diri+"era_m_zonal_meri_ERAI_500hpa.nc")

    startmon = 197901
    endmon = 201512

    f1=addfile(fils1, "r")
    time = f1->time 
    date = cd_calendar(time, 1)

    timeind = ind(date.ge.startmon.and.date.le.endmon)

    m_zonal=f1->m_zonal(timeind,:,:)
    m_meri=f1->m_meri(timeind,:,:)

    m_zonal = (/m_zonal*10^3/)
    m_meri = (/m_meri*10^3/)
    
    ; m_zonal_annual=month_to_annual(m_zonal_500, 1)
    ; m_meri_annual=month_to_annual(m_meri_500, 1)


    season = (/"MAM","JJA","SON","DJF"/)

    title = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/)

    wks=gsn_open_wks("pdf", "local-HC&WK-clm")
    ; wks1=gsn_open_wks("pdf", "local-hc-clm")

    plot=new(8,graphic)
    ; plot1=new(4,graphic)
    str = ispan(1,8,1)+""
    str1 = ispan(1,8,1)+"a"
    str2 = ispan(1,8,1)+"b"

do i=0,7
;;;;
    if i.le.3 then 
     m_zonal_annual=month_to_season(m_zonal, season(i))
     m_meri_annual=month_to_season(m_meri, season(i))
    else 
      m_zonal_annual=month_to_season(m_zonal, season(i-4))
      m_meri_annual=month_to_season(m_meri, season(i-4))
    end if
    dims=dimsizes(m_zonal_annual)
;;;;;

 ;;;trend 
    m_zonal_rc = dim_avg_n_Wrap(m_zonal_annual, 0)
    m_meri_rc = dim_avg_n_Wrap(m_meri_annual, 0)
 ;;;;



;;;画图 
 
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True

    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn =False

    res@gsnLeftStringFontHeightF = 0.02


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -6.0
    res@cnMaxLevelValF       =   6.0
    res@cnLevelSpacingF      =  1.0

    res@cnLineLabelsOn = False

    res@lbLabelBarOn = False

    res@cnFillPalette = "BlRe"
  
    res@mpLimitMode           = "LatLon"

    ; res@mpMinLatF             = -5
    ; res@mpMaxLatF             = 30
    ; res@mpMinLonF             = 100
    ; res@mpMaxLonF             = 170

    res@mpMinLatF             = -60
    res@mpMaxLatF             = 60
    res@mpMinLonF             = 0
    res@mpMaxLonF             = 360

    res@mpCenterLonF = (/res@mpMinLonF + res@mpMaxLonF /)*0.5


    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF=0.015
    res@tmYLLabelFontHeightF=0.015
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0


    res@tmXBValues =ispan(0,360,30)
    res@tmXBLabels = (/"0","30E","60E","90E","120E","150E","180","150W","120W","90W","60W","30W","0"/)
    res@tmXBMinorValues = ispan(0,360,10)

    res@tmYLValues = ispan(-60,60,30)
    res@tmYLLabels = (/"60S","30S","0","30N","60N"/)
    res@tmYLMinorValues = ispan(-60,60,10)

    res@tmXTOn=False
    res@tmYROn=False

    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.9
    ; res@vpHeightF = 0.2 
    ;;; *************************************************

  
    res@gsnRightString = ""
   

    if i.ge.4 then 
      res@gsnLeftString = title(i)+" "+season(i-4)  +" local WC climatology"
      plot(i)= gsn_csm_contour_map(wks,m_zonal_rc,res)
     else 
      res@gsnLeftString = title(i)+" "+season(i)  +" local HC climatology"
      plot(i)= gsn_csm_contour_map(wks,m_meri_rc,res)
    end if



  ; ;;添加区域
   pgres = True
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
  
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

   
   ; if i.le.3 then 
   ; pgx = (/125,125/)
   ; pgy = (/60,-60/)
   ; plot@$str1(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy , pgres)
   ; pgx = (/160,160/)
   ; pgy = (/60,-60/)
   ; plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy , pgres)
   ; else 
   ;  pgx = (/125,125/)
   ;  pgy = (/60,-60/)
   ;  plot@$str1(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy , pgres)
   ;  pgx = (/160,160/)
   ;  pgy = (/60,-60/)
   ;  plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy , pgres)


 
  delete(m_zonal_rc)
  delete(m_meri_rc)
  
end do
  pres=True
  pres@gsnPanelLabelBar = True
  pres@lbLabelFontHeightF = 0.02
  pres@pmLabelBarOrthogonalPosF = -0.02
  gsn_panel(wks, plot, (/4,2/), pres)


end 

