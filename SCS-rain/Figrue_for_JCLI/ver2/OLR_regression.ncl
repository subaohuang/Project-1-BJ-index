;;********************************************************
;;计算SST倾向以及SST值与GPCP降水的趋势
;;其中热通量的方向都是向上为正，向下为负
;;的相关系数场，
;;NCEP-NCAR的精度为1.87×1.9
;;
;;                                         lnshuheng 2020.02.18
;**********************************************************
function chazhi(data)
 begin
   xi = data&lon
   yi = data&lat
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
   data_chazhi&lat = yo
   data_chazhi&lon = xo
   copy_VarAtts(data, data_chazhi)
 
   return data_chazhi
end 
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
undef("remove_enso_signal")
function remove_enso_signal(nino3_4,  data)
  begin  
   ;;;;;
      ; data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_de  = data
      data_rc = regCoef_n(nino3_4,data_de, 0, 0)
      copy_VarCoords(data(1,:,:), data_rc)
      data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
      m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
      nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))
  
      data_relate_enso  = nino3_4_3d*data_rc_3d + m 
      data_new = data 
      data_new = (/data - data_relate_enso/)
  return data_new
end 


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
  order  = (/"(a)","(b)","(c)","(d)"/)
 

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ;;;
    delete(date)
    delete(timeind)

    ssta_de = ssta
    ; ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_in = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
    printVarSummary(ssta_in)

    nino3_4 = month_to_season(ssta_in, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)


    ; nino3_4_new = nino3_4                ;;;;前期ENSO
    ; nino3_4_new(0:33) = (/nino3_4(1:)/)
    ; nino3_4_new(34) =ssta_index(35*12-1) 

    nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
;;;

;;;area mean 
      rad    = 4.0*atan(1.0)/180.0
      re     = 6371220.0
      rr     = re*rad

      dlon   = abs(lon(2)-lon(1))*rr
      dx     = dlon*cos(lat*rad)
      ;          
      jlat  = dimsizes( lat )
      dy     = new ( jlat, typeof(dx))      ; close enough
      
      dy(0)  = abs(lat(2)-lat(1))*rr
      dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5   
      dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr

      area   = dx*dy                                ; cell area function of latitude only
      copy_VarCoords(lat, area)
      printVarSummary(area)
;;;;

  

;;;read OLR data
    ; f6 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201701.nc", "r")
    ; ;;;region
    ; lat1 = -90
    ; lat2 = 90
    ; lon1 = 0
    ; lon2 = 360

    ; date:=cd_calendar(f6->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; olr = short2flt(f6->olr(timeind,{lat1:lat2},{lon1:lon2}))


    ; olr_chazhi = chazhi(olr)

    ; system("/bin/rm -f olr_chazhi_1979-2015.nc")
    ; ncout1 = addfile("olr_chazhi_1979-2015.nc", "c")
    ; ncout1->olr_chazhi = olr_chazhi

    f6 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/season/olr_chazhi_1979-2015.nc", "r")
    ;;;region
    lat1 = -90
    lat2 = 90
    lon1 = 0
    lon2 = 360

    date:=cd_calendar(f6->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    olr_chazhi = f6->olr_chazhi


;;;;
    wks   = gsn_open_wks("pdf", "OLR回归到WEP的SST的回归系数")
  
    plot=new(4,graphic)
    plot1 = plot

    str = (/"1","2","3","4"/)
 
    iplot = new(4,graphic)

    season =  (/"MAM","JJA","SON","DJF"/)
;;;;;



do i =0,3

  ;;;region
   lat1 = -30
   lat2 = 30
   lon1 = 30
   lon2 = 180
  ;;;

 
  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season(i))
  printVarSummary(ssta_annual)

  ;;;;remove_enso_signal 
  ssta_annual = (/remove_enso_signal(nino3_4, ssta_annual)/)
  ;;;;


  ssta_annual_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), ssta_annual, 0, 0)


  ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{-20:20},{40:100}), area({-20:20}), 1.0, 0) , 0, 0)
  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:15},{125:160}), area({0:15}), 1.0, 0) , 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum

  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)

 
  ;;;;;要分析的变量
  olr_annual = month_to_season(olr_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))

  ;;;;remove_enso_signal 
  olr_annual = (/remove_enso_signal(nino3_4, olr_annual)/)
 
  ;;;;;;;;


 

  ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:15},{125:160}), area({0:15}), 1.0, 0) 


  ;;;normalized before regression 
  ssta_index = (/dim_standardize_n(ssta_index, 1, 0)/)
  ;;;;




  ;;;correlation with sst
    r_sst_p := regCoef_n(ssta_index, olr_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_p)


    print(ssta_index_rc)

    ; r_sst_p = (/r_sst_p*ssta_index_rc/)
   ;;;;;
   




    print("*******显著性检验*****************")

    ;;;t test 
    prob_p1= rc_test(r_sst_p)

    printVarSummary(prob_p1)

    prob_p1 = where(prob_p1.eq.0, 2, prob_p1)
    ;;;;

    




  print("*********开始画图********")

;;**************************************************************
;画图
;;***************************************************************


  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@cnFillPalette="BlueWhiteOrangeRed"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = lat1
  ; res@mpMaxLatF             = lat2
  ; res@mpMinLonF             = lon1
  ; res@mpMaxLonF             = lon2
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             = 100
  ; res@mpMaxLonF             = 170

  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 290



  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.025
  res@tmYLLabelFontHeightF=0.025
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@gsnScalarContour        = True
  res@mpFillOn = False




  res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res@cnMinLevelValF      = -2.5       ; set min contour level
  res@cnMaxLevelValF      = 2.5      ; set max contour level
  res@cnLevelSpacingF     = 0.5       ; set contour spacing
  


  ; res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  ; res@cnMinLevelValF      = -1.0       ; set min contour level
  ; res@cnMaxLevelValF      = 1.0      ; set max contour level
  ; res@cnLevelSpacingF     = 0.1       ; set contour spacing


  ; res@cnFillPalette = "precip_diff_12lev"
  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
  ; res@cnLevels := (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)/2.0
  ; ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2



  ; res@vcGlyphStyle="CurlyVector"
  ; res@vcMonoLineArrowColor = True
  ; res@vcLineArrowColor = "black"
  ; res@vcMinDistanceF           = 0.025           ; thin vectors

  ; res@vcRefLengthF             = 0.08     ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  ; res@vcVectorDrawOrder="PostDraw"

  ; res@vcLineArrowThicknessF=1.5
  ; ; res@vcLineArrowHeadMaxSizeF=0.01
  ; ; res@vcLineArrowHeadMinSizeF=0.01
  ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  ; res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcRefMagnitudeF          = 0.4
  ; ; res@vcMinMagnitudeF          =100.0


  res@gsnLeftStringFontHeightF = 0.025
 
  res@gsnLeftString = order(i)+" "+ season(i)+ " linear congruency of OLR trend with TIOST"

  plot(i) = gsn_csm_contour_map(wks, r_sst_p, res)

  

   ; ;;;画气候态的等值线
   ;  res@gsnRightString=""
   ;  res@gsnLeftString=""
   ;  res@cnLevelSelectionMode = "ExplicitLevels"

   ;  ; res@cnLevelSelectionMode = "AutomaticLevels"
   ;  res@cnLinesOn=True
   ;  res@cnLevels := fspan(-0.6,0.2,13)
   ;  res@cnLineLabelsOn=True
   ;  res@cnLineLabelBackgroundColor="White"
   ;  res@cnLineLabelInterval = 1
   ;  res@cnLineThicknessF = 2.
   ;  res@cnFillOn = False
   ;  res@gsnContourNegLineDashPattern=1
   ;  res@gsnContourZeroLineThicknessF=1.0
   ;  res@gsnLeftString  = ""
   ;  res@gsnRightString =""
   ;  res@cnInfoLabelOn = False    

   ;  plot1(i) = gsn_csm_contour(wks, r_sst_w, res)
   ;  ;;;  
   ;  overlay(plot(i), plot1(i))
    
 

  ;;;添加indo-pacific warm pool 
    ; ;;添加区域
   ;;warm pool
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1

    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=0
    str3=unique_string("dum3")
    str4=unique_string("dum4")
 
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
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(i) = gsn_csm_contour(wks,prob_p1,res_pro)
   
     
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
   


    ;;;  
    
      overlay(plot(i), iplot(i))
   
      delete(lat1)
      delete(lat2)
      delete(lon1)
      delete(lon2)
end do

  pres = True
  pres@gsnPanelLabelBar  = True
  pres@lbLabelFontHeightF = 0.015

  gsn_panel(wks, plot, (/2,2/), pres)
  
 
 print("done")

end