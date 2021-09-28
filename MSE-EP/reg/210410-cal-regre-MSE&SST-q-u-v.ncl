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

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  copy_VarCoords(data(1,:,:), rc1)
  
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

begin
  yearstart = 197901
  yearend = 201912

  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360

    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")

    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    MSE = f1->MSE(timeind,:,{lat1:lat2},{lon1:lon2})

  ;;;MSE 
  plev = MSE&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

  MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)

  MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

  MSE_col = (/dtrend_msg_n(ispan(1,dimsizes(MSE_col&time),1),MSE_col,False,True,0)/)

  ; MSE_col = (/runave_n_Wrap(MSE_col, 9*12, 0, 0)/)
  MSE_col = (/MSE_col/10^6/)

  MSE_col_annual = month_to_annual(MSE_col, 1)
  MSE_col_annual!0 = "time"

  MSE_col_MAM = month_to_season(MSE_col, "MAM")
  MSE_col_JJA = month_to_season(MSE_col, "JJA")
  MSE_col_SON = month_to_season(MSE_col, "SON")
  MSE_col_DJF = month_to_season(MSE_col, "DJF")

  MSE_nep_mam = dim_avg_n_Wrap(MSE_col_MAM(:,{0:20},{210:280}), (/1,2/))

  
  ;;;读入HadiSST海温数据
  diri="/home/yangsong3/data-observation/linshh/data/sst/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
  ; fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

  f1=addfile(fils1, "r")
  date:=cd_calendar(f1->time, 1)

  timeind:=ind(date.le.yearend.and.date.ge.yearstart)

  sst  = lonFlip(f1->sst(timeind,{lat1:lat2},:))
  
  ssta=(/rmMonAnnCycTLL(sst)/)
  copy_VarCoords(sst, ssta)
  
  ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
;
  ; ssta = (/runave_n_Wrap(ssta, 9*12, 0, 0)/)

  ssta_mam = month_to_season(ssta, "MAM")

  ssta_djf = month_to_season(ssta, "DJF")

  ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:)  = ssta_djf(1:dimsizes(ssta_djf&time)-1,:,:)
  ssta_djf(dimsizes(ssta_djf&time)-1,:,:)  = ssta(dimsizes(ssta&time)-1,:,:)
  ;;;
  

  ;;;dsst 
  dssta = center_finite_diff_n(ssta, 1, False, 1, 0)
  copy_VarCoords(ssta, dssta)

  ; dssta = (/runave_n_Wrap(dssta, 9*12, 0, 0)/)
  ;;;

    ;;;read data
    path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    var_name = (/"uwnd","vwnd","slp"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
    
    ;;;region
    lat1 = -50
    lat2 = 50
    lon1 = 0
    lon2 = 360

    f1 = addfile(filein(0),"r")

    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.yearstart.and.date.le.yearend)

    u =  f1->uwnd(timeind,{1000},{lat1:lat2},{lon1:lon2})

    f2 = addfile(filein(1),"r")
    v = f2->vwnd(timeind,{1000},{lat1:lat2},{lon1:lon2})

    ; f3 = addfile(filein(2),"r")
    ; date :=cd_calendar(f3->time, 1)
    ; timeind :=ind(date.ge.startmon.and.date.le.endmon)
    ; slp = f3->slp(timeind,{lat1:lat2},{lon1:lon2})      
    ; printVarSummary(slp)

    ; slp  = (/slp*100.0/)
    u =(/rmMonAnnCycTLL(u)/)
    v =(/rmMonAnnCycTLL(v)/)
    u = (/dtrend_msg_n(ispan(1,dimsizes(u&time),1),u,False,True,0)/)
    v = (/dtrend_msg_n(ispan(1,dimsizes(u&time),1),v,False,True,0)/)

;;;


;;;read GPCP data 
    f4 = addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")

    ; f2=addfile("/home/yangsong3/data-observation/CPC_monthly_precip/precip.mon.mean.nc", "r")
    date:=cd_calendar(f4->time, 1)
    timeind:=ind(date.ge.yearstart.and.date.le.yearend)
    p=f4->precip(timeind,:,:)

    p = (/rmMonAnnCycTLL(p)/)

    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
  


    season = (/"MAM","JJA","SON","DJF"/)
    plot = new(4,graphic)
    plot2 = new(4,graphic)
    iplot=new(4,graphic)

    MSE_nep_mam = (/dim_standardize_n_Wrap(MSE_nep_mam, 1, 0)/)
    MSE_nep_mam = (/dtrend_msg_n(ispan(1,dimsizes(MSE_nep_mam&time),1),MSE_nep_mam,False,True,0)/)


  ; rc_mse = regCoef_n(ispan(1,dimsizes(MSE_nep_mam),1), MSE_nep_mam, 0, 0)
  
  ; rc_mse = rc_mse*dimsizes(MSE_nep_mam)
; wks   = gsn_open_wks("pdf", "DJF SST异常回归到春季MSE")
  wks   = gsn_open_wks("pdf", "ERSST-dSST异常回归到春季MSE-"+yearstart+"-"+yearend)


do i=0,3
   dssta_mam = month_to_season(dssta, season(i))
   u_mam = month_to_season(u, season(i))
   v_mam = month_to_season(v, season(i))
   p_mam = month_to_season(p, season(i))


   if season(i).eq."DJF" then
    dssta_mam(0:dimsizes(dssta_mam&time)-2,:,:)  = (/dssta_mam(1:dimsizes(dssta_mam&time)-1,:,:)/)
    dssta_mam(dimsizes(dssta_mam&time)-1,:,:)  = (/dssta(dimsizes(dssta&time)-1,:,:)/)
    
    u_mam(0:dimsizes(u_mam&time)-2,:,:)  = (/u_mam(1:dimsizes(u_mam&time)-1,:,:)/)
    u_mam(dimsizes(u_mam&time)-1,:,:)  = (/u(dimsizes(u&time)-1,:,:)/)
    
    v_mam(0:dimsizes(v_mam&time)-2,:,:)  = (/v_mam(1:dimsizes(v_mam&time)-1,:,:)/)
    v_mam(dimsizes(v_mam&time)-1,:,:)  = (/v(dimsizes(v&time)-1,:,:)/)
    
    p_mam(0:dimsizes(p_mam&time)-2,:,:)  = (/p_mam(1:dimsizes(p_mam&time)-1,:,:)/)
    p_mam(dimsizes(p_mam&time)-1,:,:)  = (/p(dimsizes(p&time)-1,:,:)/)

   end if
    
   
    rc_mam := regCoef_n(MSE_nep_mam, dssta_mam, 0, 0)
    copy_VarCoords(dssta_mam(0,:,:), rc_mam)

    rc_u_mam := regCoef_n(MSE_nep_mam, u_mam, 0, 0)
    copy_VarCoords(u_mam(0,:,:), rc_u_mam)


    rc_v_mam := regCoef_n(MSE_nep_mam, v_mam, 0, 0)
    copy_VarCoords(v_mam(0,:,:), rc_v_mam)

    rc_p_mam := regCoef_n(MSE_nep_mam, p_mam, 0, 0)
    copy_VarCoords(p_mam(0,:,:), rc_p_mam)


    ; rc_mam = (/rc_mam*rc_mse/)

    prob_mam := rc_test(rc_mam)
    prob_u_mam := rc_test(rc_u_mam)
    prob_v_mam := rc_test(rc_v_mam)
    prob_p_mam := rc_test(rc_p_mam)

    rc_u_mam = where(prob_u_mam.le.0.1.or.prob_v_mam.le.0.1, rc_u_mam, rc_u_mam@_FillValue)
    rc_v_mam = where(prob_u_mam.le.0.1.or.prob_v_mam.le.0.1, rc_v_mam, rc_v_mam@_FillValue)
 
;;;;
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  ; res@lbLabelBarOn=False

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 150
  res@mpMaxLonF             = 290
  

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -0.1
  ; res@cnMaxLevelValF       =   0.1
  ; res@cnLevelSpacingF      =  0.02


   res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  res@cnLevels := (/-8,-6,-4,-2,2,4,6,8/)*0.02
  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_blu_red"

  res@gsnLeftString=" Regression of "+season(i)+" SST tendency & 1000UV  onto MAM MSE  " 
  ; res@gsnLeftString=  "The components of "+""+season(i)+""+" ssta tendency trend linked to MAM MSE  " 

 
  ;;;;;;;LineArrow
  res@gsnScalarContour        = True
  res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; res@vcGlyphStyle="CurlyVector"

  res@vcMinDistanceF           = 0.02           ; thin vectors
  ; res@vcMinMagnitudeF          =0.04
  res@vcRefMagnitudeF          = 1.0
  res@vcRefLengthF             = 0.08     ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF=1.5
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  ; res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcVectorDrawOrder="PostDraw"
  res@mpOutlineOn = True
  res@mpCountyLineColor = "black"
  res@mpFillOn =False
  ; res@mpFillColor = "Transparent"
  ; res@mpFillDrawOrder = "PreDraw"

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
  plot(i) = gsn_csm_vector_scalar_map(wks,rc_u_mam,rc_v_mam,rc_mam,res)
  ; print(rc_p_mam({-20:20},{150:290}))
  ; exit()
  res2 = res
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnFillOn = False
  res2@cnLinesOn = True
  res2@cnLevels := ispan(4, 20, 4)*0.1
  res2@cnLineLabelStrings =   res2@cnLevels  + ""
  res2@cnLineLabelDensityF = 1.3

  res2@cnLineColor = "green"
  res2@gsnContourNegLineDashPattern = 1
  
  res2@cnMonoLineColor  = True
  res2@cnLineThicknessF = 2
  res2@cnMonoLineThickness  = True
  res2@gsnLeftString = ""
  res2@cnInfoLabelOn = False
  ; res2@cnLineLabelsOn = True
  res2@cnLineLabelInterval = 1
  ; res2@cnLineLabelBackgroundColor = "White"
  
   plot2(i) = gsn_csm_contour(wks, rc_p_mam, res2)

  overlay(plot(i), plot2(i))
   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 
   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

  
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 210
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

     str2 = ispan(21, 24, 1) +""
     str3 = ispan(25, 28, 1) +""

     ; plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

   lat1 = 0
   lat2 = -10
   lon1 = 270
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    pgres@gsLineColor="red"


   plot@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)    
 
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    ; i=0
  ;    ; do i=0,3
  ;   iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
  ;    overlay(plot(i), iplot(i))
  ; ; end do/

    ;;;  
end do 
pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01
gsn_panel(wks, (/plot(0),plot(1)/), (/2,1/), pres)

end


