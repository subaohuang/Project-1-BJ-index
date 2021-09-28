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

function part_corr(var1,var2,var3)
begin 
  ;;;remove var3 influence 
  r12 = escorc_n(var1, var2, 0, 0)
  ; copy_VarCoords(var2(0,:,:), r12)  ;;;这是二维
  r13 = escorc_n(var1, var3, 0, 0)  ;;;这个是1dim
  r23 = escorc_n(var2, var3, 0, 0)   ;;;这个是2dim
  r12_3 = (/(r12 - r13*r23)/sqrt((1- r13^2)*(1- r23^2))/)
  rc_mam = r12
  rc_mam = (/r12_3/)
return rc_mam
end 


begin
    ;;;time span
    startmon =198001
    endmon = 202012


    f2=addfile("/home/yangsong3/data-observation/NCEP_NCAR/monthly/prate.sfc.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p_r=f2->prate(timeind,{-90:90},{0:360})
    p_r = (/p_r*86400/)
   
    ;;;处理数据
    ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

    yo = fspan(-90, 90, 73)
    xo = fspan(0.0, 357.5, 144)
    p = linint2_Wrap(p_r&lon, p_r&lat, p_r, True, xo, yo, 0)

    delete(p&LAT)
    delete(p&LON)
    p!1 = "lat"
    p!2 = "lon"
    p&lat = yo
    p&lon = xo
    p = (/rmMonAnnCycTLL(p)/)

   
;;;;
    f3=addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/evapr_oaflux_chazhi.nc", "r")
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    eva=f3->evapr(timeind,:,:)
    eva = (/eva*10/365/) ;;convert to mm/day
   


    eva_new = linint2_Wrap(eva&lon, eva&lat, eva, True, xo, yo, 0)

    delete(eva_new&LAT)
    delete(eva_new&LON)
    eva_new!1 = "lat"
    eva_new!2 = "lon"
    eva_new&lat = yo
    eva_new&lon = xo

    eva_new = (/rmMonAnnCycTLL(eva_new)/)

    p_e  = p
    e_p = (/eva_new-p/)
;;;;


    ;;;read era data

    f4 = addfile("/home/yangsong3/data-observation/NCEP_NCAR/monthly/shum.mon.mean.nc","r")
    date :=cd_calendar(f4->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    q = f4->shum(timeind,{1000:300},:,:)
    q = q/1000.0
    ;;density of water is 10^3 kg/m3, intergrated qv units: mm/s


  
    pver=q&level
    pver=pver*100.0
    pver@units="Pa"

    dp=dpres_plevel(pver,100000.0,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(q),dp,1)          ; confirm dp to caculate sum(var*dp)
    q_col = (/dim_sum_n(dp_cfm*q,1)/)               ; sum(var*dp)
    copy_VarCoords(q(:,1,:,:), q_col)

    q_col = (/center_finite_diff_n(q_col, 1,False, 0, 0)/)
;;;;    



; ;;;时间序列
;     evp_annual = month_to_annual(eva, 1)
;     p_annual  = month_to_annual(p, 1)
;     q_annual  = month_to_annual(q_col, 1)
    
;     evp_12 = split_12(eva)
;     p_12 = split_12(p)
;     q_12 = split_12(q_col)


;     evp_annual = (/dim_avg_n_Wrap(evp_12(0:4,:,:,:), 0)/)
;     p_annual = (/dim_avg_n_Wrap(p_12(0:4,:,:,:), 0)/)
;     q_annual = (/dim_avg_n_Wrap(q_12(0:4,:,:,:), 0)/)


;     p_e = p_annual
;     ; p_e = (/p_annual - evp_annual/)
;     copy_VarCoords(p_annual, p_e)
    
;     ;;;
;     lat1 = 0
;     lat2 = 20
;     lon1 = 250
;     lon2 = 280

;     p_e_t = dim_avg_n_Wrap(p_e(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
;     p_annual_t = dim_avg_n_Wrap(p_annual(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
;     evp_annual_t = dim_avg_n_Wrap(evp_annual(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
;     q_annual_t = dim_avg_n_Wrap(q_annual(:,{lat1:lat2},{lon1:lon2}), (/1,2/))
; ;;;

; ;;;corr 
; printVarSummary(q_annual)
; printVarSummary(p_annual)




   ; corr = escorc_n(q_col, p,0, 0)
   corr = escorc_n(q_col, p,0, 0)

   copy_VarCoords(q_col(1,:,:), corr)

   prob_mam = rtest(corr, dimsizes(q_col(:,1,1)), 0)

   copy_VarCoords(q_col(1,:,:), prob_mam)

 ;==================================
wks   = gsn_open_wks("pdf", "q与蒸发-降水的单点相关-JFMAM")
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
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn=True

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  

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
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2
  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
  res@lbLabelStride = 1

  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_blu_red"

  res@gsnScalarContour        = True
 


  ; ; res@vcGlyphStyle="CurlyVector"
  ; res@vcGlyphStyle="LineArrow"
  ; res@vcMonoLineArrowColor = True
  ; res@vcLineArrowColor = "black"
  ; res@vcMinDistanceF           = 0.015          ; thin vectors

  ; res@vcRefLengthF             = 0.05    ; define length of v
  ; res@vcVectorDrawOrder="PostDraw"

  ; ; resv@vcMonoLineArrowColor = True
  ; ; resv@vcLineArrowColor = "blue"
  ; res@vcLineArrowThicknessF=1.2
  ; ; res@vcLineArrowHeadMaxSizeF=0.01
  ; ; res@vcLineArrowHeadMinSizeF=0.01
  ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  ; res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcRefMagnitudeF          = 0.2*10
  ; ; res@vcMinMagnitudeF          =.3

  ; res@vcRefAnnoString1On = True
  ; res@vcRefAnnoString1 = "2.0 m/s"
  ; res@vcRefAnnoOrthogonalPosF  = -1.21              ; move ref vector down -1.0在图的右上
  ; res@vcRefAnnoParallelPosF = 1.02
  ; res@vcRefAnnoFontHeightF = 0.013
  ; res@vcRefAnnoBackgroundColor = "Transparent"
  ; res@vcRefAnnoPerimOn = False


   ; res@gsnLeftString="Corr (NEP 1000-500 integrated Lq & SST UV850)"  
   res@gsnLeftString="Corr (NEP 1000-300 integrated Lq & precipitation)"  
    ; res@gsnLeftString="Corr (NEP 1000-500 integrated Lq & SST UV850, remove Nino3.4 & IND SST)"  
  ; res@gsnLeftString="Corr (NEP 1000-500 integrated Lq & SST UV850, remove Nino3.4 & NTA SST)"  

  res@gsnLeftStringFontHeightF = 0.015
  ; res@gsnLeftString = "The components of JFMAM shum & 850UV trend linked to TNA"
  ; res@gsnLeftStringOrthogonalPosF  = -0.13
 
  plot = new(4,graphic)
  plot(0) = gsn_csm_contour_map(wks, corr, res)
  ; plot(0) = gsn_csm_vector_scalar_map(wks,rc_mam_u,rc_mam_v,rc_mam, res)


  str1 = unique_string("dum1")  ; "unique_string" will return a unique    


  pgres=True
  pgres@gsLineColor="black"

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

  str2 = ispan(21, 24, 1) +""
  do i =0,3 
   plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
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
    iplot = new(4,graphic)
     i=0
    iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)
  
     overlay(plot(i), iplot(i))
  ; end do

    ;;;  
   
; pres = True
; pres@gsnPanelLabelBar  = True
; pres@pmLabelBarOrthogonalPosF = -0.01
draw(plot(0))
frame(wks)
; gsn_panel(wks, plot, (/2,2/), pres)

end


