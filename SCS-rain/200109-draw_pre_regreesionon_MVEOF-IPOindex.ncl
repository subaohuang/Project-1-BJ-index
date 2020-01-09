function  rc_test(ssta_rc,id1,id2)
begin
   dim = dimsizes(id1)
   varname = getVarDimNames(id2)
  
   ;;有效自由度 
   r1_r = esacr(id1, 1)  
   r2_r = esacr(id2(lat|:,lon|:,year|:), 1)

   r1_o = r1_r(1)
   r2 = r2_r(:,:,1)
   f = dimsizes(r2)

   r1  = onedtond(r1_o, dimsizes(r2))

   v = dim*(/1-r1*r2/)/(/1+r1*r2/) 

 
   ssta_rc_new = onedtond(ssta_rc, dimsizes(ssta_rc))
   tval = sqrt(v-2)*(/ssta_rc_new/)/(sqrt(1- ssta_rc_new^2))

   df   = v-2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end

begin
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    p=f2->precip(timeind,:,:)
    p=(/rmMonAnnCycTLL(p)/)
    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)
    p_annual = month_to_annual(p(time|:,lat|:,lon|:), 1)
    p_annual = (/runave_n_Wrap(p_annual, 13, 1, 0)/)
    ;;注意点
     
    ;;;read MV- EOF PC1 data

    PC = asciiread("/home/ys17-19/lsh/Project/SCS-rain/sst_Q_pre_PC2_1958_2013.ascii", 56, "float")
    
    id1 = 1979-1958

    IPO = PC(id1:)
    printVarSummary(IPO)
    
    IPO=(/dim_standardize_n_Wrap(IPO, 0, 0)/)

    sst_index_rc = regline(ispan(1, dimsizes(IPO), 1), IPO)


    IPO = (/dtrend_msg_n(ispan(1,dimsizes(IPO),1),IPO ,False,True,0)/)
    
    ; sst_index = IPO
    sst_index = runave_n_Wrap(IPO, 13, 1, 0)

    dim = dimsizes(IPO)
 


 
    ;;;回归pre
     p_rc = regCoef_n(sst_index, p_annual, 0, 0)
     copy_VarCoords(p_annual(1,:,:), p_rc)


     p_rc = (/p_rc*sst_index_rc/)



    ; ;;
    
    ;;检验

   ;;;显著性检验,注意是用相关系数进行回归系数的显著性检验
    p_r = escorc_n(sst_index, p_annual, 0, 0)

    prob_p = rc_test(p_r,sst_index,p_annual)
    copy_VarCoords(p_annual(1,:,:), prob_p)


    
    ;;;
  wks = gsn_open_wks("pdf", "降水一致性到MVEOF-PC2-IP0指数")

  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnLeftString="(a) OLR"
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
  res@mpMinLonF             = 90
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

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.7

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.08
  res@cnMaxLevelValF       =   0.08
  res@cnLevelSpacingF      =  0.01

  res@cnFillPalette="MPL_BrBG"
  res@gsnLeftString="GPCP Precipitation trend congruent with IPO index"
  plot = gsn_csm_contour_map(wks, p_rc, res)
   
  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3","4"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   i = 0
   ; do i = 0,3
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   ; end do

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
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(1,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_p,res_pro)

    i = 0
    overlay(plot(i), iplot(i))
      
   draw(plot)
   frame(wks)
end



    ;;;




    

