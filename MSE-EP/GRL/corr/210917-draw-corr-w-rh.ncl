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
begin
  yearstart = 198001
  yearend = 202012

  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360

  f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_necp1.nc", "r")
  MSE = f1->MSE

; ;;;;read data
    var= (/"omega.mon.mean.nc","shum.mon.mean.nc","rhum.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)
    print(date(ind(date.eq.yearstart)))
    delete(date)     
    w = f1->omega(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    
    w = (/-50*w/)

    q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    q= q/1000.0

    rh = f3->rhum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})


;;;
;;;vertical mean  
  plev = q&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)
   
  level = 850
  q_col =  dim_avg_wgt_n_Wrap(q, dz, 1, 1)
   ; q_col = q(:,{level},:,:)
  
  q_col = (/rmMonAnnCycTLL(q_col)/)
  q_col = (/dtrend_msg_n(ispan(1,dimsizes(q_col&time),1),q_col,False,True,0)/)
  

  w_col =  dim_avg_wgt_n_Wrap(w, dz, 1, 1)
    ; w_col = w(:,{level},:,:)
  w_col = (/rmMonAnnCycTLL(w_col)/)
  w_col = (/dtrend_msg_n(ispan(1,dimsizes(w_col&time),1),w_col,False,True,0)/)

  rh_col =  dim_avg_wgt_n_Wrap(rh, dz, 1, 1)
  ; rh_col = rh(:,{level},:,:)
  rh_col = (/rmMonAnnCycTLL(rh_col)/)
  rh_col = (/dtrend_msg_n(ispan(1,dimsizes(rh_col&time),1),rh_col,False,True,0)/)
;;;;


;;;;annual mean 
  q_col_12 = split_12(q_col)
  q_col_annual = month_to_annual(q_col, 1)
  q_col_annual = (/dim_avg_n_Wrap(q_col_12(0:4,:,:,:), 0)/)

  rh_col_12 = split_12(rh_col)
  rh_col_annual = month_to_annual(rh_col, 1)
  rh_col_annual = (/dim_avg_n_Wrap(rh_col_12(0:4,:,:,:), 0)/)

  w_col_12 = split_12(w_col)
  w_col_annual = month_to_annual(w_col, 1)
  w_col_annual = (/dim_avg_n_Wrap(w_col_12(0:4,:,:,:), 0)/)


   ; corr = escorc_n(q_col, p,0, 0)
   corr = escorc_n(rh_col_annual, w_col_annual,0, 0)
   ; corr = escorc_n(q_col_annual, w_col_annual,0, 0)

   copy_VarCoords(q_col_annual(1,:,:), corr)

   prob_mam = rtest(corr, dimsizes(q_col(:,1,1)), 0)

   copy_VarCoords(q_col_annual(1,:,:), prob_mam)

 ;==================================
wks   = gsn_open_wks("pdf", "整层积分相对湿度与垂直速度的单点相关")
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

  res@gsnLeftString="Corr (vertically mean omega & RH)"  

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


