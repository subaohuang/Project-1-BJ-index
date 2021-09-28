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
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data,0,0)
  copy_VarCoords(data(1,:,:), rc1)
  rc1 = rc1*(year2-year1+1)
  return rc1
end 

function mse_diff(data,year1,year2,year3,year4)
begin
  ; rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data({year1:year2},:,:),0,0)
  ; copy_VarCoords(data(1,:,:), rc1)
  ; rc1 = rc1*(year2-year1+1)
  rc2 = dim_avg_n_Wrap(data({year1:year2},:,:),0)
  rc2 = dim_avg_n_Wrap(data({year3:year4},:,:),0) - dim_avg_n_Wrap(data({year1:year2},:,:),0)
  return rc2
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

   lat1 = -50
   lat2 = 50
   lon1 = 0
   lon2 = 360

   startmon = 198001
   endmon = 202012
   ;;read  eva data 
    year = ispan(1980,2020,1) 
   ;;;;
    f2=addfile("/home/yangsong3/data-observation/ERA5-monthly/surface/total_precipitation_mon_1x1_1979_2020.nc", "r")
  
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p = short2flt(f2->tp(timeind,0,{lat1:lat2},{lon1:lon2}))
    p = (/rmMonAnnCycTLL(p)/)

    p = 1000*p  ;;; convert to mm/day

;;;;

    f3=addfile("/home/yangsong3/data-observation/ERA5-monthly/surface/evaporation_mon_1x1_1979_2020.nc", "r")
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    eva= short2flt(f3->e(timeind,{lat1:lat2},{lon1:lon2}))
    eva = (/eva*1000/) ;;convert to mm/day
    eva = (/rmMonAnnCycTLL(eva)/)
;;;;
   printVarSummary(eva)
 

;;;时间序列
    evp_annual = month_to_annual(eva, 1)
    p_annual  = month_to_annual(p, 1)
    
    evp_12 = split_12(eva)
    p_12 = split_12(p)


    eva_annual = (/dim_avg_n_Wrap(evp_12(0:4,:,:,:), 0)/)
    p_annual = (/dim_avg_n_Wrap(p_12(0:4,:,:,:), 0)/)
    


    e_p_annual = p_annual
    e_p_annual = (/eva_annual - p_annual/)
    copy_VarCoords(p_annual, e_p_annual)
    
 ;;;计算趋势
  ;;trend year 
  year1 = 1980
  year2 = 2020

  rc_e_p  = mse_trend(e_p_annual, year1, year2)
  rc_p  = mse_trend(p_annual, year1,year2)
  rc_eva  = mse_trend(eva_annual, year1,year2)

  copy_VarCoords(p(1,:,:), rc_e_p)
  copy_VarCoords(rc_e_p,rc_p)
  copy_VarCoords(rc_e_p,rc_eva)
  
  prob_mam = rc_test(rc_e_p)
  prob_mam_p = rc_test(rc_p)
  prob_mam_eva = rc_test(rc_eva)
  
  printVarSummary(prob_mam)
  printVarSummary(prob_mam_p)
  printVarSummary(prob_mam_eva)
;;;



; wks   = gsn_open_wks("pdf", ""+data+"-MSE以及降水的趋势-"+year1+"-"+year2)
  wks   = gsn_open_wks("pdf", "ERA5-蒸发-降水的趋势-"+year1+"-"+year2)
; wks   = gsn_open_wks("pdf", ""+data+"-降水的趋势-"+year1+"-"+year2)
; 
; wks   = gsn_open_wks("pdf", "MSE的change-"+year1+"-"+year2)
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=0.015
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.2

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn= True
  ; res@lbLabelBarOn=True

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
  res@lbAutoManage = False
  res@lbLabelFontHeightF = 0.3



  res@tmXBValues = (/0,60,120,180,240,300,360/) 
  res@tmXBLabels = (/"0","60E","120E","180","120W","60W","0"/) 

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015


  ; res@tmYLValues = ispan(-20, 20, 20) 
  ; res@tmYLLabels = (/"20S","0","20N"/) 

  ; res@tmYLMinorValues = ispan(-20,20,5)
  

  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.5

  ; res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour level
  ; res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)
  res@lbLabelFontHeightF = 0.01
  res@lbLabelStride = 1
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "BlueDarkRed18"
  ; res@cnFillPalette = "BlueWhiteOrangeRed"

  res@mpFillOn = True
  res@mpFillColor = "White"
  res@mpFillDrawOrder ="PostDraw"

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpWidthF = 0.6
  ; res@vpHeightF = 0.1

    ; res@gsnLeftString="Annual mean MSE trend " 
  plot = new(3,graphic)

  ; res@tiMainString = data+ " MSE trend " + startyear + "-" + endyear
  res@gsnLeftString="JFMAM Pre trend " 
  plot(0) = gsn_csm_contour_map(wks, rc_p, res)

   ; res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)*0.5
     ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.25
  res@gsnLeftString="JFMAM Eva trend " 
  plot(1) = gsn_csm_contour_map(wks, rc_eva, res)

   ; res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)
    res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.5
  res@gsnLeftString="JFMAM E-P trend " 
  plot(2) = gsn_csm_contour_map(wks, rc_e_p, res)
  
  
  printMinMax(rc_p, 1)
 


 

  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.0
  pgres@gsLineDashPattern=0

  ;;;画box
  lat1 = 0
  lat2 = 20
  lon1 = 250
  lon2 = 280

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;
  add = ispan(1, 3, 1)+""  
  do i=0,2
  plot@$add(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do
  ;;;打点
  ; ; end do
  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnFillOn = True

  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLevels = (/0,0.1/)
  res_pro@lbLabelBarOn = False
  res_pro@cnMonoFillColor = False
  res_pro@cnFillColors = (/"Transparent","White","Transparent"/)


  res_pro@cnLinesOn=False
  ; res_pro@cnMonoLineColor =False
  ; res_pro@cnLineColors = ((/"Transparent","black","Transparent"/))

  res_pro@cnMonoFillPattern = True
  res_pro@cnFillPattern = 17
  res_pro@cnFillDotSizeF = 0.002
  
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True
  iplot = plot
  iplot(0) = gsn_csm_contour(wks,prob_mam_p,res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_mam_eva,res_pro)
     iplot(2) = gsn_csm_contour(wks,prob_mam,res_pro)
   do i=0,2
  overlay(plot(i), iplot(i))
  end do
    ;;  
   
pres = True
pres@gsnPanelLabelBar  = False
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,1/), pres)

end


