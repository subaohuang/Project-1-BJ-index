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

; ;;;;read data
;     var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc"/)

;     diri = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"

;     f1 = addfile(diri+var(0), "r")
;     f2 = addfile(diri+var(1), "r")
;     f3 = addfile(diri+var(2), "r")

;     time = f1->time 
;     date = cd_calendar(time, 1)
;     delete(time)

;     timeind :=ind(date.le.yearend.and.date.ge.yearstart)
;     print(date(ind(date.eq.yearstart)))
;     delete(date)     
;     t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

;     t=t+273.15

;     q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

;     q= q/1000.0

;     hgt = f3->hgt(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
; ;;;

; Cp = 1004
; Lv = 2.5*10^6
; g = 9.8

;  MSE = t
;  MSE = (/Cp*t + Lv*q + g*hgt/)
 
;  systemfunc("/bin/rm MSE_necp1.nc")
;  out = addfile("MSE_necp1.nc", "c")
;  out->MSE =  MSE 
;  exit()
    data = "ncep1"

    ; f1 = addfile("MSE_era5.nc", "r")
    f1 = addfile("/home/ys17-19/lsh/Project/MSE-EP/MSE_"+data+"_50-20.nc", "r")

    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    startmon = 198001
    endmon = 202012

    startyear = 1980
    endyear = 2020

    timeind :=ind(date.le.endmon.and.date.ge.startmon)

    MSE = f1->MSE(timeind,:,:,:)


;;;MSE 
    plev = MSE&level
    ptop = min(plev)
    psfc = 1013.0
    dz = dpres_plevel(plev, psfc, ptop, 0)
    MSE_col =  dim_sum_wgt_n_Wrap(MSE, dz, 1, 1)
    MSE_col = (/rmMonAnnCycTLL(MSE_col)/)

    ; MSE_nep = dim_avg_n_Wrap(MSE_col, (/1,2/))
    MSE_nep = MSE_col
    MSE_nep = (/MSE_col/10^6/)
;;;;
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p=f2->precip(timeind,:,:)

    p = (/rmMonAnnCycTLL(p)/)
;;;;

    f3=addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/evapr_oaflux_chazhi.nc", "r")
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    eva=f3->evapr(timeind,:,:)
    eva = (/eva*10/365/) ;;convert to m/day
    eva = (/rmMonAnnCycTLL(eva)/)
;;;;


; ;;;e-p
;     e_p = p
;     e_p = p-eva
;     p = (/p-eva/)
; ;;;


;;;;;;split to 12 month  is a array size (12,yearnum)
  
     MSE_nep_12 = split_12(MSE_nep)
     p_12 = split_12(p)
   
;;;;;


;;;等到年平均的数据
    MSE_nep_annual = month_to_annual(MSE_nep, 1)
    MSE_nep_annual!0 = "time"
    MSE_nep_annual&time = ispan(startyear,endyear,1)

    p_annual = month_to_annual(p, 1)
    p_annual!0 = "time"
    p_annual&time = ispan(startyear,endyear,1)
    
    p_annual = (/dim_avg_n(p_12(0:4,:,:,:),0)/)
    MSE_nep_annual = (/dim_avg_n(MSE_nep_12(0:4,:,:,:),0)/)
;;;;;



 ;;;计算趋势
  ;;trend year 
  year1 = startyear
  year2 = endyear

  rc_mse  = mse_trend(MSE_nep_annual, year1,year2)
  rc_p  = mse_trend(p_annual, year1,year2)
 
  
  prob_mam = rc_test(rc_mse)
  ; prob_mam = rc_test(rc_p)
  printVarSummary(prob_mam)
;;;


; wks   = gsn_open_wks("pdf", ""+data+"-MSE以及降水的趋势-"+year1+"-"+year2)
  wks   = gsn_open_wks("pdf", ""+data+"-MSE以及蒸发-降水的趋势-"+year1+"-"+year2)
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
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -2.0
  ; res@cnMaxLevelValF       =   2.0
  ; res@cnLevelSpacingF      =  0.5

  res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour level
  res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)
  res@lbLabelFontHeightF = 0.01
  res@lbLabelStride = 1
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueDarkRed18"
  res@cnFillPalette = "BlueWhiteOrangeRed"

  res@mpFillOn = True
  res@mpFillColor = "White"
  res@mpFillDrawOrder ="PostDraw"

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpWidthF = 0.6
  ; res@vpHeightF = 0.1

    ; res@gsnLeftString="Annual mean MSE trend " 
  plot = new(1,graphic)

  ; res@tiMainString = data+ " MSE trend " + startyear + "-" + endyear
  res@gsnLeftString="JFMAM MSE&PRE trend " 
  plot = gsn_csm_contour_map(wks, rc_mse, res)
  ; plot = gsn_csm_contour_map(wks, rc_p, res)
  
  printMinMax(rc_p, 1)
 

  res2 = res
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnFillOn = False
  res2@cnLinesOn = True
  ; res2@cnLevels := ispan(-4, 4, 1)
  ; res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour level
  res2@cnLevels = (/-4.0,-2.0,-1.0,1.0,2.0,4.0/)

  res2@cnLineLabelStrings =   res2@cnLevels  + ""
  res2@cnLineLabelDensityF = 1.3

  ; res2@cnLineColor = "green"
  res2@gsnContourNegLineDashPattern = 1
  ; res2@cnLineColors = (/"orange4","orange4","orange4","palegreen4","palegreen4","palegreen4"/)

  res2@cnLineColors = (/"palegreen4","palegreen4","palegreen4","blue4","blue4","blue4"/)
   
  res2@cnMonoLineColor  = False
  res2@cnLineThicknessF = 1.5
  res2@cnMonoLineThickness  = True
  res2@gsnLeftString = ""
  res2@cnInfoLabelOn = False
  ; res2@cnLineLabelsOn = True
  res2@cnLineLabelInterval = 1
  ; res2@cnLineLabelBackgroundColor = "White"

  plot2 = gsn_csm_contour(wks, rc_p, res2)

  overlay(plot, plot2)

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

  plot@add = gsn_add_polyline(wks, plot, pgx_pool, pgy_pool, pgres)

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
 
  iplot = gsn_csm_contour(wks,prob_mam,res_pro)

  overlay(plot, iplot)

    ;;  
   
; pres = True
; pres@gsnPanelLabelBar  = True
; pres@pmLabelBarOrthogonalPosF = -0.01
draw(plot)
frame(wks)
; gsn_panel(wks, plot, (/2,2/), pres)

end


