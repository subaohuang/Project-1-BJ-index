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
  yearstart = 190001
  yearend = 201212

  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360

;;;;read data
   var= (/"air.mon.mean.nc","shum.mon.mean.nc","hgt.mon.mean.nc"/)

    diri = "/home/yangsong3/data-observation/20thCRV2/"

    f1 = addfile(diri+var(0), "r")
    f2 = addfile(diri+var(1), "r")
    f3 = addfile(diri+var(2), "r")

    time1 = f1->time 
    time2 = f2->time 
    time3 = f3->time 

    date1 = cd_calendar(time1, 1)
    date2 = cd_calendar(time2, 1)
    date3 = cd_calendar(time3, 1)
    

    timeind1 :=ind(date1.le.yearend.and.date1.ge.yearstart)
    timeind2 :=ind(date2.le.yearend.and.date2.ge.yearstart)
    timeind3 :=ind(date3.le.yearend.and.date3.ge.yearstart)
   
    t = f1->air(timeind1,{1000:100},{lat1:lat2},{lon1:lon2})

    ; t=t+273.15

    q = f2->shum(timeind2,{1000:100},{lat1:lat2},{lon1:lon2})
     ;;;convert to kg/kg

    hgt = f3->hgt(timeind3,{1000:100},{lat1:lat2},{lon1:lon2})
;;;

 
 

Cp = 1004
Lv = 2.5*10^6
g = 1.0

printVarSummary(t)
printVarSummary(q)
printVarSummary(hgt)
MSE = t
MSE = (/Cp*t + Lv*q + g*hgt/)


 
 ; system("/bin/rm MSE_era5.nc")
 out = addfile("MSE_20CR_1900-2012.nc", "c")
 out->MSE =  MSE 
 exit()
 ; f1 = addfile("MSE_era5.nc", "r")
 ; MSE = f1->MSE



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

 MSE_nep_annual = month_to_annual(MSE_nep, 1)
 MSE_nep_annual!0 = "time"

 MSE_nep_MAM = month_to_season(MSE_nep, "MAM")
 MSE_nep_JJA = month_to_season(MSE_nep, "JJA")
 MSE_nep_SON = month_to_season(MSE_nep, "SON")
 MSE_nep_DJF = month_to_season(MSE_nep, "DJF")

 
 ;;;计算趋势
  MSE_nep_annual&time = ispan(1950,2019,1)
  MSE_nep_MAM&time = ispan(1950,2019,1)
  MSE_nep_JJA&time = ispan(1950,2019,1)
  MSE_nep_SON&time = ispan(1950,2019,1)
  MSE_nep_DJF&time = ispan(1950,2019,1)


  ;;;trend year 
  year1 = 1991
  year2 = 2019

  rc_annual  = mse_trend(MSE_nep_annual, year1,year2)
  rc_mam  = mse_trend(MSE_nep_MAM, year1,year2)
  rc_jja  = mse_trend(MSE_nep_JJA, year1,year2)
  rc_son  = mse_trend(MSE_nep_SON, year1,year2)
  rc_djf  = mse_trend(MSE_nep_DJF, year1,year2)

  
  plot = new(4,graphic)
  prob_mam = rc_test(rc_mam)
  prob_jja = rc_test(rc_jja)
  prob_son = rc_test(rc_son)
  prob_djf = rc_test(rc_djf)

  printVarSummary(prob_mam)
;;;;

    
wks   = gsn_open_wks("pdf", "MSE的趋势-"+year1+"-"+year2)
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
  res@lbLabelBarOn=False
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


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.5



  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "BlueDarkRed18"

  ; res@gsnLeftString="Annual mean MSE trend " 
  plot = new(4,graphic)
  res@gsnLeftString="MAM MSE trend " 
  plot(0) = gsn_csm_contour_map(wks, rc_mam, res)
  res@gsnLeftString="(b)JJA MSE trend " 
  plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  res@gsnLeftString="(c)SON MSE trend " 
  plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  res@gsnLeftString="(d)DJF MSE trend " 
  plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 


   pgres=True
   pgres@gsLineColor="black"

   pgres@gsLineThicknessF=1.0
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
   do i =0,3 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 


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

    iplot=new(4,graphic)

    
 
    iplot(0) = gsn_csm_contour(wks,prob_mam,res_pro)    
    iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)    
    iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)    
    iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)    

    ; iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)
    ; iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)
    ; iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)

    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

   
     do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
     overlay(plot(i), iplot(i))
  end do

    ;;;  
   
pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01
; draw(plot(0))
; frame(wks)
gsn_panel(wks, plot, (/2,2/), pres)

end


