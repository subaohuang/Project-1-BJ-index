;;********************************************************
;;draw Local Walker and Hadley circulations trend 
;;
;;                                         lnshuheng 2020.02.18
;**********************************************************

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
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 202012

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -60
   lat2 = 60
   lon1 = 0
   lon2 = 360
 

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/obs-derived/msf/"
    fils1=systemfunc("ls "+diri+"m_zonal_meri_ERA5_500hpa.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.202012.and.date.ge.197901)

    wk = f1->m_zonal(timeind,:,:)
    printVarSummary(wk)
    printMinMax(wk, 1)

    wk_a = (/rmMonAnnCycTLL(wk)/)
    copy_VarCoords(wk, wk_a)
    
    hc = f1->m_meri(timeind,:,:)
    printVarSummary(hc)
    printMinMax(hc, 1)

    hc_a = (/rmMonAnnCycTLL(hc)/)
    copy_VarCoords(hc, hc_a)

    wk_annual = month_to_annual(wk_a, 1)
    hc_annual = month_to_annual(hc_a, 1)
    

  
print("*******计算趋势*****************")
;;;calculate correlation coef
    rc_wk =   regCoef_n(ispan(1, yearnum, 1), wk_annual, 0, 0)
    copy_VarCoords(wk_annual(0,:,:), rc_wk)

    rc_hc =   regCoef_n(ispan(1, yearnum, 1), hc_annual, 0, 0)
    copy_VarCoords(hc_annual(0,:,:), rc_hc)

    rc_wk = (/rc_wk*yearnum/)
    rc_hc = (/rc_hc*yearnum/)
;;;


    
;;;;

print("*******显著性检验*****************")

;;;t test 
  prob_p1= rc_test(rc_wk)
  prob_p2= rc_test(rc_hc)
  
  ; rc_wk = (/dim_avg_n(wk, 0)/)
  ; rc_hc = (/dim_avg_n(hc, 0)/)
  
 
print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "ERA5-local-WK&HC-trend-"+startmon+"-"+endmon)

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlueWhiteOrangeRed"
  ; res@cnFillPalette="BlueRed"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2


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

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.002
  res@cnMaxLevelValF       =   0.002
  res@cnLevelSpacingF      =  0.0005
   

  plot = new(2,graphic)
  res@gsnLeftString = "Obs Wk trend"
  plot(0) = gsn_csm_contour_map(wks, rc_wk, res)        
  res@gsnLeftString = "Obs HC trend"
  plot(1) = gsn_csm_contour_map(wks, rc_hc, res)        


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,4,1)

     do i = 0,1
      plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
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
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(2,graphic)

  
    iplot(0) = gsn_csm_contour(wks,prob_p1,res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_p2,res_pro)
   

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


     do i = 0,1
      iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
        overlay(plot(i), iplot(i))
     end do 

    ;;;  
    
 pres = True
 pres@gsnPanelLabelBar = True
    
  gsn_panel(wks, plot,(/2,1/), pres)

end