;;********************************************************
;;计算海洋动力过程项,SST倾向以及SST值.
;;的相关系数场，
;;NCEP-NCAR的精度为1.87×1.9,soda数据已经插值到1°乘1°
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
    diri="/home/yangsong3/data-observation/SST/HadISST/"
    fils1=systemfunc("ls "+diri+"HadISST_sst_1870_2020.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.202012.and.date.ge.197901)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)


    ssta_annual = month_to_annual(ssta(:,{lat1:lat2},{lon1:lon2}), 1)
    printVarSummary(ssta_annual)

  
print("*******计算趋势*****************")
;;;calculate correlation coef
    
     ;;;u correlation with sst
    r_sst_sst =   regCoef_n(ispan(1, yearnum, 1), ssta_annual, 0, 0)
 
    copy_VarCoords(ssta_annual(0,:,:), r_sst_sst)

    r_sst_sst = (/r_sst_sst*yearnum/)

;;;
    
;;;;

print("*******显著性检验*****************")

;;;t test 
  prob_p1= rc_test(r_sst_sst)
  copy_VarCoords(r_sst_sst, prob_p1)

  
 

print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "Hadisst的趋势"+startmon+"-"+endmon)

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlueWhiteOrangeRed"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
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
  res@cnMinLevelValF       =   -1
  res@cnMaxLevelValF       =   1
  res@cnLevelSpacingF      =  0.2
   

  plot = new(1,graphic)
  res@gsnLeftString = "SST趋势"
  plot(0) = gsn_csm_contour_map(wks, r_sst_sst, res)        


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,4,1)

     i = 0
     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)


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

    iplot=new(1,graphic)

  
    iplot(0) = gsn_csm_contour(wks,prob_p1,res_pro)
   

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數


     i = 0
     iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    

    ;;;  
    
      overlay(plot(i), iplot(i))
    

  draw(plot)
  frame(wks)


end