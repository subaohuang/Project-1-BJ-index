;;********************************************************
;;画海温的趋势，注意有表层海洋以下的海温趋势.
;;
;;soda数据已经插值到1°乘1°
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
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -5
   lat2 = 45
   lon1 = 90
   lon2 = 160
 


;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    

    f5=addfile("/home/yangsong3/data-observation/SODA/3.3.1/soda3.3.1_temp_mn_1980-2015_chazhi.nc", "r")


    indyear = (2013-1980+1)*12-1
    T = f5->temp(0:indyear,{0:2000},{lat1:lat2},{lon1:lon2})

    T_mean = dim_avg_n_Wrap(T, 2)
    T_annual = month_to_annual(T_mean, 1)


  
print("*******计算趋势*****************")
;;;calculate correlation coef

    ;;;hadv correlation with sst
    ; r_T =   regCoef_n(ispan(1,yearnum,1), T_annual, 0, 0) 
    r_T = dim_avg_n_Wrap(T_annual, 0)
    copy_VarCoords(T_annual(0,:,:), r_T)



    printVarSummary(r_T)
    
;;;;

print("*******显著性检验*****************")

;;;t test 
  ; prob_p1= rc_test(r_T)
  ; copy_VarCoords(r_T, prob_p1)




print("*********开始画图********")

    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "海洋的温度的趋势")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@trYReverse  = True
 

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

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -0.1
  ; res@cnMaxLevelValF       =   0.1
  ; res@cnLevelSpacingF      =  0.01
 

   plot=gsn_csm_contour(wks, r_T, res)
   draw(plot)
   frame(wks)
end