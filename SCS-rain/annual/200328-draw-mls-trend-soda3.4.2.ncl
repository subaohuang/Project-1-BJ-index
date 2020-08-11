;;********************************************************
;;计算SODA数据中的平流项和垂直项
;;soda数据的精度是0.5×0.5°
;;NCEP-NCAR的精度为1.87×1.9
;;
;;                                         lnshuheng 2020.02.21
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
 

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/SODA/3.4.2/"
    fils1=systemfunc("ls "+diri+"soda3.4.2_mlt_mn_1980-2015_chazhi.nc")
    f1=addfile(fils1, "r")
    indyear = (2013-1980+1)*12-1

    mlt_r = f1->mlt(0:indyear,{lat1:lat2},{lon1:lon2})
    mlt = tofloat(mlt_r)
    copy_VarCoords(mlt_r, mlt)
    delete(mlt_r)

    mlt&lat@units = "degrees_north"
    mlt&lon@units = "degrees_east"
    printVarSummary(mlt)

    mlt_clm = dim_avg_n_Wrap(mlt, 0)

    mlt_clm_region = mlt_clm({lat1:lat2},{lon1:lon2})



    mlt_annual = month_to_annual(mlt, 1)
 
    dim = 34
    mlt_rc = regCoef_n(ispan(1, 34, 1), mlt_annual, 0, 0)
    copy_VarCoords(mlt_annual(1,:,:), mlt_rc)
    mlt_rc = (/mlt_rc*34/)


   
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "scs-mlt的趋势分布-soda.3.4.2")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlueRed"
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
  res@mpMinLonF             =lon1
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
  res@cnMinLevelValF       =   -30
  res@cnMaxLevelValF       =   30
  res@cnLevelSpacingF      =  10
 

  plot = new(1,graphic)
  
  plot(0) = gsn_csm_contour_map(wks, mlt_rc, res)
  
  
  


  ; ;;;添加赤道直线
  ;   pgx=ispan(30,390,10)
  ;   pgy=new(dimsizes(pgx),typeof(pgx))
  ;   pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
  ;   str = "" + ispan(1,2,1)

  ;  do i = 0,1
  ;    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
  ;  end do


  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   str2 = "0" + ispan(1,2,1)
   i = 0
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   

   draw(plot)
   frame(wks)

end