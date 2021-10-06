;;********************************************************
;;画HadIsst的气候态和梯度
;;lnshuheng 2020.02.17
;;;
;**********************************************************
function gradient(data)
begin 
    gradLatLon  = grad_latlon_cfd (data, data&latitude, data&longitude, True, False)
                             
    T_grad_lat  = gradLatLon[0]    ; for clarity; explicitly extract variables from returned 'list' 
    T_grad_lon  = gradLatLon[1]

    T_grad_lat = (/T_grad_lat*10^5/)
    T_grad_lon = (/T_grad_lon*10^5/)
    delete(gradLatLon)
    ; replace the generic attributes with more specific attributes
    T_grad_lon@long_name = " Zonal gradient (derivative)"
    T_grad_lat@long_name = "Meridional gradient (derivative)"
    T_grad_lat@units     = "K/m"
    T_grad_lon@units     = "K/m"

return [/T_grad_lat, T_grad_lon/]
end 

begin
  startmon = 197901
  endmon = 201512
  print("yunxing")
   ;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.197901)

    sst  = lonFlip(f1->sst(timeind,::-1,:))
    printVarSummary(sst)

    sst_mam_clm = dim_avg_n_Wrap(month_to_season(sst, "MAM"), 0)
    sst_jja_clm = dim_avg_n_Wrap(month_to_season(sst, "JJA"), 0)
    sst_son_clm = dim_avg_n_Wrap(month_to_season(sst, "SON"), 0)
    sst_djf_clm = dim_avg_n_Wrap(month_to_season(sst, "DJF"), 0)

    sst_clm = [/sst_mam_clm,sst_jja_clm,sst_son_clm,sst_djf_clm/]
    ;;;;

    ;;;读入ERA SLP数据
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"slp"/)
    filein = path_wind + var_name +".mon.mean.nc"


    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)

    slp = short2flt(f1->msl(timeind,::-1,:))
    slp = (/slp/100.0/) ;; convert to hPa
    
    slp_mam_clm = dim_avg_n_Wrap(month_to_season(slp, "MAM"), 0)
    slp_jja_clm = dim_avg_n_Wrap(month_to_season(slp, "JJA"), 0)
    slp_son_clm = dim_avg_n_Wrap(month_to_season(slp, "SON"), 0)
    slp_djf_clm = dim_avg_n_Wrap(month_to_season(slp, "DJF"), 0)

    slp_clm = [/slp_mam_clm,slp_jja_clm,slp_son_clm,slp_djf_clm/]
        
    print("=====finish read data==========")

;;;;
   ;;;gradient 
     sst_mam_gradient = gradient(sst_mam_clm)
     sst_jja_gradient = gradient(sst_jja_clm)
     sst_son_gradient = gradient(sst_son_clm)
     sst_djf_gradient = gradient(sst_djf_clm)


     sst_clm_zonal_gradient = [/sst_mam_gradient[0],sst_jja_gradient[0],sst_son_gradient[0],sst_djf_gradient[0]/]
     sst_clm_meri_gradient = [/sst_mam_gradient[1],sst_jja_gradient[1],sst_son_gradient[1],sst_djf_gradient[1]/]
   ;;;;
    
   ;;;gradient 
     slp_mam_gradient = gradient(slp_mam_clm)
     slp_jja_gradient = gradient(slp_jja_clm)
     slp_son_gradient = gradient(slp_son_clm)
     slp_djf_gradient = gradient(slp_djf_clm)

     slp_clm_zonal_gradient = [/slp_mam_gradient[0],slp_jja_gradient[0],slp_son_gradient[0],slp_djf_gradient[0]/]
     slp_clm_meri_gradient = [/slp_mam_gradient[1],slp_jja_gradient[1],slp_son_gradient[1],slp_djf_gradient[1]/]
   ;;;;

    ;;;region 
    lat1 = -10
    lat2 = 50
    lon1 = 60
    lon2 = 180
    ;;;

  plot = new(12,graphic)

  wks   = gsn_open_wks("pdf", "global季节平均SLP的气候态以及梯度-philipian-79-16")

  season = (/"MAM","JJA","SON","DJF"/)
  order = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"/)
  
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "BlueDarkRed18"
  ; res@cnFillPalette = "BlRe"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False


  res@lbOrientation="Horizontal"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -30
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             =0
  ; res@mpMaxLonF             = 360

  res@mpMinLatF             = -5
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 100
  res@mpMaxLonF             = 170

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"


  res@tmXBLabelFontHeightF=0.025
  res@tmYLLabelFontHeightF=0.025
  res@tmXBMajorThicknessF=2.0
  res@tmYLMajorThicknessF=2.0
  res@tmBorderThicknessF=2.0

  res@tmXBValues = (/100,120,140,160/) 
  res@tmXBLabels = (/100,120,140,160/) +""+"E" 

  res@tmYLValues = (/0,10,20,30/) 
  res@tmYLLabels = (/"0","10N","20N","30N"/) 

  res@tmYLMinorValues = ispan(-5,30,5)
   ; res@tmYLMinorValues = (/-5,5,15,25/)
  res@gsnLeftStringFontHeightF = 0.03

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   1000.0
  res@cnMaxLevelValF       =   1020.0
  res@cnLevelSpacingF      =  2.0
  res@cnFillPalette = "NCV_jaisnd"


  ; res@cnFillPalette = "sunshine_9lev"
  res@cnInfoLabelOn = False
 
  res@gsnRightString = ""
  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7
 


  do i=0,3 

     res@cnLinesOn=True
     res@gsnLeftString= order(i) + " " + season(i) + " SLP Climatology"
     res@gsnRightString = "[hPa]"
     if i.eq.3 then 
      res@lbLabelBarOn=True
      res@lbOrientation = "Vertical"
      res@pmLabelBarWidthF = 0.2
      res@pmLabelBarHeightF = 0.4
      res@lbLabelStride = 2
     else 
      res@lbLabelBarOn=False
     end if 

     plot(i) = gsn_csm_contour_map(wks, slp_clm[i], res)
  end do 

  res@cnLinesOn =False
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.1
  res@cnFillPalette = "BlueWhiteOrangeRed"

  printMinMax(sst_clm_zonal_gradient[0], 1)


  do i=4,7 
      if i.eq.7 then 
      res@lbLabelBarOn=True
      res@lbOrientation = "Vertical"
      res@pmLabelBarWidthF = 0.2
      res@pmLabelBarHeightF = 0.4
     else 
      res@lbLabelBarOn=False
     end if 
     res@gsnLeftString= order(i) + " " + season(i-4) + " SLP Meridional grad"
     res@gsnRightString = "[10~S~-5~N~ hPa/m]"
     plot(i) = gsn_csm_contour_map(wks, slp_clm_zonal_gradient[i-4], res)
  end do 

  do i=8,11 
     res@gsnLeftString= order(i) + " " + season(i-8) + " SLP Zonal grad"
     res@gsnRightString = "[10~S~-5~N~ hPa/m]"

     if i.eq.11 then 
      res@lbLabelBarOn=True
      res@pmLabelBarWidthF = 0.2
      res@pmLabelBarHeightF = 0.4
     else 
      res@lbLabelBarOn=False
     end if 

     plot(i) = gsn_csm_contour_map(wks, slp_clm_meri_gradient[i-8], res)
  end do 


  ; ;;;添加赤道直线
  ;  pgx=ispan(30,390,10)
  ;  pgy=new(dimsizes(pgx),typeof(pgx))
  ;  pgy=0
   
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF= 1.0
   pgres@gsLineDashPattern= 1


  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF = 1.0
   pgres@gsLineDashPattern= 1
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

    str2 = ispan(1, 12, 1) + ""
    plot_add2 = new(12,graphic)
    do i = 0,11
      plot_add2(i) = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
    end do 

  pres = True
  pres@gsnPanelLabelBar  = False
  pres@gsnMaximize = False
  ; pres@lbOrientation = "Vertical"
  pres@lbOrientation = "Horizontal"
  ; pres@pmLabelBarHeightF = 0.3
  pres@lbLabelFontHeightF = 0.02
  ; pres@lbLabelStride = 2
  pres@gsnPanelRight = 0.95
  ; drawNDCGrid(wks)
  gsn_panel(wks, plot, (/3,4/), pres)


end