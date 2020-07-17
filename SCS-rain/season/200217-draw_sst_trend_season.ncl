;;********************************************************
;;画HadIsst的sst趋势
;;
;;lnshuheng 2020.02.17
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

function cal_season_trend(data,season)
  begin
    p_djf_r=month_to_season(data, season)
    dim = dimsizes(p_djf_r)
    p_rc=regCoef_n(ispan(1, dim(0), 1), p_djf_r, 0, 0)
    copy_VarCoords(p_djf_r(1,:,:), p_rc)
    p_rc=(/p_rc*dim(0)/)
  return p_rc
end 



;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.197901)

    sst=f1->sst(timeind,:,:)
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
  

   ;;;;cal season trend 
    ssta_mam = cal_season_trend(ssta,"MAM")
    ssta_jja = cal_season_trend(ssta,"JJA")
    ssta_son  = cal_season_trend(ssta,"SON")
    ssta_djf = cal_season_trend(ssta,"DJF")



    prob_mam = rc_test(ssta_mam)
    prob_jja = rc_test(ssta_jja)

    prob_son = rc_test(ssta_son)
    prob_djf = rc_test(ssta_djf)

    ssta_mam = where(prob_mam.le.0.05, ssta_mam, ssta_mam@_FillValue)
    ssta_jja = where(prob_jja.le.0.05, ssta_jja, ssta_jja@_FillValue)
    ssta_son = where(prob_son.le.0.05, ssta_son, ssta_son@_FillValue)
    ssta_djf = where(prob_djf.le.0.05, ssta_djf, ssta_djf@_FillValue)
    
    
    wks   = gsn_open_wks("pdf", "季节平均SST的趋势-SCS-79-15-global")
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  
  ; res@cnFillPalette = "NCV_jaisnd"
   res@cnFillPalette = "BlueWhiteOrangeRed"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -60
  res@mpMaxLatF             = 60
  res@mpMinLonF             =0
  res@mpMaxLonF             = 360

  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             =110
  ; res@mpMaxLonF             = 155



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
  ; res@cnMinLevelValF       =   0.2
  ; res@cnMaxLevelValF       =   0.6
  ; res@cnLevelSpacingF      =  0.05
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.1
  res@cnInfoLabelOn = False
  

  res@gsnLeftString="Spring"

  plot=new(4,graphic)

  plot(0) = gsn_csm_contour_map(wks, ssta_mam, res)
    res@gsnLeftString="Summer"
  plot(1) = gsn_csm_contour_map(wks, ssta_jja, res)
    res@gsnLeftString="Autumn"
  plot(2) = gsn_csm_contour_map(wks, ssta_son, res)
    res@gsnLeftString="Winter"
  plot(3) = gsn_csm_contour_map(wks, ssta_djf, res)

  res2 = res
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnFillOn = False
  res2@cnLinesOn = True
  res2@cnLevels = (/0.55/)
  res2@cnLineColors = "black"
  res2@cnMonoLineColor  = True
  res2@cnLineThicknessF = 5
  res2@cnMonoLineThickness  = True
  res2@gsnLeftString = ""
  res2@cnInfoLabelOn = False
  
  plot1=new(4,graphic)

  plot1(0) = gsn_csm_contour(wks, ssta_mam, res2)
  plot1(1) = gsn_csm_contour(wks, ssta_jja, res2)
  plot1(2) = gsn_csm_contour(wks, ssta_son, res2)
  plot1(3) = gsn_csm_contour(wks, ssta_djf, res2)

  do i=0,3
  overlay(plot(i), plot1(i))
  end do   

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    


  ;;;添加赤道直线
   pgx=ispan(30,390,10)
   pgy=new(dimsizes(pgx),typeof(pgx))
   pgy=0
   
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1

  str = ispan(1, 4, 1) + ""
    plot_add = new(4,graphic)
  do i = 0,3
  plot_add(i) = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   end do 
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

    str2 = ispan(1, 4, 1) + ""
     plot_add2 = new(4,graphic)
    do i = 0,3
      plot_add2(i) = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
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
    res_pro@cnFillDotSizeF = 0.003
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  
  
     ; iplot = gsn_csm_contour(wks,prob_p  ,res_pro)

    ;;;  
    ; overlay(plot, iplot)
pres = True
pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,2/), pres)

end