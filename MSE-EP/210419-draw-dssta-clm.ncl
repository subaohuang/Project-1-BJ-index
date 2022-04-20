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


function cal_clm(data,season)
  begin
   
    p_djf_r=month_to_season(data, season)
    p_rc = dim_avg_n_Wrap(p_djf_r , 0)
  return p_rc
end 


;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201012.and.date.ge.197901)

    sst  = f1->sst(timeind,:,:)
    ; sst  = short2flt(f1->sst(timeind,:,:))

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)


    dsst = center_finite_diff_n(ssta, 1, False, 1, 0)
    
    ssta = (/dsst/)
 
   ;;;;cal season trend 
    ssta_mam = cal_clm(ssta,"MAM")
    ssta_jja = cal_clm(ssta,"JJA")
    ssta_son  = cal_clm(ssta,"SON")
    ssta_djf = cal_clm(ssta,"DJF")

    printMinMax(ssta_mam, 1)
    exit()
    
    ; prob_mam = rc_test(ssta_mam)
    ; prob_jja = rc_test(ssta_jja)

    ; prob_son = rc_test(ssta_son)
    ; prob_djf = rc_test(ssta_djf)

    ; wks   = gsn_open_wks("pdf", "global季节平均SST的趋势-philipian-79-19")
    wks   = gsn_open_wks("pdf", "ERSST-global季节平均dSST的气候态")
    ; wks   = gsn_open_wks("pdf", "ERSST-global季节平均SST的趋势-99-19")
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  
 
   ; res@cnFillPalette = "BlueWhiteOrangeRed"
   ; res@cnFillPalette = "BlueDarkRed18"
   res@cnFillPalette = "BlueWhiteOrangeRed"
   ; res@cnFillPalette = "BlRe"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
   res@mpMinLatF             = -20
  res@mpMaxLatF             = 20
  res@mpMinLonF             = 150
  res@mpMaxLonF             = 290
  

  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             = 100
  ; res@mpMaxLonF             = 170


  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  ; res@tmXBLabelFontHeightF=0.015
  ; res@tmYLLabelFontHeightF=0.015

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -0.12
  ; res@cnMaxLevelValF       =   0.12
  ; res@cnLevelSpacingF      =  0.03

   res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  res@cnInfoLabelOn = False
  

  res@gsnLeftString="MAM"
  res@gsnLeftStringFontHeightF = 0.02

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

  plot=new(4,graphic)

  plot(0) = gsn_csm_contour_map(wks, ssta_mam, res)
    res@gsnLeftString="JJA"
  plot(1) = gsn_csm_contour_map(wks, ssta_jja, res)
    res@gsnLeftString="SON"
  plot(2) = gsn_csm_contour_map(wks, ssta_son, res)
    res@gsnLeftString="DJF"
  plot(3) = gsn_csm_contour_map(wks, ssta_djf, res)


   str1 = unique_string("dum1")  ; "unique_string" will return a unique    


    ; ;;;添加赤道直线
    ; pgx=ispan(30,390,10)
    ; pgy=new(dimsizes(pgx),typeof(pgx))
    ; pgy=0

    pgres=True
    pgres@gsLineColor="red"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1

    ; str = ispan(1, 4, 1) + ""
    ; plot_add = new(4,graphic)
    ; do i = 0,3
    ; plot_add(i) = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
    ; end do 
  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   ; lat1 = 0
   ; lat2 = 20
   ; lon1 = 210
   ; lon2 = 280

   lat1 = 0
   lat2 = -10
   lon1 = 270
   lon2 = 280


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
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.001
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  
    ;  iplot = new(4,graphic)
    ;  iplot(0) = gsn_csm_contour(wks,prob_mam  ,res_pro)
    ;  iplot(1) = gsn_csm_contour(wks,prob_jja  ,res_pro)
    ;  iplot(2) = gsn_csm_contour(wks,prob_son  ,res_pro)
    ;  iplot(3) = gsn_csm_contour(wks,prob_djf  ,res_pro)
    ; ;  
    ; do i=0,3
    ;  overlay(plot(i), iplot(i))
    ; end do

pres = True
pres@gsnPanelLabelBar  = True


pres@pmLabelBarOrthogonalPosF = -0.008

pres@lbLabelStride=1


gsn_panel(wks, plot, (/2,2/), pres)

end