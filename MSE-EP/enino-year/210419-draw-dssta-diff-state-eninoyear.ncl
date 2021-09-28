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


function cal_composite(data,season,ninoyear)
  begin
  if season.eq."DJF" then 
    dssta_mam = month_to_season(data, "DJF")
    dssta_mam(0:dimsizes(dssta_mam&time)-2,:,:) = dssta_mam(1:,:,:)
    dssta_mam(dimsizes(dssta_mam&time)-1,:,:) = data(dimsizes(data&time)-1,:,:)
  else 
    dssta_mam = month_to_season(data, season)
  end if
    dssta_mam_enino = dim_avg_n_Wrap(dssta_mam(ninoyear,:,:), 0)
 
  return dssta_mam_enino
end 


;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;读入HadiSST海温数据
    ; diri="/home/yangsong3/data-observation/linshh/data/sst/"
    ; fils1=systemfunc("ls "+diri+"HadISST_sst.nc")

    diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
     fils1=systemfunc("ls "+diri+"sst.mnmean.v5.nc")

    ; diri = "/home/yangsong3/data-observation/SST/NOAA_OISST/"
    ;  fils1=systemfunc("ls "+diri+"sst.mnmean.nc")
  ; ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1

    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind1=ind(date.le.201912.and.date.ge.197901)
    timeind2=ind(date.le.201012.and.date.ge.197901)

    sst  = f1->sst(timeind1,:,:)
    ; sst  = short2flt(f1->sst(timeind,:,:))
    sst_clm = clmMonTLL(f1->sst(timeind2,:,:))
    
    ssta = calcMonAnomTLL(sst, sst_clm)
    ; ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)


    dsst = center_finite_diff_n(ssta, 1, False, 1, 0)
    copy_VarCoords(ssta, dsst)
    

    ; ely=(/1952,1954,1958,1959,1964,1966,1969,1970,1973,1977,1978,1980,1983,1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016/) - 1
    ; ely1=(/1992,1995,1998/) - 1
    ; ely2=(/2003,2005,2007,2010,2015,2016/) - 1    
    

    ely1=(/1980,1983,1987,1988,1992,1995,1998/)-1
    ely2=(/2003,2005,2007,2010,2015,2016/) - 1    
    
    ind_obs_nino1 = ely1 - 1979
    ind_obs_nino2 = ely2 - 1979

    ; ind_obs_nino1 = ely1 - 1991
    ; ind_obs_nino2 = ely2 - 1991

      
   ;;;;cal season trend 
    dssta_mam1 = cal_composite(dsst,"MAM",ind_obs_nino1)
    dssta_jja1 = cal_composite(dsst,"JJA",ind_obs_nino1)
    dssta_son1  = cal_composite(dsst,"SON",ind_obs_nino1)
    dssta_djf1 = cal_composite(dsst,"DJF",ind_obs_nino1)

    dssta_mam2 = cal_composite(dsst,"MAM",ind_obs_nino2)
    dssta_jja2 = cal_composite(dsst,"JJA",ind_obs_nino2)
    dssta_son2  = cal_composite(dsst,"SON",ind_obs_nino2)
    dssta_djf2 = cal_composite(dsst,"DJF",ind_obs_nino2)

    dssta_mam2 =(/dssta_mam2 -dssta_mam1/)
    dssta_jja2 =(/dssta_jja2 -dssta_jja1/)
    dssta_son2 =(/dssta_son2 -dssta_son1/)
    dssta_djf2 =(/dssta_djf2 -dssta_djf1/)
    
    ;  dssta_mam1 = cal_composite(ssta,"MAM",ind_obs_nino1)
    ; dssta_jja1 = cal_composite(ssta,"JJA",ind_obs_nino1)
    ; dssta_son1  = cal_composite(ssta,"SON",ind_obs_nino1)
    ; dssta_djf1 = cal_composite(ssta,"DJF",ind_obs_nino1)

    ; dssta_mam2 = cal_composite(ssta,"MAM",ind_obs_nino2)
    ; dssta_jja2 = cal_composite(ssta,"JJA",ind_obs_nino2)
    ; dssta_son2  = cal_composite(ssta,"SON",ind_obs_nino2)
    ; dssta_djf2 = cal_composite(ssta,"DJF",ind_obs_nino2)
    
    ; prob_mam = rc_test(ssta_mam)
    ; prob_jja = rc_test(ssta_jja)

    ; prob_son = rc_test(ssta_son)
    ; prob_djf = rc_test(ssta_djf)

    ; wks   = gsn_open_wks("pdf", "global季节平均SST的趋势-philipian-79-19")
    ; wks   = gsn_open_wks("pdf", "ERSST-Eninoyear-DSST-evolution-79-19")
        wks   = gsn_open_wks("pdf", "ERSST-Eninoyear-dSSTa-evolution-79-19")
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
  res@lbLabelBarOn=True
  res@lbOrientation="Horizontal"
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
  ; res@cnMinLevelValF       =   -0.4
  ; res@cnMaxLevelValF       =   0.4
  ; res@cnLevelSpacingF      =  0.1

  ;   res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  0.2
   res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; res@cnLevels := (/-4,-3,-2,-1,1,2,3,4/)*0.1
  res@cnLevels := (/-2,-1.5,-1,-0.5,0.5,1,1.5,2/)*0.1
  res@cnInfoLabelOn = False
  

  res@gsnLeftString="MAM SSTA tendency ~S~o~N~C/month"
  res@gsnLeftStringFontHeightF = 0.02

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

  plot=new(8,graphic)

  plot(0) = gsn_csm_contour_map(wks, dssta_mam1, res)
    res@gsnLeftString="JJA SSTA tendency ~S~o~N~C/month"
  plot(2) = gsn_csm_contour_map(wks, dssta_jja1, res)
    res@gsnLeftString="SON SSTA tendency ~S~o~N~C/month"
  plot(4) = gsn_csm_contour_map(wks, dssta_son1, res)
    res@gsnLeftString="DJF SSTA tendency ~S~o~N~C/month"
  plot(6) = gsn_csm_contour_map(wks, dssta_djf1, res)

 res@gsnLeftString="MAM SSTA tendency ~S~o~N~C/month"
  plot(1) = gsn_csm_contour_map(wks, dssta_mam2, res)
    res@gsnLeftString="(after 2000)-(before 2000) JJA SSTA tendency ~S~o~N~C/month"
  plot(3) = gsn_csm_contour_map(wks, dssta_jja2, res)
    res@gsnLeftString="SON SSTA tendency ~S~o~N~C/month"
  plot(5) = gsn_csm_contour_map(wks, dssta_son2, res)
    res@gsnLeftString="DJF SSTA tendency ~S~o~N~C/month"
  plot(7) = gsn_csm_contour_map(wks, dssta_djf2, res)



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

    str2 = ispan(1, 8, 1) + ""
     plot_add2 = new(8,graphic)



    do i = 0,7
      plot_add2(i) = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
    end do 

   ; ;;;打点
   ;  res_pro=True
   ;  res_pro@gsnAddCyclic = False
   ;  res_pro@gsnDraw=False
   ;  res_pro@gsnFrame=False
   ;  res_pro@cnFillOn = True

   ;  res_pro@cnLevelSelectionMode = "ExplicitLevels"
   ;  res_pro@cnLevels = (/0,0.1/)
   ;  res_pro@lbLabelBarOn = False
   ;  res_pro@cnLinesOn=False
   ;  ; res_pro@cnMonoLineColor =False
   ;  ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

   ;  res_pro@cnMonoFillPattern = True
   ;  res_pro@cnFillPattern = 17
   ;  res_pro@cnFillDotSizeF = 0.001
   ;  res_pro@cnMonoFillColor = False
   ;  res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
   ;  ; res_pro@cnLevels = ispan(-6,6,1)
   ;  res_pro@cnLineLabelsOn=False
   ;  res_pro@tmYRMode = "Automatic"
   ;  res_pro@cnInfoLabelOn=False
   ;  ; res_pro@cnLineLabelBackgroundColor="Transparent"
   ;  ; res_pro@cnLineThicknessF = 3.
   ;  ; res_pro@cnFillOn        = True
  
   ;   iplot = new(4,graphic)
   ;   iplot(0) = gsn_csm_contour(wks,prob_mam  ,res_pro)
   ;   iplot(1) = gsn_csm_contour(wks,prob_jja  ,res_pro)
   ;   iplot(2) = gsn_csm_contour(wks,prob_son  ,res_pro)
   ;   iplot(3) = gsn_csm_contour(wks,prob_djf  ,res_pro)
   ;  ;  
   ;  do i=0,3
   ;   overlay(plot(i), iplot(i))
   ;  end do

pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.008
pres@lbLabelStride=1


draw(plot(3))
frame(wks)
; gsn_panel(wks, (/plot(1),plot(3),plot(5),plot(7)/), (/2,2/), pres)

end