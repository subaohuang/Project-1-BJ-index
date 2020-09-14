;;;画ERA 的0-12小时forecast月平均的通量数据的趋势
;;
;;
;;                                         lnshuheng 2020.03.10
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
;;;;

function change(lhtfl,season)
begin 
    ; lhtfl_annual = month_to_season(lhtfl, season)
    lhtfl_annual = month_to_annual(lhtfl, 1)
    lhtfl_change = lhtfl_annual(1,:,:)
    dim  = dimsizes(lhtfl_annual)
    rc = regCoef_n(ispan(1, dim(0), 1), lhtfl_annual,0,0)
    lhtfl_change = (/rc*dim(0)/)
return lhtfl_change
end     

;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   ; lat1 = -5
   ; lat2 = 45
   ; lon1 = 90
   ; lon2 = 160
   lat1 = -10
   lat2 = 20
   lon1 = 100
   lon2 = 160

   lat1 = -30
   lat2 = 30
   lon1 = 0
   lon2 = 360
;;;;


;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.198001)

    sst_r = lonFlip(f1->sst)  

    sst = sst_r(timeind,{lat1:lat2},{lon1:lon2})
    
;;;read NCEP-NCAR heat flux data
    fils3 = "/home/ys17-19/lsh/Project/SCS-rain/annual/ERAI-heatflux-1degree.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.198001)
  
    
    lhtfl = f3->lhtfl(timeind,{lat1:lat2},{lon1:lon2})
    shtfl = f3->shtfl(timeind,{lat1:lat2},{lon1:lon2})
    nswrs = f3->nswrs(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs = f3->nlwrs(timeind,{lat1:lat2},{lon1:lon2})

    lhtfl = (/-1*lhtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2   ;;;convert to positive 
    shtfl = (/-1*shtfl/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2    ;;;convert to positive
    nswrs = (/nswrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2   
    nlwrs = (/-1*nlwrs/(24*3600)/)  ;;units j/m^2 (12hour) to w/m^2   ;;;convert to positive
   ;;;;


  ;;;;choose season

  season = (/"MAM","JJA","SON","DJF"/)

  wks   = gsn_open_wks("pdf", "ATM-ocean-forcing-SST-annual-global")

  plot = new(12,graphic)
  str = "" + ispan(1,dimsizes(plot),1)
  str2 = "0" + ispan(1,dimsizes(plot),1)

 i = 0
    sst_change = change(sst,season(i))
    printVarSummary(sst_change)
    lhtfl_change = change(lhtfl,season(i))
    shtfl_change = change(shtfl,season(i))
    nswrs_change = change(nswrs,season(i))
    nlwrs_change = change(nlwrs,season(i))

    Do_change  = lhtfl_change
    Do_change = (/-1*( nswrs_change - nlwrs_change - lhtfl_change -shtfl_change)/) ;;;;Do = -Qnet
    

    lhtfl_clm  = dim_avg_n_Wrap(lhtfl, 0)
    printVarSummary(lhtfl_clm)
    
    lhtfl_change_atm = lhtfl_clm
    lhtfl_change_atm  = (/lhtfl_change - 0.06*lhtfl_clm*sst_change/)

    Qa_change = lhtfl_change_atm 
    Qa_change = (/ nswrs_change - nlwrs_change -shtfl_change - lhtfl_change_atm/)


    T_change_atm = Do_change
    T_change_atm = (/Qa_change/(0.06*lhtfl_clm)/)

    T_change_ocean = Do_change
    T_change_ocean = (/Do_change/(0.06*lhtfl_clm)/)
    
    T_total = T_change_ocean
    T_total = (/T_change_atm +T_change_ocean/)
;;;;
   
print("*********开始画图********")
;;**************************************************************
;画图
;;***************************************************************


  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  ; res@cnFillPalette="BlueWhiteOrangeRed"
  res@cnFillPalette = "NCV_jaisnd"
  res@gsnSpreadColorStart = 130
  res@gsnSpreadColorEnd = 253

  res@cnFillOn=True
  res@cnLinesOn=True
  res@cnLineLabelInterval = 2
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillDrawOrder = "PostDraw"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
   lat1 = 0
   lat2 = 15
   lon1 = 125
   lon2 = 160

   ; lat1 = -30
   ; lat2 = 30
   ; lon1 = 0
   ; lon2 = 360

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
  res@gsnRightString = ""

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   0.0
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.1





  printVarSummary(sst_change)
  printVarSummary(T_change_atm)
  printVarSummary(T_change_ocean)

  res@gsnLeftString = "SST trend "
  plot(i*3) = gsn_csm_contour_map(wks, sst_change, res)
  
  delete(res@gsnSpreadColorStart)
  delete(res@gsnSpreadColorEnd)
  res@cnFillPalette = "BlRe"
  res@cnMinLevelValF       =   -3.0
  res@cnMaxLevelValF       =   3.0
  res@cnLevelSpacingF      =  0.5


  res@gsnLeftString = "Atm related SST trend"
  plot(i*3+1) = gsn_csm_contour_map(wks, T_change_atm, res)
  res@gsnLeftString = "Ocean related SST trend"
  plot(i*3+2) = gsn_csm_contour_map(wks, T_change_ocean, res)


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
   


     plot@$str(i*3)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
     plot@$str(i*3+1)$ = gsn_add_polyline(wks, plot(i+1), pgx, pgy, pgres)
     plot@$str(i*3+2)$ = gsn_add_polyline(wks, plot(i+2), pgx, pgy, pgres)




  ; ;;;添加indo-pacific warm pool
  ;   ; ;;添加区域
  ;  ;;warm pool
  ;  pgres@gsLineThicknessF=2.0
  ;  pgres@gsLineDashPattern=0
  ;  str3=unique_string("dum3")
  ;  str4=unique_string("dum4")
   
  ;  ;;;画box
  ;  lat1 = 0
  ;  lat2 = 15
  ;  lon1 = 110
  ;  lon2 = 155

  ;  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  ;  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;  ;;;
  ; plot@$str2(i*3)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  ; plot@$str2(i*3+1)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  ; plot@$str2(i*3+2)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)


pres = True
pres@gsnPanelLabelBar = False

gsn_panel(wks, (/plot(0),plot(1),plot(2)/), (/1,3/), pres)

 


end