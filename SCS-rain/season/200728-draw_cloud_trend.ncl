   ;********************************************************
;;画HadIsst的sst趋势
;;
;;lnshuheng 2020.02.17
;**********************************************************
function  rc_test(cldc_rc)
begin
  tval = onedtond(cldc_rc@tval , dimsizes(cldc_rc))
   df   = onedtond(cldc_rc@nptxy, dimsizes(cldc_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(cldc_rc, prob)

return prob
end

function cal_season_trend(data,season)
  begin
    p_djf_r=month_to_season(data, season)
      if season.eq."DJF"
      p_djf_r(0:dimsizes(p_djf_r&time)-2,:,:) = p_djf_r(1:,:,:)
      p_djf_r(dimsizes(p_djf_r&time)-1,:,:) = data(dimsizes(data&time)-1,:,:)
    end if
    dim = dimsizes(p_djf_r)
    p_rc=regCoef_n(ispan(1, dim(0), 1), p_djf_r, 0, 0)
    copy_VarCoords(p_djf_r(1,:,:), p_rc)
    p_rc=(/p_rc*dim(0)/)
  return p_rc
end 

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
    yearin  = ispan(1979, 2014, 1)
    fils4 = "/home/yangsong3/data-observation/ERA-interim/monthly/cloud/cloud_monthly_" + yearin + ".nc"
    f4 = addfiles(fils4, "r")

    date:=cd_calendar(f4[:]->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    time = f4[:]->time(timeind)
    lat = f4[1]->latitude
    lon = f4[1]->longitude
  

    cldc = new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),float)
   
    do i = 0,dimsizes(timeind)/12-1
      cldc(i*12:(i+1)*12-1,:,:) = short2flt(f4[i]->hcc(:,:,:))
    end do
    ; ;;;插值，missing value值给去掉，只能限制longitude的5点和纬向的1点
    ; opt=(/0,5/)
    ; ;;先插值lon方向的
    ; cldc=(/linmsg_n(cldc, opt, 2)/)
    ; ;;插值lat方向的
    ; opt=(/0,1/)
    ; cldc=(/linmsg_n(cldc, opt, 1)/)
    ; delete(opt)
    ; cldc_annual  = month_to_season(cldc(:,{lat1:lat2},{lon1:lon2}),season)
    ;;;

    ; cldc_annual&year = ispan(1, 34, 1)
    ; printVarSummary(cldc_annual)

    ; cldc_annual_tropical=dim_avg_n_Wrap(month_to_annual(cldc(:,{-30:30},:),1), (/1,2/))
    ; do i=0,dimsizes(cldc_annual&year)-1
    ;     cldc_annual(i,:,:)=(/cldc_annual(i,:,:)-cldc_annual_tropical(i)/)
    ; end do
    ; ;;
;;;
 ;;;region
     lat1 = 0
     lat2 = 15
     lon1 = 110
     lon2 = 155
     ;;;


;;;;cal season trend 
    cldc_mam = cal_season_trend(cldc(:,{lat1:lat2},{lon1:lon2}),"MAM")
    cldc_jja = cal_season_trend(cldc(:,{lat1:lat2},{lon1:lon2}),"JJA")
    cldc_son  = cal_season_trend(cldc(:,{lat1:lat2},{lon1:lon2}),"SON")
    cldc_djf = cal_season_trend(cldc(:,{lat1:lat2},{lon1:lon2}),"DJF")



    prob_mam = rc_test(cldc_mam)
    prob_jja = rc_test(cldc_jja)

    prob_son = rc_test(cldc_son)
    prob_djf = rc_test(cldc_djf)

    ; cldc_mam = where(prob_mam.le.0.05, cldc_mam, cldc_mam@_FillValue)
    ; cldc_jja = where(prob_jja.le.0.05, cldc_jja, cldc_jja@_FillValue)
    ; cldc_son = where(prob_son.le.0.05, cldc_son, cldc_son@_FillValue)
    ; cldc_djf = where(prob_djf.le.0.05, cldc_djf, cldc_djf@_FillValue)
    

  
    wks   = gsn_open_wks("pdf", "季节平均high-cloud的趋势-SCS-79-15")
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@gsnAddCyclic = False
  ; res@cnFillPalette = "NCV_jaisnd"
   res@cnFillPalette = "BlRe"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -60
  ; res@mpMaxLatF             = 60
  ; res@mpMinLonF             =0
  ; res@mpMaxLonF             = 360

  res@mpMinLatF             = 0
  res@mpMaxLatF             = 15
  res@mpMinLonF             =110
  res@mpMaxLonF             = 155



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

  res@cnMinLevelValF       =   -0.05
  res@cnMaxLevelValF       =   0.05
  res@cnLevelSpacingF      =  0.01
  res@cnInfoLabelOn = False
  

  res@gsnLeftString="Spring"

  plot=new(4,graphic)

  plot(0) = gsn_csm_contour_map(wks, cldc_mam, res)
    res@gsnLeftString="Summer"
  plot(1) = gsn_csm_contour_map(wks, cldc_jja, res)
    res@gsnLeftString="Autumn"
  plot(2) = gsn_csm_contour_map(wks, cldc_son, res)
    res@gsnLeftString="Winter"
  plot(3) = gsn_csm_contour_map(wks, cldc_djf, res)

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

  plot1(0) = gsn_csm_contour(wks, cldc_mam, res2)
  plot1(1) = gsn_csm_contour(wks, cldc_jja, res2)
  plot1(2) = gsn_csm_contour(wks, cldc_son, res2)
  plot1(3) = gsn_csm_contour(wks, cldc_djf, res2)

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
  




