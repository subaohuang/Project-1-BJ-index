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

undef("remove_enso_signal")
function remove_enso_signal(data)
  begin  
  ; ;;;读入HadiSST海温数据
      diri="/home/yangsong3/data-observation/linshh/data/sst/"
      fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
      f1=addfile(fils1, "r")
      date:=cd_calendar(f1->time, 1)

      timeind:=ind(date.le.201512.and.date.ge.197901)

      sst =  lonFlip(f1->sst(timeind,:,:))
   
      ssta=(/rmMonAnnCycTLL(sst)/)
      copy_VarCoords(sst, ssta)

      ssta_de = ssta
      ; ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
      ssta_index = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
      printVarSummary(ssta_index)

      nino3_4 = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
      printVarSummary(nino3_4)

      
      ; nino3_4_new = nino3_4                ;;;;前期ENSO
      ; nino3_4_new(0:33) = (/nino3_4(1:)/)
      ; nino3_4_new(34) =ssta_index(35*12-1) 
   
      nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
  ;;;;;
      data_de = data
      ; data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_rc = regCoef_n(nino3_4,data_de, 0, 0)
      copy_VarCoords(data(1,:,:), data_rc)
      data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
      m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
      nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))
  
      data_relate_enso  = nino3_4_3d*data_rc_3d + m 
      data_new = data 
      data_new = (/data - data_relate_enso/)
  return data_new
end 

function cal_season_trend(data,season)
  begin
   
    p_djf_r=month_to_season(data, season)
     ; ;;;;remove_enso_signal
     ;     p_djf_r = (/remove_enso_signal(p_djf_r)/)
     ; ;;;;

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


function eot(sst_mam_clm)          ;;;; -14.0275*x+375.845  (26.75~28.75) 25.2828*x+-758.963(28.75~30.25)                                     
 begin 
  eot_mam = sst_mam_clm
  eot_mam = 0.0
  eot_mam =  where(sst_mam_clm.lt.26.75, 0.0, eot_mam)
  eot_mam =  where(sst_mam_clm.ge.26.75.and.sst_mam_clm.le.28.75, -14.0275*sst_mam_clm + 375.845 , eot_mam)
  eot_mam =  where(sst_mam_clm.gt.28.75.and.sst_mam_clm.le.30.25,  25.2828*sst_mam_clm -758.963, eot_mam)
  return eot_mam
 end 


function eot_cdf(sst_mam_clm,ssta,regionidx)          ;;;; -14.0275*x+375.845  (26.75~28.75) 25.2828*x+-758.963(28.75~30.25)                                     
 begin 

  eot_mam_cdf= sst_mam_clm
  eot_mam_cdf = 0.0
  
  eot_mam_cdf2 = sst_mam_clm
  eot_mam_cdf2= 0.0

do i = 0,2 
    region = regionidx(i)

    if region.eq."twnp" then 
    lat1 = -5
    lat2 = 20
    lon1 = 120
    lon2 = 170
    a = -14.0275
    b = 375.845
    c = 25.2828
    d = -758.963
    ;  a = -13.1
    ; b = 349.0
    ; c = 32.2
    ; d = -956.4
    range = (/26.5,28.75,29.75/)

    else  if region.eq."wnp" then    ;;;;Wnp region 100-170 20-30
    lat1 = 20
    lat2 = 30
    lon1 = 100
    lon2 = 170
    a = -10.1639
    b = 276.726
    c = 20.1925
    d = -612.912
    range = (/27.0,29.25,30.25/)

    else  if region.eq."scs" then    ;;;;Wnp region 100-170 20-30

    lat1 = -5
    lat2 = 20
    lon1 = 100
    lon2 = 120
    ; a = -5.32061
    ; b = 136.234
    ; c = 11.5947
    ; d = -349.058
    ; range = (/25.25,28.5,30.25/)

    ; a = -10.1639
    ; b = 276.726
    ; c = 20.1925
    ; d = -612.912
    ; range = (/27.0,29.25,30.25/)
    a = -14.0275
    b = 375.845
    c = 25.2828
    d = -758.963
    ;  a = -13.1
    ; b = 349.0
    ; c = 32.2
    ; d = -956.4
    range = (/26.5,28.75,29.75/)
    end if 
    end if 
    end if 
    ;;;;clm SST induced OLR 
    eot_mam_cdf({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).lt.range(0), 0.0, eot_mam_cdf({lat1:lat2},{lon1:lon2}))
    eot_mam_cdf({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).ge.range(0).and.sst_mam_clm({lat1:lat2},{lon1:lon2}).le.range(1), (a/2.0)*sst_mam_clm({lat1:lat2},{lon1:lon2})^2 + b*sst_mam_clm({lat1:lat2},{lon1:lon2}) , eot_mam_cdf({lat1:lat2},{lon1:lon2}))
    eot_mam_cdf({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).gt.range(1).and.sst_mam_clm({lat1:lat2},{lon1:lon2}).le.range(2),  (c/2.0)*sst_mam_clm({lat1:lat2},{lon1:lon2})^2 +d*sst_mam_clm({lat1:lat2},{lon1:lon2}), eot_mam_cdf({lat1:lat2},{lon1:lon2}))

    ;;;; change sst induced OLR 
    eot_mam_cdf2({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).lt.range(0), 0.0, eot_mam_cdf2({lat1:lat2},{lon1:lon2}))
    eot_mam_cdf2({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).ge.range(0).and.sst_mam_clm({lat1:lat2},{lon1:lon2}).le.range(1).and.(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2})).le.range(1), (a/2.0)*(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2}))^2 + \
                                                   b*(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2})) , eot_mam_cdf2({lat1:lat2},{lon1:lon2}))
    eot_mam_cdf2({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).ge.range(0).and.sst_mam_clm({lat1:lat2},{lon1:lon2}).le.range(1).and.(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2})).gt.range(1), (a/2.0)*(range(1))^2 + b*(range(1)) + \
     (c/2.0)*(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2}))^2 +d*(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2})) - ((c/2.0)*(range(1))^2 + d*(range(1))) , eot_mam_cdf2({lat1:lat2},{lon1:lon2}))

    eot_mam_cdf2({lat1:lat2},{lon1:lon2}) =  where(sst_mam_clm({lat1:lat2},{lon1:lon2}).gt.range(1).and.sst_mam_clm({lat1:lat2},{lon1:lon2}).le.range(2),  (c/2.0)*(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2}))^2 +d*(sst_mam_clm({lat1:lat2},{lon1:lon2})+ssta({lat1:lat2},{lon1:lon2})), eot_mam_cdf2({lat1:lat2},{lon1:lon2}))


  delete(region)
  delete([/lat1,lat2,lon1,lon2,a,b,c,d/])
  end do 
  eot_mam_cdf_final = eot_mam_cdf
  eot_mam_cdf_final = (/eot_mam_cdf2 - eot_mam_cdf/)

  return eot_mam_cdf_final
 end 


;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.197901)

    sst  = f1->sst(timeind,{-5:30},{100:170})

    sst_mam_clm = dim_avg_n_Wrap(month_to_season(sst, "MAM"), 0)


    sst_jja_clm = dim_avg_n_Wrap(month_to_season(sst, "JJA"), 0)
    sst_son_clm = dim_avg_n_Wrap(month_to_season(sst, "SON"), 0)
    sst_djf_clm = dim_avg_n_Wrap(month_to_season(sst, "DJF"), 0)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    

   ;;;;cal season trend 
    ssta_mam = cal_season_trend(ssta,"MAM")
    ssta_jja = cal_season_trend(ssta,"JJA")
    ssta_son  = cal_season_trend(ssta,"SON")
    ssta_djf = cal_season_trend(ssta,"DJF")
   ;;;;

    sst_plus_mam  = sst_mam_clm 
    sst_plus_jja = sst_jja_clm
    sst_plus_son = sst_son_clm
    sst_plus_djf = sst_djf_clm

    sst_plus_mam = (/sst_plus_mam + ssta_mam/)
    sst_plus_jja = (/sst_plus_jja + ssta_jja/)
    sst_plus_son = (/sst_plus_son + ssta_son/)
    sst_plus_djf = (/sst_plus_djf + ssta_djf/)
    


    ; ssta_jja = (/ssta_jja - ssta_mam/)
    ; ssta_son = (/ssta_son- ssta_mam/)
    ; ssta_djf = (/ssta_djf -ssta_mam/)
    ; ssta_mam = 0.0

    prob_mam = rc_test(ssta_mam)
    prob_jja = rc_test(ssta_jja)

    prob_son = rc_test(ssta_son)
    prob_djf = rc_test(ssta_djf)

    ; eot_mam = eot(sst_mam_clm)
    ; eot_jja = eot(sst_jja_clm)   
    ; eot_son = eot(sst_son_clm)
    ; eot_djf = eot(sst_djf_clm)

    
    olr_rela_mam  = ssta_mam
    olr_rela_jja  = ssta_mam
    olr_rela_son  = ssta_mam
    olr_rela_djf  = ssta_mam


    ; olr_rela_mam = (/eot_mam*ssta_mam/)
    ; olr_rela_jja = (/eot_jja*ssta_jja/)
    ; olr_rela_son = (/eot_son*ssta_son/)
    ; olr_rela_djf = (/eot_djf*ssta_djf/)

    regionidx_want = (/"twnp","wnp","scs"/)

    olr_rela_mam = (/eot_cdf(sst_mam_clm,ssta_mam,regionidx_want)/)
    olr_rela_jja = (/eot_cdf(sst_jja_clm,ssta_jja,regionidx_want)/)
    olr_rela_son = (/eot_cdf(sst_son_clm,ssta_son,regionidx_want)/)
    olr_rela_djf = (/eot_cdf(sst_djf_clm,ssta_djf,regionidx_want)/)


    print(avg(olr_rela_mam({0:15},{125:160})))
    print(avg(olr_rela_jja({0:15},{125:160})))
    print(avg(olr_rela_son({0:15},{125:160})))
    print(avg(olr_rela_djf({0:15},{125:160})))

  ; out = addfile("olr_change_related_2_SST.nc","c")
  ; out->olr_rela_mam  = olr_rela_mam
  ; out->olr_rela_jja  = olr_rela_jja
  ; out->olr_rela_son  = olr_rela_son
  ; out->olr_rela_djf  = olr_rela_djf
  ; exit()


    print(avg(olr_rela_mam({0:15},{125:160})))


    ; wks   = gsn_open_wks("pdf", "global季节平均SST的趋势-philipian-79-16")
    wks   = gsn_open_wks("pdf", "OLR-change-related-to-SSTchange-79-15-two-region")
;;**************************************************************
;画图
;;***************************************************************

    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic = False


    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"
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
    ; res@mpMinLatF             = -30
    ; res@mpMaxLatF             = 30
    ; res@mpMinLonF             =0
    ; res@mpMaxLonF             = 360

    res@mpMinLatF             = -5
    res@mpMaxLatF             = 30
    res@mpMinLonF             = 100
    res@mpMaxLonF             = 170

    ; res@mpMinLatF             = 0
    ; res@mpMaxLatF             = 15
    ; res@mpMinLonF             = 125
    ; res@mpMaxLonF             = 160

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

    res@gsnLeftStringFontHeightF = 0.025


  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -12.0
  ; res@cnMaxLevelValF       =   12.0
  ; res@cnLevelSpacingF      =  2.0
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnFillPalette= "BlueWhiteOrangeRed"
    res@cnLevels := (/-12.0,-10.0,-8.0,-6.0,-4.0,-2.0,2.0,4.0,6.0,8.0,10.0,12.0/)
 


    res@cnInfoLabelOn = False


    res@gsnLeftString="(a)  MAM SST warming induced OLR change"
    res@gsnRightString = ""


  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

   plot=new(4,graphic)


  ; plot(0) = gsn_csm_contour_map(wks, olr_rela_mam({-0.5:15.5},{124.5:160.5}), res)
  ; res@gsnLeftString="(b) JJA OLR change related to SST change"
  ; plot(1) = gsn_csm_contour_map(wks, olr_rela_jja({-0.5:15.5},{124.5:160.5}), res)
  ; res@gsnLeftString="(c) SON OLR change related to SST change"
  ; plot(2) = gsn_csm_contour_map(wks, olr_rela_son({-0.5:15.5},{124.5:160.5}), res)
  ; res@gsnLeftString="(d) DJF OLR change related to SST change"
  ; plot(3) = gsn_csm_contour_map(wks, olr_rela_djf({-0.5:15.5},{124.5:160.5}), res)

 plot(0) = gsn_csm_contour_map(wks, olr_rela_mam, res)
  res@gsnLeftString="(b) JJA SST warming induced OLR change"
  plot(1) = gsn_csm_contour_map(wks, olr_rela_jja, res)
  res@gsnLeftString="(c) SON SST warming induced OLR change"
  plot(2) = gsn_csm_contour_map(wks, olr_rela_son, res)
  res@gsnLeftString="(d) DJF SST warming induced OLR change"
  plot(3) = gsn_csm_contour_map(wks, olr_rela_djf, res)


  res2 = res
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnFillOn = False
  res2@cnLinesOn = True
  res2@cnLevels := fspan(26.0, 30.0, 9)
  res2@cnLineLabelStrings =   res2@cnLevels  + ""
  res2@cnLineLabelDensityF = 1.3


  res2@cnLineColors = "black"
  res2@cnMonoLineColor  = True
  res2@cnLineThicknessF = 2
  res2@cnMonoLineThickness  = True
  res2@gsnLeftString = ""
  res2@cnInfoLabelOn = False
  res2@cnLineLabelsOn = True
  res2@cnLineLabelInterval = 1
  res2@cnLineLabelBackgroundColor = "White"
  
  plot1=new(4,graphic)

  plot1(0) = gsn_csm_contour(wks, sst_mam_clm, res2)
  plot1(1) = gsn_csm_contour(wks, sst_jja_clm, res2)
  plot1(2) = gsn_csm_contour(wks, sst_son_clm, res2)
  plot1(3) = gsn_csm_contour(wks, sst_djf_clm, res2)

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
  ; do i = 0,3
  ; plot_add(i) = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
  ;  end do 
  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF = 3.0
   pgres@gsLineDashPattern= 1
   pgres@gsLineColor = "black"
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
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  
     iplot = new(4,graphic)
     iplot(0) = gsn_csm_contour(wks,prob_mam  ,res_pro)
     iplot(1) = gsn_csm_contour(wks,prob_jja  ,res_pro)
     iplot(2) = gsn_csm_contour(wks,prob_son  ,res_pro)
     iplot(3) = gsn_csm_contour(wks,prob_djf  ,res_pro)
    ;  
    do i=0,3
     overlay(plot(i), iplot(i))
    end do

pres = True
pres@gsnPanelLabelBar  = True
pres@gsnMaximize = True
; pres@lbOrientation = "Vertical"
pres@lbOrientation = "Horizontal"
; pres@pmLabelBarHeightF = 0.3
pres@lbLabelStrings =  (/-12.0,-10.0,-8.0,-6.0,-4.0,-2.0,2.0,4.0,6.0,8.0,10.0,12.0/)

pres@lbLabelFontHeightF = 0.015
pres@lbAutoManage = False
pres@lbLabelAutoStride =False
pres@lbLabelStride = 1
; pres@lbLabelStride = 2
gsn_panel(wks, plot, (/2,2/), pres)

end