;;;;;;
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
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;

;;;
    fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
    f3 = addfiles(fils3, "r")

    date=cd_calendar(f3[:]->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)
    time = f3[:]->time(timeind)
    lat = f3[1]->latitude
    lon = f3[1]->longitude
  
    ;;;region 
      lat1 = -90
      lat2 = 90
      lon1 = 0
      lon2 = 360
    ;;;

      
    ; div_r = new((/dimsizes(time),7, dimsizes(lat),dimsizes(lon)/),float)
      
    ; do i = 0,dimsizes(timeind)/12-1
    ; div_r(i*12:(i+1)*12-1,:,:,:) = short2flt(f3[i]->d(:,{1000:850},:,:))
    ; end do
    
    div_r = new((/dimsizes(time), dimsizes(lat),dimsizes(lon)/),float)
      
    do i = 0,dimsizes(timeind)/12-1
    div_r(i*12:(i+1)*12-1,:,:) = short2flt(f3[i]->d(:,{200},:,:))
    end do
    
    

    ; div = dim_avg_n_Wrap(div_r(:,:,{-10:30},{100:180}),1)
    div = div_r
    ; div = dim_avg_n_Wrap(div_r,1)

    div = (/div*10^5/)
    div = (/-1*div/)


    div_mam_clm = dim_avg_n_Wrap(month_to_season(div, "MAM"), 0)
    div_jja_clm = dim_avg_n_Wrap(month_to_season(div, "JJA"), 0)
    div_son_clm = dim_avg_n_Wrap(month_to_season(div, "SON"), 0)
    div_djf_clm = dim_avg_n_Wrap(month_to_season(div, "DJF"), 0)




    diva=(/rmMonAnnCycTLL(div)/)

    copy_VarCoords(div, diva)
    

   ;;;;cal season trend 
    diva_mam = cal_season_trend(diva,"MAM")
    diva_jja = cal_season_trend(diva,"JJA")
    diva_son  = cal_season_trend(diva,"SON")
    diva_djf = cal_season_trend(diva,"DJF")
    ;;;;

  ; wks   = gsn_open_wks("pdf", "global季节平均SST的趋势-philipian-79-16")
    ; wks   = gsn_open_wks("pdf", "西太低层风场幅合季节平均div的气候态-79-15-850")
    wks   = gsn_open_wks("pdf", "西太低层风场幅合季节平均div的气候态-79-15-200")
    ; wks   = gsn_open_wks("pdf", "全球风场幅合季节平均div的气候态-79-15")
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@gsnAddCyclic = False



  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueDarkRed18"
  ; res@cnFillPalette = "BlRe"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=True

  res@cnLineLabelsOn=True
  res@lbLabelBarOn = False
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



  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"



  res@tmXBLabelFontHeightF=0.025
  res@tmYLLabelFontHeightF=0.025
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0

  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 

  res@gsnLeftStringFontHeightF = 0.03


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levelW
  res@cnMinLevelValF       =   -0.3
  res@cnMaxLevelValF       =   .3
  res@cnLevelSpacingF      =  0.05
  ; res@cnFillPalette = "NCV_jaisnd"
  ; res@cnFillPalette = "sunshine_9lev"
  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; res@cnLevels = (/-2.5,-1.5,-0.5,0.5,1.5,2.5/)
  ; ; res@cnFillColors= (/4,5,6,7,8,9,10,11,12,13,14/) - 2

  res@cnInfoLabelOn = False


  res@gsnLeftString="(a)MAM Div clm"
  ; res@gsnRightString = "10~S~-6~N~S^-1"


  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

  plot=new(4,graphic)

  plot(0) = gsn_csm_contour_map(wks, div_mam_clm, res)
  res@gsnLeftString="(b) JJA Div clm"
  plot(1) = gsn_csm_contour_map(wks, div_jja_clm, res)
  res@gsnLeftString="(c) SON Div clm"
  plot(2) = gsn_csm_contour_map(wks, div_son_clm, res)
  res@gsnLeftString="(d) DJF Div clm"
  plot(3) = gsn_csm_contour_map(wks, div_djf_clm, res)

 pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=0

  str = ispan(1, 4, 1) + ""
    plot_add = new(4,graphic)
  ; do i = 0,3
  ; plot_add(i) = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
  ;  end do 
  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
;;画box
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


pres = True
pres@gsnPanelLabelBar  = True
; pres@gsnMaximize = True
; pres@lbOrientation = "Vertical"
pres@lbOrientation = "Horizontal"
; pres@pmLabelBarHeightF = 0.3
pres@lbLabelFontHeightF = 0.02
; pres@lbLabelStride = 2
gsn_panel(wks, plot, (/2,2/), pres)

end