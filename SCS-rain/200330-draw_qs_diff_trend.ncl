;;********************************************************
;;计算ERA的10-m风趋势
;;
;;lnshuheng 2020.03.17
;**********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/crop.ncl"
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
begin

;==========================================================
;;;read data
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/NCEP_DOE/monthly/shum.2m.mon.mean.nc"
    
  
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201312.and.date.ge.198001)

    sst=lonFlip(f1->sst(timeind,{-90:90},{-180:180}))
    printVarSummary(sst)

   

    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
    
   ;;;region
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155
   ;;;


    ;;;time span
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    ;;


    f1 = addfile(path_wind,"r")

    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)

    shum_r = short2flt(f1->shum(timeind,0,{-90:90},:))

  
    print("=====finish read data==========")

    

    print("*****************开始水平插值************************")  
    xo  =  sst&longitude
    yo  =  sst&latitude
    shum = linint2_Wrap(shum_r&lon, shum_r&lat, shum_r, \
                  True, xo, yo, 0)

    
    printVarSummary(shum)
   ;;;

  
  ;;;计算大气的饱和比湿
     es = sst
     es =satvpr_temp_fao56(sst,(/0,2/))             ;;;先计算饱和水气压,输出的单位是kpa
     p = 100 ;; 100kpa,地面气压
     qs = es
   
     qs =   0.622*es/(p - 0.378*es)
     printVarSummary(qs)
     
     
     q_diff = shum
     printVarSummary(shum)
     q_diff = (/shum - qs/)


   ;;;处理成年数据 
      q_diff_annual = month_to_annual(q_diff(:,{lat1:lat2},{lon1:lon2}), 1)
   ;;;

    
   ;;趋势
   q_rc = regCoef_n(ispan(1, yearnum, 1), q_diff_annual, 0, 0)
   copy_VarCoords(q_diff_annual(1,:,:), q_rc)
   q_rc = (/q_rc*yearnum/)
   ;;;
   delete(shum)


    q_diff_scs = dim_avg_n_Wrap(q_diff_annual(:,:,:), (/1,2/))
    wks = gsn_open_wks("pdf", "era-q随时间的变化")
    res = True
    plot= gsn_csm_xy(wks, ispan(1980,2013,1), q_diff_scs, res)
    draw(plot)
    frame(wks)
    exit()



   ;;;显著性检验
    
    prob_u = rc_test(q_rc)
 

   ;;;
 
   ; d_rc = where(prob_d.le.0.1, d_rc, d_rc@_FillValue)
   ; sta_s_rc = where(prob_sta_s.le.0.1, sta_s_rc, sta_s_rc@_FillValue)
   ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "绝对湿度差在南海区域的趋势")

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlRe"

  res@cnFillOn=True
  res@cnLinesOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -80.0
  ; res@cnMaxLevelValF       =   80.0
  ; res@cnLevelSpacingF      =  10.0
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.001
  res@cnMaxLevelValF       =   0.001
  res@cnLevelSpacingF      =  0.0001

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
  
  res@mpFillOn =False

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.8
  plot = new(1, graphic)
  plot = gsn_csm_contour_map(wks,q_rc,res)

  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
    i =0
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   

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
    iplot = new(1,graphic)
  
    iplot = gsn_csm_contour(wks,prob_u,res_pro)
    overlay(plot, iplot)
    
    draw(plot)
    frame(wks)


end
    
  







    
   