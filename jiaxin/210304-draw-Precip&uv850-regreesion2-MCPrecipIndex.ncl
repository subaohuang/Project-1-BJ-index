;;********************************************************
;;计算SST倾向以及SST值与GPCP降水的趋势
;;其中热通量的方向都是向上为正，向下为负
;;的相关系数场，
;;NCEP-NCAR的精度为1.87×1.9
;;
;;                                         lnshuheng 2020.02.18
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

undef("remove_enso_signal")
function remove_enso_signal(nino3_4,  data)
  begin  
   ;;;;;
      ; data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_de  = data
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


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 201612

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
  order  = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/)
 

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,::-1,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ;;;
    delete(date)
    delete(timeind)

    ssta_de = ssta
    ; ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_in = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
    printVarSummary(ssta_in)

    nino3_4 = month_to_season(ssta_in, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)


    ; nino3_4_new = nino3_4                ;;;;前期ENSO
    ; nino3_4_new(0:33) = (/nino3_4(1:)/)
    ; nino3_4_new(34) =ssta_index(35*12-1) 

    nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
;;;

;;;area mean 
      rad    = 4.0*atan(1.0)/180.0
      re     = 6371220.0
      rr     = re*rad

      dlon   = abs(lon(2)-lon(1))*rr
      dx     = dlon*cos(lat*rad)
      ;          
      jlat  = dimsizes( lat )
      dy     = new ( jlat, typeof(dx))      ; close enough
      
      dy(0)  = abs(lat(2)-lat(1))*rr
      dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5   
      dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr

      area   = dx*dy                                ; cell area function of latitude only
      copy_VarCoords(lat, area)
      printVarSummary(area)
;;;;

  

;;read GPCP precipitation data
    
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")   
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p=f2->precip(timeind,:,:)
    printVarSummary(p)

    p = (/rmMonAnnCycTLL(p)/)
;;;

;;read era data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"uwnd","vwnd"/)
    filein = path_wind + var_name +".mon.mean.nc"

    ;;;region 
    lat1 = -90
    lat2 = 90
    lon1 = 0
    lon2 = 360
    ;;;

    f1 = addfile(filein(0),"r")

    date :=cd_calendar(f1->time, 1)

    timeind :=ind(date.ge.197901.and.date.le.201612)

    u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    v = short2flt(f2->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    u!1 = "lat"
    u!2 = "lon"

    v!1 = "lat"
    v!2 = "lon"
    printVarSummary(u)


    u = (/rmMonAnnCycTLL(u)/)
    v = (/rmMonAnnCycTLL(v)/)

;;;;


;;;;
    ; wks   = gsn_open_wks("pdf", "CMAP&olr回归到印度洋SST的回归系数")
    wks   = gsn_open_wks("pdf", "Precip&850uv regression on GPCP MC precip")
    wks1   = gsn_open_wks("pdf", "HadISST&850uv regression on GPCP MC precip")
  
    plot=new(8,graphic)
    plot1 = plot

  
 
    iplot = new(8,graphic)

    season =  (/"MAM","JJA","SON","DJF"/)
;;;;;



do i =0,3

  ;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
  ;;;

 
  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season(i))
  printVarSummary(ssta_annual)

  ;;;remove_enso_signal 
  ssta_annual = (/remove_enso_signal(nino3_4, ssta_annual)/)
  ;;;

  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
 
  ;;;;;要分析的变量
  p_annual = month_to_season(p(:,{lat1:lat2},{lon1:lon2}), season(i))
  p_MAM = month_to_season(p(:,{lat1:lat2},{lon1:lon2}), "MAM")


  ;;;remove_enso_signal 
  p_annual = (/remove_enso_signal(nino3_4, p_annual)/)
  p_MAM = (/remove_enso_signal(nino3_4, p_MAM)/)
  ;;;;;;;


 
  ;;;;;要分析的变量
  u_annual = month_to_season(u(:,{lat1:lat2},{lon1:lon2}), season(i))
  ;;;remove_enso_signal 
  u_annual = (/remove_enso_signal(nino3_4, u_annual)/)
  ;;;
  
  ;;;;;要分析的变量
  v_annual = month_to_season(v(:,{lat1:lat2},{lon1:lon2}), season(i))

  ;;;remove_enso_signal 
  v_annual = (/remove_enso_signal(nino3_4, v_annual)/)
  ;;;



  ;;;;去趋势
  p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)
  u_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),u_annual,False,True,0)/)
  v_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),v_annual,False,True,0)/)
  
  
  ; MCPRE_index = dim_avg_n_Wrap( p_annual(:,{-10:10},{90:150}), (/1,2/))
  MCPRE_index = dim_avg_n_Wrap( p_MAM(:,{-10:10},{90:150}), (/1,2/))


  ;;;Precip 
    r_pre := regCoef_n(MCPRE_index, p_annual, 0, 0)
    copy_VarCoords(p_annual(0,:,:), r_pre)
  ;;;;;

  ;;;sst
    r_sst := regCoef_n(MCPRE_index, ssta_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst)

    ;;;correlation with u
    r_u := regCoef_n(MCPRE_index, u_annual, 0, 0)
    copy_VarCoords(u_annual(0,:,:), r_u)

    ;;;;;

    ;;;correlation with v
    r_v := regCoef_n(MCPRE_index, v_annual, 0, 0)
    copy_VarCoords(v_annual(0,:,:), r_v)

    ;;;;;

  

    print("*******显著性检验*****************")

    ;;;t test 
    prob_p1= rc_test(r_pre)
    prob_p2= rc_test(r_u)
    prob_p3= rc_test(r_v)
    prob_p4= rc_test(r_sst)

    ; r_u = where(prob_p2.le.0.05.or.prob_p3.le.0.05, r_u, r_u@_FillValue)
  
  print("*********开始画图********")

;;**************************************************************
;画图
;;***************************************************************


  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  ; res@cnFillPalette="MPL_BrBG"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn= False

  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbLabelStride=1
  res@lbOrientation="horizontal"

   
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmBorderThicknessF=1.0

  res@tmXBLabelFontHeightF=0.02
  res@tmYLLabelFontHeightF=0.02
  res@tmXBMajorThicknessF=1.5
  res@tmYLMajorThicknessF=1.5
  res@tmBorderThicknessF = 1.5

  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 
  
  ; res@tmYLValues = (/-30,-20,-10,0,10,20,30/) 
  ; res@tmYLLabels = (/"30S","20S","10S","0","10N","20N","30N"/) 

  ; res@tmYLMinorValues = ispan(-30,30,5)

  res@gsnScalarContour        = True
  res@mpFillOn = False


  res@cnFillPalette = "precip_diff_12lev"
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "CBR_drywet"
  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels


  res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5/)

  res@lbLabelStride=1


  res@cnExplicitLabelBarLabelsOn = True
  ; res@lbLabelStrings =  res@cnLevels
  res@vcExplicitLabelBarLabelsOn = True
  res@stExplicitLabelBarLabelsOn = True
  res@lbLabelFontHeightF = 0.016
  res@lbAutoManage = False

  
  ; res@cnFillPalette = "BlRe"
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

  ; res@cnMinLevelValF       =   -0.5
  ; res@cnMaxLevelValF       =   0.5
  ; res@cnLevelSpacingF      =  0.1

  res@vcGlyphStyle="LineArrow"
  res@vcMonoLineArrowColor = True
  res@vcLineArrowColor = "black"
  res@vcMinDistanceF           = 0.02          ; thin vectors

  res@vcRefLengthF             = 0.025    ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  ; res@vcVectorDrawOrder="PostDraw"

  
  ; res@vcLineArrowThicknessF=1.2
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  ; res@vcLineArrowHeadMinSizeF = 0.006

  res@vcRefMagnitudeF          = 1.0
  ; res@vcMinMagnitudeF          = 0.01

  ; res@vcRefAnnoString1On = True
  ; res@vcRefAnnoString1 = "2 m/s"
  ; res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
  ; res@vcRefAnnoParallelPosF = 1.02
  ; res@vcRefAnnoFontHeightF = 0.013

  ; res@vcRefAnnoBackgroundColor = "Transparent"
  ; res@vcRefAnnoPerimOn = False

  res@gsnLeftString = order(i)+" "+ season(i)+ " linear regression  Precip & 850UV trend on MCPRE_index"
  res@gsnLeftStringFontHeightF = 0.0185

  ; res@gsnLeftStringOrthogonalPosF  = -0.12

  
  plot(i) = gsn_csm_vector_scalar_map(wks,r_u, r_v, r_pre, res)
  printVarSummary(r_pre)
  printVarSummary(r_u)
  printVarSummary(r_v)



  
  ; delete(res@cnLevelSelectionMode)
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillColors)
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "MPL_bwr"

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  res@lbLabelStride=1

  res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res@cnMinLevelValF      = -1.0       ; set min contour level
  res@cnMaxLevelValF      = 1.0     ; set max contour level
  res@cnLevelSpacingF     = 0.1       ; set contour spacing

  res@gsnLeftStringOrthogonalPosF  = 0.06


  res@gsnLeftString = order(i+4)+" "+ season(i)+ " linear regression  SST  on MCPRE_index"
  plot(i+4) = gsn_csm_contour_map(wks1, r_sst, res)   

  ; plot(i) = gsn_csm_vector_scalar_map(wks,r_sst_u,r_sst_v , r_sst_w, res)
  

    ;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(i) = gsn_csm_contour(wks,prob_p1,res_pro)

    iplot(i+4) = gsn_csm_contour(wks1,prob_p4,res_pro)
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    iplot(i+4) = gsn_contour_shade(iplot(i+4) , 0,0.1,opt)
    
  
    overlay(plot(i), iplot(i))
    overlay(plot(i+4), iplot(i+4))
 
    delete(lat1)
    delete(lat2)
    delete(lon1)
    delete(lon2)
end do

  pres = True
  pres@gsnMaximize = True
  
  pres@gsnPanelBottom = 0.05

  gsn_panel(wks, plot(0:3), (/2,2/), pres)

  gsn_panel(wks1, plot(4:7), (/2,2/), pres)
  
; system("convert -geometry 1800x1800 -density 600x600 -trim  GPCP&olr回归到印度洋SST的回归系数.png  GPCP&olr回归到印度洋SST的回归系数_high.png")     

 
 print("done")

end