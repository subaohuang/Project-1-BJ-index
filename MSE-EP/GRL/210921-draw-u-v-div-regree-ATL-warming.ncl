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


function partial_reg(data1,data2,data3) ;;;data1:index data2: varible data3:effect want to remove
begin 
    ;;data1 :x1 data2 :y ; data3: x2
      r12 = escorc_n(data1, data2, 0, 0) ;;;二维 ry1
      r13 = escorc_n(data1, data3, 0, 0) ;;;二维 r12
      r23 = escorc_n(data3, data2, 0, 0) ;;;二维 ry2
      std2 = dim_stddev_n_Wrap(data2, 0) ;;;二维 sy
      std1 = dim_stddev_n_Wrap(data1, 0) ;;;一个数字 s1
      std1_2d = std2
      std1_2d = (/std1/)

      reg12 = (/r12 - r23*r13/)/(1-r13^2)*(std2/std1_2d)
      corr = (/r12 - r23*r13/)/sqrt((1-r23^2)*(1-r13^2))
      copy_VarCoords(data2(1,:,:), reg12)
      ;偏相关的有效自由度
        ; r1_auto = esacr_n(data1,1,0)
        ; r2_auto = esacr_n(data2,1,0)
        
        ; one = r2_auto(0,:,:)
        ; one =1
        ; N = dimsizes(data1&time)*(one- r1_auto(0)*r2_auto(0,:,:))/(one+r1_auto(0)*r2_auto(0,:,:))
      N = dimsizes(data1)
      ; N = equiv_sample_size(data1,0.05,0)

      prob_mam = rtest(corr, N-1, 0)
      copy_VarCoords(data2(0,:,:), prob_mam)
      reg12@corr = corr
      reg12@prob = prob_mam

return reg12
end 

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 202012

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

    sst=lonFlip(f1->sst(timeind,:,:))
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

  

; ;;;read GPCP precipitation data
    
;     f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")   
;     date:=cd_calendar(f2->time, 1)
;     timeind:=ind(date.ge.startmon.and.date.le.endmon)
;     p=f2->precip(timeind,:,:)
;     printVarSummary(p)

;     p_chazhi = p
;     p_chazhi=(/rmMonAnnCycTLL(p_chazhi)/)

; ;;;; 

  path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
  var_name = (/"uwnd.mon.mean.nc","vwnd.mon.mean.nc","shum.mon.mean.nc"/)

  filein = path_wind + var_name 

  ;;;region 
  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360
  ;;;

  ; ;;;global region 
  ; lat1 = -50
  ; lat2 = 50
  ; lon1 = 0
  ; lon2 = 360
  ; ;;;


  startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
  endyear_r = str_split_by_length(endmon, 4 )

  startyear = stringtoint(startyear_r(0))
  endyear = stringtoint(endyear_r(0))

  yearnum = endyear - startyear +1


  f1 = addfile(filein(0),"r")

  date:=cd_calendar(f1->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)


  u = f1->uwnd(timeind,{1000:850},{lat1:lat2},{lon1:lon2})


  f2 = addfile(filein(1),"r")
  date:=cd_calendar(f2->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  v = f2->vwnd(timeind,{1000:850},{lat1:lat2},{lon1:lon2})

  dv = uv2dv_cfd(u,v,u&lat,v&lon, 3)  
  copy_VarCoords(v, dv)
  dv = (/dv*-1*10^6/)
  p_chazhi = dim_avg_n_Wrap(dv, 1) 
  
  u_chazhi = u(:,{1000},:,:)
  v_chazhi = v(:,{1000},:,:)

  u_chazhi = dim_avg_n_Wrap(u(:,{1000:850},:,:),1)
  v_chazhi = dim_avg_n_Wrap(v(:,{1000:850},:,:),1)



; print("=====finish read data==========")

;;;;;;split to 12 month  is a array size (12,yearnum)
  
     p_chazhi_12 = split_12(p_chazhi)
     u_chazhi_12 = split_12(u_chazhi)
     v_chazhi_12 = split_12(v_chazhi)
     ssta_12 = split_12(ssta)
   
;;;;;



   season = "MAM"

  ;;;region
   lat1 = -50
   lat2 = 50
   lon1 = 0
   lon2 = 360
  ;;;
 
 
  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season)
  
  ssta_annual = (/dim_avg_n(ssta_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)

  ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:30},{280:360}), area({0:30}), 1.0, 0) , 0, 0)
  ; ssta_index_rc  =ssta_index_rc*yearnum

  ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:20},{250:280}), area({0:20}), 1.0, 0) , 0, 0)
  ; ssta_index_rc  =ssta_index_rc*yearnum

 
  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:30},{40:110}), area({0:30}), 1.0, 0) , 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum

  print("*******ssta_index_rc="+ssta_index_rc)


  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
 
  ;;;;;要分析的变量
  p_annual = month_to_season(p_chazhi(:,{lat1:lat2},{lon1:lon2}), season)
  p_annual = (/dim_avg_n(p_chazhi_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)


  ;;;;;要分析的变量
  u_annual = month_to_season(u_chazhi(:,{lat1:lat2},{lon1:lon2}), season)
  u_annual = (/dim_avg_n(u_chazhi_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)



  ;;;;;要分析的变量
  v_annual = month_to_season(v_chazhi(:,{lat1:lat2},{lon1:lon2}),season)
  v_annual = (/dim_avg_n(v_chazhi_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)
  ;;;;



  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:30},{280:360}), area({0:30}), 1.0, 0)
  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:20},{250:280}), area({0:20}), 1.0, 0)  
   ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:30},{40:110}), area({0:30}), 1.0, 0)  
  ; ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{0:20},{250:280}), (/1,2/))
  ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{-5:5},{190:240}), (/1,2/))
  ; ssta_local_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
  ; ssta_local_index = ssta_index

  ;;;;去趋势

  p_annual_rc = regCoef_n(ispan(1,dimsizes(p_annual&time),1), p_annual , 0, 0)
  copy_VarCoords(p_annual(1,:,:), p_annual_rc)
  p_annual_rc = (/p_annual_rc*yearnum/)

  p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)
  u_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),u_annual,False,True,0)/)
  v_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),v_annual,False,True,0)/)
  

  ;;;correlation with sst
    ; r_sst_p := regCoef_n(ssta_index, p_annual, 0, 0)
    ; copy_VarCoords(p_annual(0,:,:), r_sst_p)
    ; r_sst_p = (/r_sst_p*ssta_index_rc/)
    r_sst_p := partial_reg(ssta_index,p_annual,ssta_local_index)
    r_sst_p = (/r_sst_p*ssta_index_rc/)
   ;;;

    ;;;correlation with sst
    r_sst_u :=partial_reg(ssta_index,u_annual,ssta_local_index)
    r_sst_v :=partial_reg(ssta_index,v_annual,ssta_local_index)
    r_sst_u = (/r_sst_u*ssta_index_rc/)
    r_sst_v = (/r_sst_v*ssta_index_rc/)
    
    ;;;;;


   printMinMax(r_sst_p, 1)
   
  


    print("*******显著性检验*****************")

    ;;;t test 
    prob_p1= r_sst_p@prob
    prob_p2= r_sst_u@prob
    prob_p3= r_sst_v@prob
    copy_VarCoords(p_annual(0,:,:), prob_p1)
    copy_VarCoords(u_annual(0,:,:), prob_p2)
    copy_VarCoords(v_annual(0,:,:), prob_p3)
   
    ; prob_p1 = where(p_percent.ge.0.9, 0.01, 2)
    ; ; prob_p2 = where(u_percent.ge.0.5, 0.01, 2)
    ; ; prob_p3 = where(v_percent.ge.0.5, 0.01, 2)
    ; prob_p4 = where(olr_percent.ge.0.9, 0.01, 2)


    ; prob_p1 = where(p_percent.le.-0.9, 0.06, prob_p1)
    ; ; prob_p2 = where(u_percent.le.-0.5, 0.06, 2)
    ; ; prob_p3 = where(v_percent.le.-0.5, 0.06, 2)
    ; prob_p4 = where(olr_percent.le.-0.9, 0.06, prob_p4)



    ; prob_p1 = where(prob_p1.eq.0, 2, prob_p1)
    ;;;;


    r_sst_u = where(prob_p2.le.0.1.or.prob_p3.le.0.1, r_sst_u, r_sst_u@_FillValue)
;;;;


   percent = avg(r_sst_p({0:20},{250:280}))/avg(p_annual_rc({0:20},{250:280}))


 print("percent="+percent)


;;;;
    ; wks   = gsn_open_wks("pdf", "CMAP&olr回归到印度洋SST的回归系数")
    ; wks   = gsn_open_wks("pdf", "NCEP_NCAR shun&UV回归大西洋SST的回归系数-去掉ENSO的影响")
    ; wks   = gsn_open_wks("pdf", "NCEP_NCAR shun&UV回归到局地SST的回归系数-去掉ENSO的影响")
    ; wks   = gsn_open_wks("pdf", "NCEP_NCAR shun&UV&div regreed into ATL SST-rmENSO")
    ; wks   = gsn_open_wks("pdf", "NCEP_NCAR shun&UV&div regreed into local SST-rmENSO")
    wks   = gsn_open_wks("pdf", "NCEP_NCAR shun&UV&div regreed into IND-rmENSO")
    ; LeftString =  "JFMAM 1000-850 mean UV & DIV trends linked to NTA warming "
    ; LeftString =  "JFMAM 1000-850 mean UV & DIV trends linked to local warming "
    LeftString =  "JFMAM 1000-850 mean UV & DIV trends linked to IND warming "
    ; wks   = gsn_open_wks("pdf", "NCEP_NCAR shun&UV回归大西洋SST的回归系数-去掉local的影响")
  
    plot=new(4,graphic)
    plot1 = plot

    str = ispan(1,4,1)+ ""
    str2 = ispan(5,8,1)+ ""
 
    iplot = new(4,graphic)

;;;;;



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

    res@lbLabelBarOn= True


    res@pmLabelBarOrthogonalPosF = 0.2

    res@lbLabelStride=1
    res@lbOrientation="horizontal"
    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = lat1
    ; res@mpMaxLatF             = lat2
    ; res@mpMinLonF             = lon1
    ; res@mpMaxLonF             = lon2
    ; res@mpMinLatF             = -5
    ; res@mpMaxLatF             = 30
    ; res@mpMinLonF             = 100
    ; res@mpMaxLonF             = 170

    ; res@mpMinLatF             = -30
    ; res@mpMaxLatF             = 40
    ; res@mpMinLonF             = 200
    ; res@mpMaxLonF             = 360

    
    res@mpMinLatF             = -50
    res@mpMaxLatF             = 50
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
    res@mpFillOn = True


    ; res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "CBR_drywet"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels


    ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
    ; res@cnLevels := (/-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4/)*0.5
    res@cnLevels := (/-0.25,-0.2,-0.15,-0.1,0.1,0.15,0.2,0.25/)*3.0

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

   
    ; res@vcGlyphStyle="CurlyVector"
    res@vcGlyphStyle="LineArrow"
    res@vcMonoLineArrowColor = True
    res@vcLineArrowColor = "black"
    res@vcMinDistanceF           = 0.01          ; thin vectors

    res@vcRefLengthF             = 0.025    ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcVectorDrawOrder="PostDraw"

    ; resv@vcMonoLineArrowColor = True
    ; resv@vcLineArrowColor = "blue"
    res@vcLineArrowThicknessF=1.2
    ; res@vcLineArrowHeadMaxSizeF=0.01
    ; res@vcLineArrowHeadMinSizeF=0.01
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006

    res@vcRefMagnitudeF          = 0.2*10
    res@vcMinMagnitudeF          =.1
    res@gsnLeftStringFontHeightF = 0.015
    ; res@gsnLeftStringFontHeightF = 0.02

    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "2.0 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.2             ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013

    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False

    res@gsnLeftString = LeftString
    ; res@gsnLeftString = "The components of JFMAM shum & 850UV trend linked to TNA"

    res@gsnLeftStringOrthogonalPosF  = -0.12

  i = 0
  plot(i) = gsn_csm_vector_scalar_map(wks,r_sst_u,r_sst_v , r_sst_p, res)

  
  ; delete(res@cnLevelSelectionMode)
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillColors)


   
  ;;;添加indo-pacific warm pool 
    ; ;;添加区域
   ;;warm pool
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF= 2.0
    pgres@gsLineDashPattern= 1

    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=1
    str3=unique_string("dum3")
    str4=unique_string("dum4")
 
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 250
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)

    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   ;;;
   pgres@gsLineColor="red"

   lat1 = 10
   lat2 = 30
   lon1 = 280
   lon2 = 330

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   
    ; plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   
   


    ;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(i) = gsn_csm_contour(wks,prob_p1,res_pro)

   
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
    
         
        ; res_pro=True
        ; res_pro@gsnDraw=False
        ; res_pro@gsnFrame=False
        ; res_pro@cnFillOn = True

        ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
        ; res_pro@cnLevels = (/0,0.05,1.0/)
        ; res_pro@lbLabelBarOn = False

        ; res_pro@cnLinesOn=False
        ; res_pro@cnFillDrawOrder = "PostDraw"

        ; res_pro@cnMonoFillPattern = True
        ; res_pro@cnFillPattern = 17
        ; res_pro@cnFillDotSizeF = 0.0028
        ; res_pro@cnMonoFillColor = False
        ; res_pro@cnFillColors = (/"Transparent","red","blue","Transparent"/)
        ; ; res_pro@cnLevels = ispan(-6,6,1)
        ; res_pro@cnLineLabelsOn=False
        ; res_pro@tmYRMode = "Automatic"
        ; res_pro@cnInfoLabelOn=False

        ; iplot(i*2) = gsn_csm_contour(wks,prob_p1,res_pro)

        ; iplot(i*2+1) = gsn_csm_contour(wks,prob_p4,res_pro)
    
    

    ;;;  
    
      overlay(plot(i), iplot(i))
      
   
      delete(lat1)
      delete(lat2)
      delete(lon1)
      delete(lon2)


  pres = True
  pres@gsnMaximize = True
  
  pres@gsnPanelBottom = 0.05

  ; gsn_panel(wks, plot, (/2,2/), pres)
  draw(plot(0))
  frame(wks)
  
; system("convert -geometry 1800x1800 -density 600x600 -trim  GPCP&olr回归到印度洋SST的回归系数.png  GPCP&olr回归到印度洋SST的回归系数_high.png")     

 
 print("done")

end