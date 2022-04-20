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
  var_name = (/"omega.mon.mean.nc"/)

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


  q = f1->omega(timeind,{1000:100},{lat1:lat2},{lon1:lon2})
  q= q*-50

  plev = q&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

  q_col =  dim_avg_wgt_n_Wrap(q, dz, 1, 1)
  q_col = (/rmMonAnnCycTLL(q_col)/)
  ; q_col = (/dtrend_msg_n(ispan(1,dimsizes(q_col&time),1),q_col,False,True,0)/)
  ; q_col_nep = dim_avg_n_Wrap(q_col(:,{0:20},{210:280}), (/1,2/))
  p_chazhi = q_col
  
print("=====finish read data==========")



;;;;;;split to 12 month  is a array size (12,yearnum)
  
     p_chazhi_12 = split_12(p_chazhi)
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

  ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0) , 0, 0)
  ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{-30:30},{40:110}), area({-30:30}), 1.0, 0) , 0, 0)
  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:20},{250:280}), area({0:20}), 1.0, 0) , 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum


  print("*******ssta_index_rc="+ssta_index_rc)


  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
 
  ;;;;;要分析的变量
  p_annual = month_to_season(p_chazhi(:,{lat1:lat2},{lon1:lon2}), season)
  p_annual = (/dim_avg_n(p_chazhi_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)
  
  p_annual_rc = regCoef_n(ispan(1,yearnum,1), p_annual , 0, 0)
  p_annual_rc = p_annual_rc*yearnum
  copy_VarCoords(p_annual(1,:,:), p_annual_rc)
   
  prob_obs = rc_test(p_annual_rc)


  ;;;;




  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
    ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{-30:30},{40:110}), area({-30:30}), 1.0, 0)
  ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:20},{250:280}), area({0:20}), 1.0, 0)  
  ; ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{0:20},{250:280}), (/1,2/))
  ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{-5:5},{190:240}), (/1,2/))
  ; ssta_local_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
  ; ssta_local_index = ssta_index

  ;;;;去趋势
  p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)

  ;;;correlation with sst
    ; r_sst_p := regCoef_n(ssta_index, p_annual, 0, 0)
    ; copy_VarCoords(p_annual(0,:,:), r_sst_p)
    ; r_sst_p = (/r_sst_p*ssta_index_rc/)
    r_sst_p := partial_reg(ssta_index,p_annual,ssta_local_index)
    r_sst_p = (/r_sst_p*ssta_index_rc/)
   ;;;;
 
   percent  = avg(r_sst_p({0:20},{250:280}))/avg(p_annual_rc({0:20},{250:280}))



   printMinMax(r_sst_p, 1)
   

    print("*******显著性检验*****************")

    ;;;t test 
    prob_p1= r_sst_p@prob
    copy_VarCoords(p_annual(0,:,:), prob_p1)

    


  print("*********开始画图********")
  ;;;;
    ; wks   = gsn_open_wks("pdf", "NCEP_NCAR omega trend and regression to IND SST")
     ; wks   = gsn_open_wks("pdf", "NCEP_NCAR omega trend and regression to ATL SST")
    wks   = gsn_open_wks("pdf", "NCEP_NCAR omega trend and regression to local SST")
  
    plot=new(2,graphic)
    plot1 = plot

    str = ispan(1,2,1)+ ""
    str2 = ispan(3,4,1)+ ""
 
    iplot = new(2,graphic)

;;;;;


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
  ; res@lbLabelBarOn= False

 
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbLabelStride=1
  res@lbOrientation="horizontal"
  res@mpLimitMode           = "LatLon"

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

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.5
  res@tmYLMajorThicknessF=1.5
  res@tmBorderThicknessF = 1.5

  res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

  res@cnMinLevelValF       =   -1.6
  res@cnMaxLevelValF       =   1.6
  res@cnLevelSpacingF      =  0.4


  res@cnExplicitLabelBarLabelsOn = True
  ; res@lbLabelStrings =  res@cnLevels
  res@vcExplicitLabelBarLabelsOn = True
  res@stExplicitLabelBarLabelsOn = True
  res@lbLabelFontHeightF = 0.016
  res@lbAutoManage = False


  ; res@gsnLeftString = "The components of omega trend linked to ATL "
    ; res@gsnLeftString = "The components of JFMAM shum & 850UV trend linked to TNA"

  ; res@gsnLeftStringOrthogonalPosF  = -0.12

   res@gsnLeftString = "(a) Trends of vertically averaged w"
   res@gsnLeftStringFontHeightF = 0.015

   plot(0) = gsn_csm_contour_map(wks,p_annual_rc, res)
     
   ; res@gsnLeftString = "(b) Trends linearly congruent with ATL warming"+" "+sprintf("%4.2f",percent*100)+"%"
   ; res@gsnLeftString = "(c) Trends linearly congruent with local warming"+" "+sprintf("%4.2f",percent*100)+"%"
   res@gsnLeftString = "(c) Trends linearly congruent with IND"+" "+sprintf("%4.2f",percent*100)+"%"
   plot(1) = gsn_csm_contour_map(wks,r_sst_p, res)
  
  
  ; delete(res@cnLevelSelectionMode)
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillColors)



      ;;;添加indo-pacific warm pool 
      ; ;;添加区域
      ;;warm pool
    do i=0,1
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
     end do 



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

  
    iplot(0) = gsn_csm_contour(wks,prob_obs,res_pro)
    iplot(1) = gsn_csm_contour(wks,prob_p1,res_pro)
   
   
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數
    
    do i=0,1
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
    end do
         
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
    do i=0,1
      overlay(plot(i), iplot(i))
    end do
   
      delete(lat1)
      delete(lat2)
      delete(lon1)
      delete(lon2)


  ; pres = True
  ; pres@gsnMaximize = True
  ; pres@gsnPanelLabelBar = True
  ; pres@gsnPanelBottom = 0.05

  ; gsn_panel(wks, plot, (/2,1/), pres)
draw(plot(1))
frame(wks)
 
 print("done")

end