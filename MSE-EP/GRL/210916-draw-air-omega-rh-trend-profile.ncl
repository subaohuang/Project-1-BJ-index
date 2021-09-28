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

    data_12 = new((/12,yearnum,dim(1)/),float)
    do i=0,11   
    data_12(i,:,:) = data(i::12,:)
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

  
  path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
  var_name = (/"rhum.mon.mean.nc","omega.mon.mean.nc","air.mon.mean.nc"/)

  filein = path_wind + var_name 

  ;;;region 
  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360
  ;;;

  ;;;;
  startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
  endyear_r = str_split_by_length(endmon, 4 )

  startyear = stringtoint(startyear_r(0))
  endyear = stringtoint(endyear_r(0))

  yearnum = endyear - startyear +1


  f1 = addfile(filein(0),"r")
  date:=cd_calendar(f1->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  
  rh = f1->rhum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
  rh = rh/100.0

  f2 = addfile(filein(1),"r")
  date:=cd_calendar(f2->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  w = f2->omega(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
  w = w *-50

  f3 = addfile(filein(2),"r")
  date:=cd_calendar(f3->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  t = f3->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
  
  
;;;;;
  print("=====finish read data==========")
;;;


;;; 
   w_nep = dim_avg_n_Wrap(w(:,:,{0:20},{250:280}), (/2,3/))
   rh_nep = dim_avg_n_Wrap(rh(:,:,{0:20},{250:280}), (/2,3/))
   t_nep = dim_avg_n_Wrap(t(:,:,{0:20},{250:280}), (/2,3/))
;;;;

  
;;;;;;split to 12 month  is a array size (12,yearnum)
   w_nep_12 = split_12(w_nep)
   rh_nep_12 = split_12(rh_nep)
   t_nep_12 = split_12(t_nep)
;;;;;

;;;annual 
  w_annual =dim_avg_n(w_nep_12(0:4,:,:),0)
  copy_VarCoords(w_nep_12(1,:,:), w_annual)

  rh_annual = dim_avg_n(rh_nep_12(0:4,:,:),0)
  copy_VarCoords(rh_nep_12(1,:,:), rh_annual)

  t_annual = dim_avg_n(t_nep_12(0:4,:,:),0)
  copy_VarCoords(t_nep_12(1,:,:), t_annual)
;;;;


;;;;trend 
   w_rc = regCoef_n(ispan(1,yearnum,1),w_annual, 0, 0)
   w_rc = w_rc*yearnum
   copy_VarCoords(w_annual(1,:), w_rc)

   rh_rc = regCoef_n(ispan(1,yearnum,1),rh_annual, 0, 0)
   rh_rc = rh_rc*yearnum
   copy_VarCoords(rh_annual(1,:), rh_rc)

   t_rc = regCoef_n(ispan(1,yearnum,1),t_annual, 0, 0)
   t_rc = t_rc*yearnum
   copy_VarCoords(t_annual(1,:), t_rc)


;;;

  plev = w&level


  plev = w&level
  ptop = min(plev)
  psfc = 1013.0
  dz = dpres_plevel(plev, psfc, ptop, 0)

  rh_rc_col =  dim_avg_wgt_n_Wrap(rh_rc, dz, 1, 0)


  ; q_col = (/rmMonAnnCycTLL(q_col)/)
  ; q_col = (/dtrend_msg_n(ispan(1,dimsizes(q_col&time),1),q_col,False,True,0)/)
  ; q_col_nep = dim_avg_n_Wrap(q_col(:,{0:20},{210:280}), (/1,2/))
  ; p_chazhi = q_col





;   ssta_annual = (/dim_avg_n(ssta_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)
; ;;;contruct sst index 
;   season = "MAM"
;   ;;;region
;    lat1 = -50
;    lat2 = 50
;    lon1 = 0
;    lon2 = 360
;   ;;;
 
;   ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season)
  
;   ssta_annual = (/dim_avg_n(ssta_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)

;   ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0) , 0, 0)
;   ssta_index_rc  =ssta_index_rc*yearnum
;   print("*******ssta_index_rc="+ssta_index_rc)
;   ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
 
;   ;;;;;要分析的变量
;   p_annual = month_to_season(p_chazhi(:,{lat1:lat2},{lon1:lon2}), season)
;   p_annual = (/dim_avg_n(p_chazhi_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)
  
;   p_annual_rc = regCoef_n(ispan(1,yearnum,1), p_annual , 0, 0)
;   p_annual_rc = p_annual_rc*yearnum
;   copy_VarCoords(p_annual(1,:,:), p_annual_rc)
   
;   prob_obs = rc_test(p_annual_rc)
; ;;;;
  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
  ; ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:20},{250:280}), area({0:20}), 1.0, 0)  
  ; ; ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{0:20},{250:280}), (/1,2/))
  ; ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{-5:5},{190:240}), (/1,2/))
  ; ; ssta_local_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
  ; ; ssta_local_index = ssta_index

  ; ;;;;去趋势
  ; p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)

  ; ;;;correlation with sst
  ;   ; r_sst_p := regCoef_n(ssta_index, p_annual, 0, 0)
  ;   ; copy_VarCoords(p_annual(0,:,:), r_sst_p)
  ;   ; r_sst_p = (/r_sst_p*ssta_index_rc/)
  ;   r_sst_p := partial_reg(ssta_index,p_annual,ssta_local_index)
  ;   r_sst_p = (/r_sst_p*ssta_index_rc/)
  ;  ;;;;
 
  ;  percent  = avg(r_sst_p({0:20},{250:280}))/avg(p_annual_rc({0:20},{250:280}))
  ;  printMinMax(r_sst_p, 1)
  ;   print("*******显著性检验*****************")

  ;   ;;;t test 
  ;   prob_p1= r_sst_p@prob
  ;   copy_VarCoords(p_annual(0,:,:), prob_p1)
;;;;;


;;;;

  print("*********开始画图********")
  ;;;;
    wks   = gsn_open_wks("pdf", "w_rc_profile")
    wks1  = gsn_open_wks("pdf", "relativehum_rc_profile")
    wks2  = gsn_open_wks("pdf", "t_rc_profile")
  
    plot=new(1,graphic)
    plot1 = plot
    plot2 = plot

    str = ispan(1,2,1)+ ""
    str2 = ispan(3,4,1)+ ""
 
    iplot = new(2,graphic)

;;;;;


;;**************************************************************
;画图
;;***************************************************************
;;;;
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False

    resxy@xyLineThicknessF  =  (/2.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True

    resxy@xyLineColors       :=  (/"black"/)          ; change line color
    ; resxy@xylineColors :=  (/"red","blue","green","black"/)
    resxy@xyMarkLineModes    :=  (/"MarkLines"/)

    resxy@xyDashPatterns := (/0/)
    resxy@xyMonoMarker        = False
    resxy@xyMarkers           := (/16/)
    resxy@xyMarkerSizes = (/0.01/)
    ; resxy@xyMarkerColors :=  (/"blue","blue","red","red"/)
    resxy@xyMarkerColors :=  (/"black"/)

    resxy@gsnLeftStringFontHeightF = 0.02
    resxy@gsnLeftString = "Trends of vertical velocity"

    resxy@tmXBLabelFontHeightF=0.025
    resxy@tmYLLabelFontHeightF=0.025
    resxy@tmXBMajorThicknessF=1.0
    resxy@tmYLMajorThicknessF=1.0

    resxy@tmXBMode  =  "Explicit"
    resxy@trXMaxF =  1.4
    resxy@trXMinF = 0.
    resxy@tmXBValues = ispan(0,14,2)*0.1
    resxy@tmXBLabels = ispan(0,14,2)*0.1
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues = ispan(0,14,1)*0.1
    resxy@tiYAxisString = ""
   
    resxy@trYReverse        = True                       ; reverse Y-axis
    resxy@gsnXRefLine = 0.0

    plot = gsn_csm_xy(wks, w_rc ,plev ,resxy)
    

    resxy@gsnLeftString = "Trends of relative humidity"
    resxy@trXMaxF =  0.1
    resxy@trXMinF = -0.1
    resxy@tmXBValues := ispan(-10,10,4)*0.01
    resxy@tmXBLabels := ispan(-10,10,4)*0.01
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues := ispan(-10,10,2)*0.01

    plot1 = gsn_csm_xy(wks1, rh_rc ,plev ,resxy)

    resxy@gsnLeftString = "Trends of temperature"

    resxy@trXMaxF =  1.2
    resxy@trXMinF = 0
    resxy@tmXBValues := ispan(0,12,2)*0.1
    resxy@tmXBLabels := ispan(0,12,2)*0.1
    resxy@tmXBMinorOn = True
    resxy@tmXBMinorValues := ispan(0,12,1)*0.1

    plot2 = gsn_csm_xy(wks2, t_rc ,plev ,resxy)

    
    ;;添加legend 
    gres = True
    gres@YPosPercent = 20.    ; expressed as %, 0->100, sets position of top border of legend 
          ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 55    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = resxy@xyLineColors ; colors chosen
    lineres@lgLineThicknesses = 1.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.
    lineres@lgDashIndexes   = resxy@xyDashPatterns
    lineres@lgMonoMarkerIndex = False 
    lineres@lgMarkerIndexes = resxy@xyMarkers
    lineres@lgItemTypes =  resxy@xyMarkLineModes
    lineres@lgMarkerIndex = resxy@xyMarkers
                        ; expressed as %, 0->100, length of line
    textres = True
    textres@lgLabels = (/"Vertical velocity (W)"/)
    textres@lgLabelFontHeightF = 0.02
    plot = simple_legend(wks,plot,gres,lineres,textres)
    textres@lgLabels = (/"Relative humidity (Rh)"/)
    plot1 = simple_legend(wks1,plot1,gres,lineres,textres)
      textres@lgLabels = (/"Temperature"/)
    plot2 = simple_legend(wks2,plot2,gres,lineres,textres)

draw(plot)
draw(plot1)
frame(wks)
frame(wks1)
draw(plot2)
frame(wks2)

end 
