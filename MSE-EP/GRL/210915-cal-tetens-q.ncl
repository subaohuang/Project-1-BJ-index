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
 
  path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
  var_name = (/"air.mon.mean.nc","shum.mon.mean.nc","rhum.mon.mean.nc"/)

  filein = path_wind + var_name 

  ;;;region 
  lat1 = -50
  lat2 = 50
  lon1 = 0
  lon2 = 360
  ;;;

  startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
  endyear_r = str_split_by_length(endmon, 4 )

  startyear = stringtoint(startyear_r(0))
  endyear = stringtoint(endyear_r(0))

  yearnum = endyear - startyear +1

  f1 = addfile(filein(0),"r")

  date:=cd_calendar(f1->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  t = f1->air(timeind,{1000:300},{lat1:lat2},{lon1:lon2}) ;;; ℃

  f2 = addfile(filein(1),"r")
  date:=cd_calendar(f2->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  q = f2->shum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
  q = q/1000.0

  f3 = addfile(filein(2),"r")
  date:=cd_calendar(f3->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  rh = f3->rhum(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
  
  rh = rh/100.0

   rh_clm = avg(rh)
   t_clm = avg(t)


;;;;;;;;;;;;
    e = t
    printVarSummary(rh)
    printVarSummary(t)

    e = (/6.112*exp(17.67*t/(t+243.5))*rh/)

    plev = q&level
    plev_3d  = conform_dims(dimsizes(q), plev, 1)

    qm = q
    qm = (/0.622*e/plev_3d/)

   
    e := (/6.112*exp(17.67*t/(t+243.5))*rh_clm/)

    qm_fixrh = q
    qm_fixrh = (/0.622*e/plev_3d/)

    
    e := (/6.112*exp(17.67*t_clm/(t_clm+243.5))*rh/)

    qm_fixt = q
    qm_fixt = (/0.622*e/plev_3d/)
;;;


;;;;;Vertical mean 

    pver = q&level
    ; pver = pver*100.0
    pver@units = "Pa"

    dp = dpres_plevel(pver,1000.0,min(pver),0)              ; caculate dp
    ; dp_cfm = conform_dims(dimsizes(q),dp,1)          ; confirm dp to caculate sum(var*dp)

    q_col = dim_avg_wgt_n_Wrap(q, dp, 1, 1)               ; sum(var*dp)
    copy_VarCoords(q(:,1,:,:), q_col)

    t_col = dim_avg_wgt_n_Wrap(t, dp, 1, 1)               ; sum(var*dp)
    copy_VarCoords(t(:,1,:,:), t_col)

    rh_col = dim_avg_wgt_n_Wrap(rh, dp, 1, 1)               ; sum(var*dp)
    copy_VarCoords(rh(:,1,:,:), rh_col)


     ; q_col = q(:,{850},:,:)
    q_col_12 = split_12(q_col)
   
    ; t_col = t(:,{850},:,:)
    t_col_12 = split_12(t_col)
   ; rh_col = rh(:,{850},:,:)
    rh_col_12 = split_12(rh_col)


    qm_col = dim_avg_wgt_n_Wrap(qm, dp, 1, 1)               ; sum(var*dp)
    copy_VarCoords(qm(:,1,:,:), qm_col)

    qm_fixrh_col = dim_avg_wgt_n_Wrap(qm_fixrh, dp, 1, 1)               ; sum(var*dp)
    copy_VarCoords(qm_fixrh(:,1,:,:), qm_fixrh_col)

    qm_fixt_col = dim_avg_wgt_n_Wrap(qm_fixt, dp, 1, 1)               ; sum(var*dp)
    copy_VarCoords(qm_fixt(:,1,:,:), qm_fixt_col)
    




;;;annual 
   ; qm_col = qm(:,{850},:,:)

   q_col_annual = month_to_annual(q_col, 1)
   q_col_annual = (/dim_avg_n_Wrap(q_col_12(0:4,:,:,:), 0)/)

   t_col_annual = month_to_annual(t_col, 1)
   t_col_annual = (/dim_avg_n_Wrap(t_col_12(0:4,:,:,:), 0)/)

   rh_col_annual = month_to_annual(rh_col, 1)
   rh_col_annual = (/dim_avg_n_Wrap(rh_col_12(0:4,:,:,:), 0)/)



   qm_col_12 = split_12(qm_col)
   
   qm_col_annual = month_to_annual(qm_col, 1)
   qm_col_annual = (/dim_avg_n_Wrap(qm_col_12(0:4,:,:,:), 0)/)


   ; qm_fixrh_col = qm_fixrh(:,{850},:,:)
   qm_fixrh_col_12 = split_12(qm_fixrh_col)
   
   qm_fixrh_col_annual = month_to_annual(qm_fixrh_col, 1)
   qm_fixrh_col_annual = (/dim_avg_n_Wrap(qm_fixrh_col_12(0:4,:,:,:), 0)/)


   ; qm_fixt_col = qm_fixt(:,{850},:,:)
   qm_fixt_col_12 = split_12(qm_fixt_col)
   
   qm_fixt_col_annual = month_to_annual(qm_fixt_col, 1)
   qm_fixt_col_annual = (/dim_avg_n_Wrap(qm_fixt_col_12(0:4,:,:,:), 0)/)
;;;;;


;;;nep area mean 

   q_col_annual_nep = dim_avg_n_Wrap(q_col_annual(:,{0:20},{250:280}), (/1,2/))
   t_col_annual_nep = dim_avg_n_Wrap(t_col_annual(:,{0:20},{250:280}), (/1,2/))
   rh_col_annual_nep = dim_avg_n_Wrap(rh_col_annual(:,{0:20},{250:280}), (/1,2/))

  
   qm_col_nep =  dim_avg_n_Wrap(qm_col_annual(:,{0:20},{250:280}), (/1,2/))

   qm_fixrh_col_nep =  dim_avg_n_Wrap(qm_fixrh_col_annual(:,{0:20},{250:280}), (/1,2/))

   qm_fixt_col_nep =  dim_avg_n_Wrap(qm_fixt_col_annual(:,{0:20},{250:280}), (/1,2/))

  
   ; qm_diff = avg(qm_col_nep(20:)) - avg(qm_col_nep(0:19))
   ; qm_fixrh_diff = avg(qm_fixrh_col_nep(20:)) - avg(qm_fixrh_col_nep(0:19))
   ; qm_fixt_diff = avg(qm_fixt_col_nep(20:)) - avg(qm_fixt_col_nep(0:19))
   ; q_diff = avg(q_col_annual_nep(20:)) - avg(q_col_annual_nep(0:19))

    ; print(qm_diff)
   ; print(qm_fixrh_diff)
   ; print(q_diff)
   ; print(qm_fixt_diff)
    
   ; print("percent1="+qm_fixrh_diff/qm_diff)
   ; print("percent2="+qm_fixt_diff/qm_diff)
   
   qm_diff = dim_avg_n_Wrap(qm_col_annual(20:,:,:),0) - dim_avg_n_Wrap(qm_col_annual(0:19,:,:),0)
   copy_VarCoords(qm_col_annual(1,:,:), qm_diff)

   qm_fixrh_diff = dim_avg_n_Wrap(qm_fixrh_col_annual(20:,:,:),0) - dim_avg_n_Wrap(qm_fixrh_col_annual(0:19,:,:),0)
   copy_VarCoords(qm_fixrh_col_annual(1,:,:), qm_fixrh_diff)

   qm_fixt_diff = dim_avg_n_Wrap(qm_fixt_col_annual(20:,:,:),0) - dim_avg_n_Wrap(qm_fixt_col_annual(0:19,:,:),0)
   copy_VarCoords(qm_fixt_col_annual(1,:,:), qm_fixt_diff)

   q_diff =dim_avg_n_Wrap(q_col_annual(20:,:,:),0) - dim_avg_n_Wrap(q_col_annual(0:19,:,:),0)
   copy_VarCoords(q_col_annual(1,:,:), q_diff)
;;;;;


;;; Lin Shuheng' method 
  ; q_2_test =  (/6.112*exp(17.67*(t_clm+t_rc)/(t_clm+t_rc+243.5))*(rh_clm+rh_rc)/)
  ; q_change_test = q_2_test - q_clm_clm
  ; a = 17.67
  ; b = 243.5
  ; t_clm_all  = t_clm

  ; h = a*b/(t_clm_all+b)^2

  ; dq_q = q_diff/avg(q_col_annual_nep(0:19))
  ; ; dq_q = qm_diff/avg(qm_col_nep(0:19))

  ; h_dt = h*(avg(t_col_annual_nep(20:)) - avg(t_col_annual_nep(0:19)))
  
  ; rh_diff = avg(rh_col_annual_nep(20:)) - avg(rh_col_annual_nep(0:19))

  ; rh_drh = rh_diff/avg(rh_col_annual_nep(0:19))

  ;  print("percent3="+h_dt/dq_q)
  ;  print("percent4="+rh_drh/dq_q)
  a = 17.67
  b = 243.5

  t_clm_2d = dim_avg_n_Wrap(t_col_annual, 0)

  t_clm_all  = t_clm_2d

  h = a*b/(t_clm_all+b)^2



  ; dq_q = q_diff/dim_avg_n_Wrap(q_col_annual(0:19,:,:))
  ; ; dq_q = qm_diff/avg(qm_col_nep(0:19))
  t_diff = dim_avg_n_Wrap(t_col_annual(20:,:,:), 0) - dim_avg_n_Wrap(t_col_annual(0:19,:,:), 0)

  q_dt = q_diff
  q_dt = h*(t_diff)*dim_avg_n_Wrap(q_col_annual(0:19,:,:),0)

  
  rh_diff = dim_avg_n_Wrap(rh_col_annual(20:,:,:),0) - dim_avg_n_Wrap(rh_col_annual(0:19,:,:),0)
  
  q_drh = q_diff
  q_drh =  dim_avg_n_Wrap(q_col_annual(0:19,:,:), 0)*rh_diff/dim_avg_n_Wrap(rh_col_annual(0:19,:,:),0)
  
;;


;;;;trend 
    q_rc = regCoef_n(ispan(1,yearnum,1), q_col_annual, 0,0)
    q_rc = (/q_rc*yearnum/)

    copy_VarCoords(q_col_annual(1,:,:), q_rc)

    t_rc = regCoef_n(ispan(1,yearnum,1), t_col_annual, 0,0)
    t_rc = (/t_rc*yearnum/)

    copy_VarCoords(t_col_annual(1,:,:), t_rc)

    rh_rc = regCoef_n(ispan(1,yearnum,1), rh_col_annual, 0,0)
    rh_rc = (/rh_rc*yearnum/)
    
    copy_VarCoords(rh_col_annual(1,:,:), rh_rc)

    q_dt = q_rc
    q_dt = h*(t_rc)*dim_avg_n_Wrap(q_col_annual,0)
    
    q_drh = q_rc
    q_drh =  dim_avg_n_Wrap(q_col_annual, 0)*rh_rc/dim_avg_n_Wrap(rh_col_annual,0)

    

;;;;

  
;;;;regional mean 
  ; q_test_mean = dim_avg_n_Wrap(q_test(:,{850},{0:20},{250:280}), (/1,2/))

  ; q_mean = dim_avg_n_Wrap(q(:,{850},{0:20},{250:280}), (/1,2/))

  ; q_test_clm = avg(q_test_mean)
  ; q_clm  = avg(q_mean)
;;;;

q_diff = q_diff*10^3
qm_diff = qm_diff*10^3

qm_fixt_diff = qm_fixt_diff*10^3
qm_fixrh_diff = qm_fixrh_diff*10^3

q_drh = q_drh*10^3 
q_dt = q_dt*10^3 

; percent1 = avg(q_dt({0:20},{250:280}))/avg(q_diff({0:20},{250:280}))
; percent2 = avg(q_drh({0:20},{250:280}))/avg(q_diff({0:20},{250:280}))

percent1 = avg(q_dt({0:20},{250:280}))/avg(10^3*q_rc({0:20},{250:280}))
percent2 = avg(q_drh({0:20},{250:280}))/avg(10^3*q_rc({0:20},{250:280}))

print("percent1="+percent1)
print("percent2="+percent2)


;;;

;==================================
; wks   = gsn_open_wks("pdf", "q_change pattern-wenxiu method")
wks   = gsn_open_wks("pdf", "q_change pattern-shuheng method")
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn=False

  res@lbLabelStride=2

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
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.1

  ; res@cnLevelSelectionMode  = "ExplicitLevels"
  ; res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
  res@lbLabelStride = 1

  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "BlueWhiteOrangeRed"
  res@gsnLeftString="(a) q change (NCEP_NCAR q data)"    
  res@gsnLeftStringFontHeightF = 0.015
  ; res@gsnLeftStringOrthogonalPosF  = -0.13
 
  plot = new(4,graphic)

  plot(0) = gsn_csm_contour_map(wks, q_diff, res)
  res@gsnLeftString="(b) q change (based on C-C relation)"    
  plot(1) = gsn_csm_contour_map(wks, qm_diff, res)


  ; res@gsnLeftString="(c) q change (fix relative humidity)"    
  ; plot(2) = gsn_csm_contour_map(wks, qm_fixrh_diff, res)
  ; res@gsnLeftString="(d) q change (fix temperature)"    
  ; plot(3) = gsn_csm_contour_map(wks, qm_fixt_diff, res)


  res@gsnLeftString="(c) q change (induced by t change)"    
  plot(2) = gsn_csm_contour_map(wks, q_dt, res)
  res@gsnLeftString="(d) q change (induced by Rh change)"    
  plot(3) = gsn_csm_contour_map(wks, q_drh, res)


  str1 = unique_string("dum1")  ; "unique_string" will return a unique    



  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=2.0
  pgres@gsLineDashPattern=0


  ;;;画box
  lat1 = 0
  lat2 = 20
  lon1 = 250
  lon2 = 280

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;

  str2 = ispan(21, 24, 1) +""
  do i =0,3 
   plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 

  ;  ;;;打点
  ;   res_pro=True
  ;   res_pro@gsnAddCyclic = False
  ;   res_pro@gsnDraw=False
  ;   res_pro@gsnFrame=False
  ;   res_pro@cnFillOn = True

  ;   res_pro@cnLevelSelectionMode = "ExplicitLevels"
  ;   res_pro@cnLevels = (/0,0.01/)
  ;   res_pro@lbLabelBarOn = False
  ;   res_pro@cnLinesOn=False
  ;   ; res_pro@cnMonoLineColor =False
  ;   ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

  ;   res_pro@cnMonoFillPattern = True
  ;   res_pro@cnFillPattern = 17
  ;   res_pro@cnFillDotSizeF = 0.002
  ;   res_pro@cnMonoFillColor = False
  ;   res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
  ;   ; res_pro@cnLevels = ispan(-6,6,1)
  ;   res_pro@cnLineLabelsOn=False
  ;   res_pro@tmYRMode = "Automatic"
  ;   res_pro@cnInfoLabelOn=False
  ;   ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ;   ; res_pro@cnLineThicknessF = 3.
  ;   ; res_pro@cnFillOn        = True
  ;   iplot = new(4,graphic)
  ;    i=0
  ;   iplot(i) = gsn_csm_contour(wks,prob_mam,res_pro)
  ;    overlay(plot(i), iplot(i))
  ; ; end do
pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01
gsn_panel(wks, plot, (/2,2/), pres)

end

















