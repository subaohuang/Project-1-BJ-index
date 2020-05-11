;;********************************************************
;;画GPCP的时间序列图

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
;;;region
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155
 
 
;;;read Q1 monthly data from linshuheng
  fils = systemfunc("ls /home/ys17-19/lsh/Project/SCS-rain/Q_time_era.nc")
  f1 = addfile(fils, "r")
  date = cd_calendar(f1->time, 1)

  timeind = ind(date.le.201612.and.date.ge.197901)

  q1 = f1->Q_time(timeind,{1000:200},{lat1:lat2},{lon1:lon2})
  level = f1->level({1000:200})
  
  printVarSummary(q1)

  q1_annual = month_to_annual(q1, 1)
  ; q1_annual = month_to_season(q1, "MAM")

  dim = dimsizes(q1_annual)

  q1_annual_rc = regCoef_n(ispan(1,dim(0),1), q1_annual,0,0)
  copy_VarCoords(q1_annual(1,:,:,:), q1_annual_rc)  
  q1_annual_rc = q1_annual_rc*dim(0)
;;;;

;;;;read GPCP data
  f2  = addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
  date := cd_calendar(f2->time, 1)
  timeind := ind(date.ge.197901.and.date.le.201612)
  p = f2->precip(timeind,{lat1:lat2},{lon1:lon2})
  p_annual = month_to_annual(p, 1)
;;;


;;;;read CFSR diabatic and deep convection data 

    f2=addfile("/home/ys17-19/lsh/Project/SCS-rain/diabatic_and_deepcon_heating_CFSR_ANN_SCS_timeseries.nc", "r")
  
    deepcon_q_time = f2->deepcon_q_time(:,{100000:20000})
    diab_q_time = f2->diab_q_time(:,{100000:20000})
    lev = f2->lev({100000:20000})
    printVarSummary(deepcon_q_time)
  
    ; printMinMax(diab_q_time, 1)
    deepcon_q_time = deepcon_q_time*86400.0
    diab_q_time = diab_q_time*86400.0
;;;;

;;get region mean as index     
   q1_scs = dim_avg_n_Wrap(q1_annual, (/2,3/))
   p_scs = dim_avg_n_Wrap(p_annual, (/1,2/))
;;;

;;;垂直平均
  
  dp=dpres_plevel(lev,100000,min(lev),0)              ; caculate dp
  dp_cfm=conform_dims(dimsizes(deepcon_q_time),dp,1)          ; confirm dp to caculate sum(var*dp)

  deepcon_q_wgt = (/dim_sum_n(dp_cfm*deepcon_q_time,1)/dim_sum(dp)/)               ; sum(var*dp)
  diab_q_wgt = (/dim_sum_n(dp_cfm*diab_q_time,1)/dim_sum(dp)/)               ; sum(var*dp)
;;;

;;;垂直平均—q1
    dp:=dpres_plevel(level,1000,min(level),0)              ; caculate dp
    dp_cfm:=conform_dims(dimsizes(q1_scs),dp,1)          ; confirm dp to caculate sum(var*dp)
    q1_wgt = (/dim_sum_n(dp_cfm*q1_scs,1)/dim_sum(dp)/)               ; sum(var*dp)
;;;

printVarSummary(diab_q_wgt)

   ; ssta_index_rc1 = regCoef_n(ispan(0,18,1), ssta_index(0:18), 0, 0)
   ; ssta_index_rc2 = regCoef_n(ispan(0,yearnum-18-1,1), ssta_index(18:), 0, 0)
   
   ; line1 = ssta_index_rc1*ispan(0,18,1) + ssta_index_rc1@yintercept
   ; line2 = ssta_index_rc2*ispan(0,yearnum-17,1) + ssta_index_rc2@yintercept


   q1_wgt_std = dim_standardize(q1_wgt, 1)
   p_scs_std = dim_standardize(p_scs, 1)
   diab_q_wgt_std = dim_standardize(diab_q_wgt, 1)
   deepcon_q_wgt_std = dim_standardize(deepcon_q_wgt, 1)

  data_all = new((/4,dimsizes(timeind)/12/),"float")
  data_all(0,:)=q1_wgt_std
  data_all(1,:)=p_scs_std
  data_all(2,0:31)=diab_q_wgt_std
  data_all(3,0:31)=deepcon_q_wgt_std

  




print("***********画时间序列***************")

;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "q1-CFSR_diabatic-deepconvection-GPCP降水的时间序列-2016")
  plot=new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF =  endyear
  resxy@trXMinF = startyear
  resxy@tmXBValues = ispan(startyear,endyear,5)
  resxy@tmXBLabels = ispan(startyear,endyear,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(startyear,endyear,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknesses =  (/ 4,4,4,4/)          ; make second line thicker
  
  resxy@xyLineColors       =  (/"black","gray","red","blue"/)          ; change line color
  resxy@xyMonoDashPattern = False
  resxy@xyDashPatterns =(/0,0,0,0/)
  resxy@xyMarkLineModes    =  (/"lines","lines","lines","lines"/)
  
  time = ispan(startyear,endyear,1)
  plot = gsn_csm_xy(wks,time,data_all, resxy)



  ; ;;;添加两条趋势线
  ; res = True
  ; res@gsLineColor = "red"
  
  ;  dum1 = gsn_add_polyline(wks, plot, ispan(1979,1997,1), line1, res)
  ;  dum2 = gsn_add_polyline(wks, plot, ispan(1998,2013,1), line2, res)


  ;;添加legend 
  gres = True
  gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                    ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 5    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

  lineres = True
  ; lineres@lgLineColors = (/"black"，/) ; line colors
  lineres@lgLineColors       =   resxy@xyLineColors 
  lineres@lgLineThicknesses = resxy@xyLineThicknesses                 ; line thicknesses
  lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
  lineres@lgDashIndexes = resxy@xyDashPatterns

  textres = True
  textres@lgLabels = (/"Q1","Precipitation","Diabatic heating","Deep convection latent heating"/)     ;加线的标注

  plot = simple_legend(wks,plot,gres,lineres,textres)   
  draw(plot)
  frame(wks)

  end
