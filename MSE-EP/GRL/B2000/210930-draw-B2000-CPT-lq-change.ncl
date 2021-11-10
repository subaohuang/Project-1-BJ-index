function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      ; simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
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

begin

    yearstart = 25601
    yearend = 28012

    lat1 = -50
    lat2 = 50
    lon1 = 0
    lon2 = 360

;;;;read contrl data
    modelname = "B2000_alt_north_year_CTRL"
    diri = "/home/yangsong3/data-observation/linshh/CESM-data/B/"+modelname+"/"
    filename= modelname+".cam.h1.0251-0280_chazhi.nc"

    f1 = addfile(diri+filename, "r")
   
    time = f1->time 
    date = cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)


    plev = f1->lev({1000:300})


    t_ct = f1->T(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q_ct = f1->Q(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    ; q = q/1000.0  ;;;convert to kg/kg
    hgt_ct = f1->Z3(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

;;;


;;;;read data
    ; diri = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year/"
    ; diri = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_db_year_80_20/"
    diri = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_db_year_80_20/"
    filename="/lsh_B2000_alt_north_db_year_80_20.cam.h1.0251-0280_chazhi.nc"
    modelname = "lsh_B2000_alt_north_db_year_80_20"

     wks = gsn_open_wks("pdf","JFMAM-CpT和Lq的变化-CT-"+modelname+"-1980-2020")          ; send graphics to PNG file

    f1 = addfile(diri+filename, "r")
   
    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)

    t_sen = f1->T(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    q_sen = f1->Q(timeind,{1000:300},{lat1:lat2},{lon1:lon2})

    ; q = q/1000.0  ;;;convert to kg/kg
    hgt_sen = f1->Z3(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
;;;


;;;;
    Cp = 1004
    Lv = 2.5*10^6
    g = 9.8
    
    hgt_ct = (/hgt_ct*9.8/)
    t_ct  = Cp*t_ct
    q_ct  = Lv*q_ct
    
    mse_ct = hgt_ct
    mse_ct = hgt_ct + t_ct + q_ct

    hgt_sen = (/hgt_sen*9.8/)
    t_sen  = Cp*t_sen
    q_sen  = Lv*q_sen
    season  = (/"MAM"/)
    
    mse_sen = hgt_ct
    mse_sen = hgt_sen + t_sen + q_sen


;;;;


;;;MSE 
    ptop = min(plev)
    psfc = 1013.0
    dz = plev
    dz = dz
    dz = dpres_plevel(plev, psfc, ptop, 0)
    dz&lev = plev
    printVarSummary(dz)
    
;;;;
 

;;;;;;;;;;;;;;;;;;;
 t_ct_col =  dim_sum_wgt_n_Wrap(t_ct(:,{1000:500},:,:), dz({1000:500}), 1, 1)
 

 t_ct_col_12 = split_12(t_ct_col)


 q_ct_col =  dim_sum_wgt_n_Wrap(q_ct(:,{1000:500},:,:), dz({1000:500}), 1, 1)
 

 q_ct_col_12 = split_12(q_ct_col)


 mse_ct_col =  dim_sum_wgt_n_Wrap(mse_ct(:,{1000:300},:,:), dz({1000:300}), 1, 1)
 

 mse_ct_col_12 = split_12(mse_ct_col)

;;;;;


;;;
 t_sen_col =  dim_sum_wgt_n_Wrap(t_sen(:,{1000:500},:,:), dz({1000:500}), 1, 1)
 

 t_sen_col_12 = split_12(t_sen_col)


 q_sen_col =  dim_sum_wgt_n_Wrap(q_sen(:,{1000:500},:,:), dz({1000:500}), 1, 1)
 

 q_sen_col_12 = split_12(q_sen_col)


 mse_sen_col =  dim_sum_wgt_n_Wrap(mse_sen(:,{1000:300},:,:), dz({1000:300}), 1, 1)
 
 mse_sen_col_12 = split_12(mse_sen_col)

;;;;


 t_ct_col_annual = month_to_season(t_ct_col, "MAM")
 q_ct_col_annual = month_to_season(q_ct_col, "MAM")
 mse_ct_col_annual = month_to_season(mse_ct_col, "MAM")



 
 t_ct_col_annual = (/dim_avg_n_Wrap(t_ct_col_12(0:4,:,:,:), 0)/)
 q_ct_col_annual = (/dim_avg_n_Wrap(q_ct_col_12(0:4,:,:,:), 0)/)
 mse_ct_col_annual = (/dim_avg_n_Wrap(mse_ct_col_12(0:4,:,:,:), 0)/)


 t_sen_col_annual = month_to_season(t_sen_col, "MAM")
 q_sen_col_annual = month_to_season(q_sen_col, "MAM")
 mse_sen_col_annual = month_to_season(mse_sen_col, "MAM")
 
 t_sen_col_annual = (/dim_avg_n_Wrap(t_sen_col_12(0:4,:,:,:), 0)/)
 q_sen_col_annual = (/dim_avg_n_Wrap(q_sen_col_12(0:4,:,:,:), 0)/)
 mse_sen_col_annual = (/dim_avg_n_Wrap(mse_sen_col_12(0:4,:,:,:), 0)/)
 
 

 rc_A_mam  = t_ct_col_annual(1,:,:)
 rc_A_mam  = (/dim_avg_n_Wrap(t_sen_col_annual,0) - dim_avg_n_Wrap(t_ct_col_annual,0)/)


 rc_B_mam  = q_ct_col_annual(1,:,:)
 rc_B_mam  = (/dim_avg_n_Wrap(q_sen_col_annual,0) - dim_avg_n_Wrap(q_ct_col_annual,0)/)

 rc_A_mam = (/rc_A_mam/10^6/)
 rc_B_mam = (/rc_B_mam/10^6/)


 rc_C_mam  = mse_ct_col_annual(1,:,:)
 rc_C_mam  = (/dim_avg_n_Wrap(mse_sen_col_annual,0) - dim_avg_n_Wrap(mse_ct_col_annual,0)/)
 
 rc_C_mam = (/rc_C_mam/10^6/)

 prob_A = t_test(t_ct_col_annual,t_sen_col_annual)
 prob_B = t_test(q_ct_col_annual,q_sen_col_annual)
 prob_C = t_test(mse_ct_col_annual,mse_sen_col_annual)
  

 printVarSummary(rc_A_mam)
 printVarSummary(rc_B_mam)
  printVarSummary(rc_C_mam)
 


 

;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=0.01
  res@cnFillOn=True
  res@cnLinesOn=False
  res@pmLabelBarOrthogonalPosF = 0.3

  res@cnLineLabelsOn=False
  ; res@lbLabelBarOn=True
  res@lbLabelBarOn= False
  ; res@lbLabelBarOn=True

  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = -20
  res@mpMaxLatF             = 20
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  

 
  res@tmYLValues = (/-10,0,10,20/) 
  res@tmYLLabels = (/"10S","0","10N","20N"/) 

  res@tmYLMinorValues = ispan(-15,20,5)


  res@mpFillOn = True
  ; res@mpLandFillColor = "gray"
  ; res@mpOceanFillColor = "Transparent"
  ; res@mpMonoFillColor = True
  res@mpFillDrawOrder = "PostDraw"


 
  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  res@lbBoxLinesOn              = True
  res@lbBoxLineThicknessF       = 0.


  res@tmXBLabelFontHeightF=0.01
  res@tmYLLabelFontHeightF=0.01
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0


  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.8*2
  res@cnMaxLevelValF       =   0.8*2
  res@cnLevelSpacingF      =  0.2*2

  ; res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour level
  ; res@cnLevels = (/-2.0,-1.5,-1.0,1.0,1.5,2.0/)

  res@lbLabelFontHeightF = 0.01
  res@lbLabelStride = 1
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueDarkRed18"
  res@cnFillPalette = "BlueWhiteOrangeRed"

;   res@mpShapeMode = "FreeAspect"
;   res@vpHeightF = 0.2
;   res@vpWidthF = 0.7
; ; 
  ; res@gsnLeftString="Annual mean MSE trend " 
  plot = new(3,graphic)
  res@gsnRightString = ""

  res@gsnLeftString="(a) JFMAM 1000-300hPa integrated MSE (ATL warming - CTRL) " 
  plot(0) = gsn_csm_contour_map(wks, rc_C_mam, res)
  res@gsnLeftString="(b) JFMAM 1000-300hPa integrated CpT (ATL warming - CTRL) " 
  plot(1) = gsn_csm_contour_map(wks, rc_A_mam, res)
  res@gsnLeftString="(c) JFMAM 1000-500hPa integrated Lq (ATL warming - CTRL) " 
  plot(2) = gsn_csm_contour_map(wks, rc_B_mam, res)

  str1 = unique_string("dum1")  ; "unique_string" will return a unique    

  pgres=True
  pgres@gsLineColor="black"

  pgres@gsLineThicknessF=1.0
  pgres@gsLineDashPattern=0

  ;;;画box
  lat1 = 0
  lat2 = 20
  lon1 = 250
  lon2 = 280

  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
  ;;;
  add  = ispan(1,3,1)+""
  do i=0,2
  plot@$add(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do
  ;;;打点

  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnFillOn = True

  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLevels = (/0,0.1/)
  res_pro@lbLabelBarOn = False
  res_pro@cnMonoFillColor = False
  res_pro@cnFillColors = (/"Transparent","White","Transparent"/)


  res_pro@cnLinesOn=False
  ; res_pro@cnMonoLineColor =False
  ; res_pro@cnLineColors = ((/"Transparent","black","Transparent"/))

  res_pro@cnMonoFillPattern = True
  res_pro@cnFillPattern = 17
  res_pro@cnFillDotSizeF = 0.002
  
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True
   
  iplot = plot 
  iplot(0) = gsn_csm_contour(wks,prob_C,res_pro)
  iplot(1) = gsn_csm_contour(wks,prob_A,res_pro)
  iplot(2) = gsn_csm_contour(wks,prob_B,res_pro)
  
  do i=0,2
  overlay(plot(i), iplot(i))
  end do

    ;;  

  pres = True
  pres@gsnPanelLabelBar  = True
  pres@pmLabelBarOrthogonalPosF = -0.02
  ; pres@pmLabelBarParallelPosF = 0.02
  pres@pmLabelBarHeightF = 0.05
  pres@pmLabelBarWidthF = 0.6

  gsn_panel(wks, plot, (/3,1/), pres)

end



