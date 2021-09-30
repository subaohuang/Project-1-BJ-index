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

    u_ct = f1->U(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    v_ct = f1->V(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    
    div_ct = u_ct
    div_ct = uv2dvF(u_ct,v_ct)

    delete(u_ct)
    delete(v_ct)

    
;;;


;;;;read data
    diri = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year/"
    filename="B2000_alt_north_year.cam.h1.0251-0280_chazhi.nc"

    f1 = addfile(diri+filename, "r")
   
    time := f1->time 
    date := cd_calendar(time, 1)
    delete(time)

    timeind :=ind(date.le.yearend.and.date.ge.yearstart)

    u_sen = f1->U(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    v_sen = f1->V(timeind,{1000:300},{lat1:lat2},{lon1:lon2})
    
    div_sen = u_sen
    div_sen = uv2dvF(u_sen,v_sen)
    delete(u_sen)
    delete(v_sen)
 
;;;


;;;;;;;;;;;;;;;;;;;
  div_ct_col =  dim_avg_n_Wrap(div_ct(:,{1000:900},:,:),1)
  ; div_ct_col =  div_ct(:,{1000},:,:)

  div_ct_col_12 = split_12(div_ct_col)

  div_sen_col =  dim_avg_n_Wrap(div_sen(:,{1000:900},:,:),1)
  ; div_sen_col = div_sen(:,{1000},:,:)
  div_sen_col_12 = split_12(div_sen_col)
;;;;;


 div_ct_col_annual = month_to_season(div_ct_col, "MAM")
 
 
 div_ct_col_annual = (/dim_avg_n_Wrap(div_ct_col_12(0:4,:,:,:), 0)/)



 div_sen_col_annual = month_to_season(div_sen_col, "MAM")
 
 
 div_sen_col_annual = (/dim_avg_n_Wrap(div_sen_col_12(0:4,:,:,:), 0)/)

 

 rc_A_mam  = div_ct_col_annual(1,:,:)
 rc_A_mam  = (/dim_avg_n_Wrap(div_sen_col_annual,0) - dim_avg_n_Wrap(div_ct_col_annual,0)/)

 rc_A_mam = (/rc_A_mam*-1*10^6/)

 prob_A = t_test(div_ct_col_annual,div_sen_col_annual)
  



 printVarSummary(rc_A_mam)

 print("rc_A_mam="+avg(rc_A_mam({0:20},{250:280})))

 
 



  wks = gsn_open_wks("pdf","div的变化-CT-SEN-1980-2020")          ; send graphics to PNG file

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
  
  res@mpMinLatF             = -40
  res@mpMaxLatF             = 40
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0
  

 
  ; res@tmYLValues = (/-10,0,10,20/) 
  ; res@tmYLLabels = (/"10S","0","10N","20N"/) 

  ; res@tmYLMinorValues = ispan(-15,20,5)


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
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.1

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
  plot = new(1,graphic)

  res@gsnLeftString="(a) JFMAM 1000-300hPa integrated CpT (ATO - CTRL) " 
  plot = gsn_csm_contour_map(wks, rc_A_mam, res)
  
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
  add  = ispan(1,2,1)+""
  i=0
  plot@$add(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  
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
  iplot(0) = gsn_csm_contour(wks,prob_A,res_pro)
  
  
  overlay(plot(0), iplot(0))
  

    ;;  

  pres = True
  pres@gsnPanelLabelBar  = True
  pres@pmLabelBarOrthogonalPosF = -0.02
  ; pres@pmLabelBarParallelPosF = 0.02
  pres@pmLabelBarHeightF = 0.05
  pres@pmLabelBarWidthF = 0.6

  gsn_panel(wks, plot, (/2,1/), pres)

end



