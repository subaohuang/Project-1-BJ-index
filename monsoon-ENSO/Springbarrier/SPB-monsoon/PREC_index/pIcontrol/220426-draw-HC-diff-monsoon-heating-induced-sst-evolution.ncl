;;*************************************************
;; draw circulation pattern during monsoon year 
;;*************************************************************
function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)
 return confi_pv_east
 end

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum/),float)
    do i=0,11   
    data_12(i,:) = data(i::12)
    end do 
return data_12
end 
function composite_mean(ssta_month_12,ind1)
begin 
   lat1 = -5
   lat2 = 5
   lon1 = 190
   lon2 = 240
   ssta_month_str = dim_avg_n_Wrap(ssta_month_12(:,ind1,{lat1:lat2},{lon1:lon2}),(/1,2,3/))
    ; u_month_str = dim_avg_n_Wrap(u_850_month_12(:,ind2,{lat1:lat2},{lon1:lon2}),(/1,2,3/))
    ; v_month_str = dim_avg_n_Wrap(v_850_month_12(:,ind3,{lat1:lat2},{lon1:lon2}),(/1,2,3/))

return [/ssta_month_str/]
end 


;;*************************************************
;; main program 
;;*************************************************************
begin
   
  stringlesft = (/"(a) El Nino","(b) La Nina","(c) Normal"/)
  str1= (/"El Nino","La Nina","Normal"/)

  plot = new(2,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "1std-DIFF-HC-ssta-evolution")
  
   
  stss = ispan(1,2,1)+""

  data_total = new((/2,6,18/),float)


  ;;case1 = [/wet_str_warm_year,wet_weak_warm_year,dry_str_warm_year,dry_weak_warm_year,nor_str_warm_year,nor_weak_warm_year/]

  ;;case2 = [/wet_str_cold_year, wet_weak_cold_year, dry_str_cold_year, dry_weak_cold_year,nor_str_cold_year,nor_weak_cold_year/]



do k=0,1
  print("***********start*********************"+k)

  opt1="1std-"+str1(k)+""+"-diff-intensity-same-monsoon-SSTevolution.txt"
  data = asciiread(opt1, -1, "float")

  data_new = reshape(data, (/6,18/))
  data_total(k,:,:) = data_new
end do 


;;;; wet - dry 

 highhc = -1*(/(data_total(0,0,:) - data_total(0,2,:))/2.0 +  (data_total(1,0,:) - data_total(1,2,:))/2.0/)/2.0

 lowhc  = -1*(/(data_total(0,1,:) - data_total(0,3,:))/2.0 +  (data_total(1,1,:) - data_total(1,3,:))/2.0/)/2.0
;;;;


;;;;
  std_nino = 0.846
    
k = 0

;;**************************************************************
;画图
;;***************************************************************
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  res@tmXTOn         = False
  res@tmYROn = False


  res@trXMinF        = 0
  res@trXMaxF        = 17

  res@gsnYRefLine = 0.0

  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 17, 1)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M(+1)","A(+1)","M(+1)","J(+1)","J(+1)","A(+1)"/)
  res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
   res@tmYLMode = "Explicit"

   res@trYMinF        = -1.0
   res@trYMaxF        = 1.0

   res@tmYLValues := fspan(-1.0, 1.0, 11)
   res@tmYLLabels := fspan(-1.0, 1.0, 11)


  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.02

  res@tiXAxisString = "Month"
  res@tiXAxisFontHeightF = 0.02

  res@xyMonoMarkLineMode = True
  res@xyMarkLineMode = (/"MarkLines"/)
  res@xyMonoMarker = True
  res@xyMarker = (/1/)
  res@xyMonoMarkerSize = True
  res@xyMarkerSizes = (/0.03/)
  

  res@xyMarkerColors= (/"red","blue"/)

  res@xyMonoMarkerThickness = True
  res@xyMarkerThicknesses = (/4.0/)

  res@xyMonoLineThickness = True
  res@xyLineThicknesses  = (/2.0/)

  res@xyLineColors       = res@xyMarkerColors
  res@xyMonoDashPattern = True
  res@xyDashPatterns     = (/1/)

  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = stringlesft(k)
  plot(k)  = gsn_csm_xy (wks, ispan(0,17,1), (/highhc,lowhc/), res)


  respl = True
  respl@gsLineColor = "gray"
  respl@gsLineDashPattern = 1.0

  x =  1.0*res@tmXBValues
  y = x
 
  y = 0.43*std_nino
  plot@$stss(k)$ =  gsn_add_polyline(wks, plot(k),x, y, respl)


  ;;添加legend 
    gres = True
    if k.eq.0 then 
    gres@YPosPercent = 25.    ; expressed as %, 0->100, sets position of top border of legend 
    else
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    end if
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True

    textres@lgLabels = (/"High HC ","Low HC "/)

    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 
    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   

pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/2,1/), pres)
end 