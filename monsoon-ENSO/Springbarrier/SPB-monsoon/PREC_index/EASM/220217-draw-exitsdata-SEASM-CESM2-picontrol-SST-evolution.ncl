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

  plot = new(3,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "CESM2-picontrol-ssta-evolution-during-SEASM-diffmonsoon")
  
   
  stss = ispan(1,3,1)+""

do k=0,2 
  print("***********start*********************"+k)


  ssta_str_evo = new((/18/),float)
  ssta_weak_evo = new((/18/),float)
  ssta_nor_evo = new((/18/),float)

  
    
  opt1=str1(k)+""+"-EASM-diffmonsoon-SSTevolution.txt"
  data = asciiread(opt1, -1, "float")
  ssta_str_evo  = data(0:17)
  ssta_weak_evo  = data(18:35)
  ssta_nor_evo  = data(36:)

  std_nino = 0.846
    

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

  if k.eq.2 then 
   res@trYMinF        = -1.0
   res@trYMaxF        = 1.0

   res@tmYLValues := fspan(-1.0, 1.0, 11)
   res@tmYLLabels := fspan(-1.0, 1.0, 11)
  else
   res@trYMinF        = -2.0
   res@trYMaxF        = 2.0

   res@tmYLValues := fspan(-2.0, 2.0, 11)
   res@tmYLLabels := fspan(-2.0, 2.0, 11)
  end if 

  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.02

  res@tiXAxisString = "Month"
  res@tiXAxisFontHeightF = 0.02


  res@xyMarkLineModes = (/"MarkLines","MarkLines","MarkLines"/)
  res@xyMarkers = (/1,1,1/)
  res@xyMarkerSizes = (/0.03,0.03,0.03/)
  res@xyMarkerColors= (/"red","blue","black"/)



  res@xyMarkerThicknesses = (/4.0,4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0,2.0/)
  res@xyLineColors       = res@xyMarkerColors
  res@xyDashPatterns     = (/1,1,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = stringlesft(k)
  ; plot(k)  = gsn_csm_xy (wks, ispan(0,17,1), (/ssta_str_evo,ssta_weak_evo,ssta_nor_evo/), res)
  plot(k)  = gsn_csm_xy (wks, ispan(0,17,1), (/ssta_weak_evo,ssta_str_evo,ssta_nor_evo/), res)


  respl = True
  respl@gsLineColor = "gray"
  respl@gsLineDashPattern = 1.0

  x =  1.0*res@tmXBValues
  y = x

  if k.eq.0 then 
  y = 0.43*std_nino
  end if 

  if k.eq.1 then 
  y = -0.43*std_nino
  end if 
  
  if k.ne.2 then 
  plot@$stss(k)$ =  gsn_add_polyline(wks, plot(k),x, y, respl)
  end if


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

    textres@lgLabels = (/"Strong SEASM","Weak SEASM", "Normal SEASM"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 

    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   
   
end do
pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/2,2/), pres)
end 