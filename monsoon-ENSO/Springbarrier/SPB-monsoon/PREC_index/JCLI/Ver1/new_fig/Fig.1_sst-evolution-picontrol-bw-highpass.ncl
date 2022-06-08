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
   

  ntilde  = "n~H-13V2F35~D~FV-2H3~"

  stringlesft = (/"(d) El Ni"+ntilde+"o","(e) La Ni"+ntilde+"o","(f) Normal"/)
  
  str1 = (/"El Nino","La Nina","Normal"/)

  plot = new(3,graphic)
  midplot = plot

  wks = gsn_open_wks("pdf", "BW-Highpass-CESM2-picontrol-ssta-evolution-during-diffmonsoon")
  
   
  stss = ispan(1,3,1)+""

do k=0,2 
  print("***********start*********************"+k)


  ssta_str_evo = new((/18/),float)
  ssta_weak_evo = new((/18/),float)
  ssta_nor_evo = new((/18/),float)

  ; diri = "/home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/pIcontrol/"
  ; opt1 =diri + str1(k)+""+"diffmonsoon-SSTevolution-highpass.txt"

  opt1 ="highpass_bw_" + str1(k)+""+"diffmonsoon-SSTevolution.txt"

  data = asciiread(opt1, -1, "float")
  ssta_str_evo  = data(0:17)
  ssta_weak_evo  = data(18:35)
  ssta_nor_evo  = data(36:)

  std_nino = 0.567
    
;;**************************************************************
;画图
;;***************************************************************
  num_x = 16

  res = True
  
  res@gsnDraw        = False
  res@gsnFrame       = False

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  ; res@tmXTOn         = False
  ; res@tmYROn = False


  res@trXMinF        = 0
  res@trXMaxF        = num_x-1



  res@gsnYRefLine = 0.0


  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025

  res@gsnLeftStringFontHeightF   = 0.025
  res@gsnRightStringFontHeightF  = 0.025

  res@tmBorderThicknessF = 2.0

  res@tmXBMode = "Explicit"
  ; res@tmXBValues = ispan(0, num_x-1, 1)
  res@tmXBValues = ispan(0, num_x-1, 3)
  res@tmXBMinorValues = ispan(0, num_x-1, 1)

  

  res@tmXBLabels = (/"Mar","Jun","Sep","Dec","Mar(1)","Jun(1)"/)
  ; res@tmXBLabels = (/"M","M","J","A","S","N","D","F(+1)","M(+1)","M(+1)",/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J(+1)","F(+1)","M(+1)","A(+1)","M(+1)"/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J1","F1","M1","A1","M1","J1","J1","A1"/)
  ; res@tmXBLabels = (/"M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A"/)
  res@tmYLMode = "Explicit"



  if k.eq.0 then 
   ; res@trYMinF        = -0.6
   ; res@trYMaxF        = 1.6

   res@trYMinF        = -1.0
   res@trYMaxF        = 2.0

   ; res@tmYLValues := ispan(-6, 16, 2)*0.1
   ; res@tmYLLabels := ispan(-6, 16, 2)*0.1
  end if 

  if k.eq.1 then 
   ; res@trYMinF        = -1.2
   ; res@trYMaxF        = 0.4
   res@trYMinF        = -2.0
   res@trYMaxF        = 1.0

   ; res@tmYLValues := ispan(-12, 4, 2)*0.1
   ; res@tmYLLabels := ispan(-12, 4, 2)*0.1
  end if 

  if k.eq.2 then 

   res@trYMinF        = -0.6
   res@trYMaxF        = 0.6

   ; res@tmYLValues := ispan(-6, 6, 2)*0.1
   ; res@tmYLLabels := ispan(-6, 6, 2)*0.1
  end if 

  if k.eq.0 then 
  res@tiYAxisString = "SSTA [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.025
  else
   res@tiYAxisString = "" 
  end if
  ; res@tiXAxisString = "Month"
  ; res@tiXAxisFontHeightF = 0.025


  res@xyMarkLineModes = (/"MarkLines","MarkLines","MarkLines"/)
  res@xyMarkers = (/1,1,1/)
  res@xyMarkerSizes = (/0.03,0.03,0.03/)
  res@xyMarkerColors= (/"red","blue","black"/)


  res@xyMarkerThicknesses = (/4.0,4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0,2.0/)
  res@xyLineColors       = (/"red","blue","black"/)
  res@xyDashPatterns     = (/1,1,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = stringlesft(k)
  res@gsnRightString = "PiControl"

  if k.ne.2 then 
   plot(k)  = gsn_csm_xy (wks, ispan(0,15,1), (/ssta_str_evo(0:15),ssta_weak_evo(0:15),ssta_nor_evo(0:15)/), res)
  else
   plot(k)  = gsn_csm_xy (wks, ispan(0,15,1), (/ssta_str_evo(0:15),ssta_weak_evo(0:15)/), res)
  end if 


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
    ; if k.eq.0 then 
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
    ; else
    ; gres@YPosPercent = 25.    ; expressed as %, 0->100, sets position of top border of legend 
    ; end if
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors

    if k.ne.2 then 
    lineres@lgLineColors       = res@xyLineColors   ; colors chosen
    else
    lineres@lgLineColors := (/"red","blue"/)
    end if
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line
    lineres@lgDashIndexes = res@xyDashPatterns
    textres = True
    
    if k.ne.2 then 
    textres@lgLabels := (/"Strong monsoon","Weak monsoon", "Normal monsoon"/)
    else
    textres@lgLabels := (/"Strong monsoon","Weak monsoon"/)
    end if
    textres@lgLabelFontHeightF = 0.015
    ; textres@lgLabelOpacity = (/0/) 
    midplot(k) = simple_legend(wks,plot(k),gres,lineres,textres)   
   

end do
pres = True
; pres@gsnPanelLabelBar  = True
; pres@lbLabelStride=1
gsn_panel(wks, plot, (/1,3/), pres)
end 