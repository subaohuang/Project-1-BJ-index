
begin 
 
    model_name = (/"B2000_f19g16_CP_CTRL","lsh_B2000_WNP_heating_IDLthree"/)
    dim_model = dimsizes(model_name)

    yearstart = 25601
    yearend = 28012
    lat1 = -30
    lat2 = 30
    lon1 = 120
    lon2 = 290
    box_divide = 190
;;;
do k = 0,dim_model-1
;;;;
   model = model_name(k)

   if model.eq."B2000_f19g16_CP_CTRL" then 
     diri ="/home/yangsong3/data-observation/linshh/CESM-data/B/"
   else 
     diri ="/home/yangsong3/CMIP6/linshh/CESM-data/B/"
   end if

    print("*******************start calculate"+" "+model+"**********************************************")
 
    ; fils1 :=  systemfunc("ls /home/yangsong3/data-observation/linshh/CESM-data/B/"+model+"/"+model+".TEMP.fvgrid.0251-0280.nc")
    ; fils1 := systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/"+model+"/"+model+".TEMP.fvgrid.0251-0280.nc")
   
    fils1 :=  systemfunc("ls "+diri+model+"/"+model+".TEMP.fvgrid.0251-0280.nc")


    f1 := addfile(fils1, "r")

    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    time = (/time-1/) 
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)   

    tos := f1->TEMP(timeind,{0:5000},{lat1:lat2},:)
    tos!1 = "lev"
    
    tos = (/dtrend_msg_n(ispan(1,dimsizes(tos&time),1),tos,False,True,0)/)
    tos_ar = rmMonAnnCycLLLT(tos(lev|:,lat|:,lon|:,time|:))
    tos_a = tos_ar(time|:,lev|:,lat|:,lon|:)
    tos_clm = dim_avg_n_Wrap(tos, 0)  ;;气候态的海温

;;;

;;parameter 
     pi = get_pi("float")
     dx := tos&lat
     dx = 6378388.*cos(0.0174533*tos&lat)
     dx!0 = "lat"
     copy_VarCoords(tos(1,1,:,1), dx)
     dy := tos&lon
     dy = (pi/180)*1*6378388.
     dy!0 = "lon"
     copy_VarCoords(tos(1,1,1,:), dy)

     
     depth = tofloat(tos&lev)/100.0
     dz = depth

     dz(1:)= depth(1:) - depth(0:dimsizes(depth)-2)
     dz(0)= depth(0) 
;;;

;;;;计算T_clm的纬向梯度
  
     dtdxdy = grad_latlon_cfd(tos_clm, tos_clm&lat, tos_clm&lon, False, False)
     dtdx = dtdxdy[1]
     copy_VarCoords(tos_clm, dtdx)
;;;;


;;;计算东box区域dtdx平均，混合层深度取50 ，185.2E - 80w -5-5N
     dtdx_vol = wgt_volave(dtdx(:,{-5:5},{box_divide:280}), dz, 1, 1, 0)
;;;
   
;;;计算东box区域异常海温平均，混合层深度取50 ，185.2E - 80w -5-5N
     Tbox_a = wgt_volave(tos_a(:,:,{-5:5},{box_divide:280}), dz, 1, 1, 0)
;;;


print("*********************Thermal Damping********************************")
  
;;;;读入Q数据 /home/yangsong3/CMIP6/linshh/CESM-data/B/
    ; fils5 :=  systemfunc("ls /home/yangsong3/data-observation/linshh/CESM-data/B/"+model+"/"+model+".SHF.fvgrid.0251-0280.nc")
    ; fils5 :=  systemfunc("ls /home/yangsong3/CMIP6/linshh/CESM-data/B/"+model+"/"+model+".SHF.fvgrid.0251-0280.nc")
    fils5 :=  systemfunc("ls "+diri+model+"/"+model+".SHF.fvgrid.0251-0280.nc")
    f5 := addfile(fils5, "r")
   
    time:=f5->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    hfds = f5->SHF(timeind,{lat1:lat2},:)


    Cp = 3940;1996
    rho0 = 1022.4;1026
    Hm=50.0

  ; density = 1025  ;;kg/m3
  ; cp  = 3850         ;;j/kg℃
  ; coef  = density*cp
    hfds = (/hfds/(Cp*rho0*Hm)/) ;; k/s

    hfds = (/dtrend_msg_n(ispan(1,dimsizes(hfds&time),1),hfds,False,True,0)/)
    hfds_a = rmMonAnnCycTLL(hfds)
    hfds_clm = dim_avg_n_Wrap(hfds, 0)


;;;;

;;;hfds_a_box 
    hfds_a_box = hfds_a(:,{-30:30},:)

;;;alfas  hfds_a_box = -alfas*Tbox_a
    alfas  =  regCoef_n(-1*Tbox_a, hfds_a_box, 0, 0) ;; 区域洋流对热带纬向风的响应
 
;;;

;;**TD = -alfas

   if model.eq."B2000_f19g16_CP_CTRL" then 
      TD1 = -alfas
   else 
      TD2 = -alfas
   end if 
;;;;
end do 

   TD_diff = TD2 - TD1
   copy_VarCoords(hfds_a_box(1,:,:), TD_diff)

   TD_diff = 86400*365*TD_diff
   


  wks = gsn_open_wks("pdf", "TD_diff_"+model_name(1))
 
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False

    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@pmLabelBarOrthogonalPosF = 0.3

    res@cnFillOn = True 
    res@cnLinesOn = False

    ; res@tmXBMode="Explicit"
    ; res@tmYLMode="Explicit"
    

    ; res@tmXBLabelFontHeightF      = 0.01    ;font height of tick labels
    ; res@tmYLLabelFontHeightF      = 0.01

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0
    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    ; res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "BlueRed"
    ; res@cnFillPalette = "BlueDarkRed18"
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4
    res@cnLevelSpacingF      =  0.1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    ; res@mpMinLatF             = -15.0
    ; res@mpMaxLatF             = 15.0
    ; res@mpMinLonF             = 125.0
    ; res@mpMaxLonF             = 280.0

    res@mpMinLatF             = -30.0
    res@mpMaxLatF             =30.0
    res@mpMinLonF             = 110.0
    res@mpMaxLonF             = 300.0

    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2

    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    ; res@gsnRightString = ev_sst@pcvar(0) +""
    res@gsnLeftStringFontHeightF = 0.02
    ; res@vpHeightF = 0.3
    ; res@vpWidthF = 0.7
    ; res@mpShapeMode = "FreeAspect"

    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=2.0
    res@tmYLMajorThicknessF=2.0
    res@tmBorderThicknessF=2.0

    res@gsnLeftString = "SST STD DIFF (WNP_heating - CON)"
    res@gsnRightString = ""
    plot = gsn_csm_contour_map(wks,TD_diff,res)
    ; res@gsnRightString = ev_sst@pcvar(1) +""
    ; plot(1) = gsn_csm_contour_map(wks,tos_enino2,res)

    draw(plot)
    frame(wks)
    ; pres = True
    ; gsn_panel(wks, plot, (/2,1/), pres)


end





