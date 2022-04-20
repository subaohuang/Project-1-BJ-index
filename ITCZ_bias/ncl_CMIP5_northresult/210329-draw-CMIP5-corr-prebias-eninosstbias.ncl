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


begin
;读入数据
    read_name = asciiread("/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/info.txt", -1,"string")
    print(read_name)
    model_name = read_name 

    dim_model = dimsizes(model_name)

    plot = new(12,graphic)
    iplot= plot

 
    f1 = addfile("all_month_pre_and_1000uv_bias.nc", "r")
    pre_bias = f1->pre_bias
 

    pre_bias = dim_rmvmean_n(pre_bias, 0)

    dim_pre = dimsizes(pre_bias)
   

    ; Nino3_bais_CMIP5 = asciiread("Nino3_bais_CMIP5", -1, "float") 

    Nino3_bais_CMIP5 = asciiread("Lanina3_bais_0.5std_CMIP5", -1, "float") 

    Nino3_bais_CMIP5_3d = conform_dims(dimsizes(pre_bias(:,1,:,:)), Nino3_bais_CMIP5, 0)

    ; wks = gsn_open_wks("pdf", "所有月份CMIP5模式降水bias和DJFnino3bias的相关系数")
    wks = gsn_open_wks("pdf", "所有月份CMIP5模式降水bias和DJF Lanino3bias的相关系数")

    month = ispan(1, 12, 1)

do i = 0, 11
   print("*******************"+month(i)+"**********************************************")
    corr = escorc_n(pre_bias(:,i,:,:), Nino3_bais_CMIP5_3d, 0, 0)
    copy_VarCoords(pre_bias(0,1,:,:), corr)

    prob = rtest(corr,dim_pre(0),0)
    copy_VarCoords(pre_bias(0,1,:,:), prob)

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn = False
    res@pmLabelBarOrthogonalPosF = 0.2
    res@cnFillOn = True 
    res@cnLinesOn = False


    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.01    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.01

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0



    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    res@gsnLeftStringFontHeightF = 0.03
    res@gsnLeftString = month(i)
       
    ; res@vpHeightF = 0.3
    ; res@vpWidthF = 0.7
    ; res@mpShapeMode = "FreeAspect"

    ; res@lbOrientation = "Vertical"

    ; corr = where(prob.le.0.3, corr, corr@_FillValue)
    plot(i) = gsn_csm_contour_map(wks,corr,res)


; ;;;打点
;   res_pro=True
;   res_pro@gsnDraw=False
;   res_pro@gsnFrame=False
;   res_pro@cnLevelSelectionMode = "ExplicitLevels"
;   res_pro@cnLinesOn=False
;   res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
;   res_pro@cnMinLevelValF      = 0.00        ; set min contour level
;   res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
;   res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing
;   ; res_pro@cnLevels = ispan(-6,6,1)
;   res_pro@cnLineLabelsOn=False
;   res_pro@tmYRMode = "Automatic"
;   res_pro@cnInfoLabelOn=False
;   ; res_pro@cnLineLabelBackgroundColor="Transparent"
;   ; res_pro@cnLineThicknessF = 3.
;   ; res_pro@cnFillOn        = True
 
;   iplot(i) = gsn_csm_contour(wks,prob,res_pro)

;   opt                      = True
;   opt@gsnShadeFillType     = "pattern"
;   opt@gsnShadeMid     = 17         ; stipple
;   ; opt@gsnShadeHigh   = "green"
;   opt@gsnShadeFillScaleF   = 0.2         ; add extra density
;   opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger
;   iplot(i) = gsn_contour_shade(iplot(i), 0,0.1, opt)
; ;;;  

; overlay(plot(i), iplot(i))
end do 

pres = True
pres@gsnPanelLabelBar  =True
pres@pmLabelBarOrthogonalPosF = - 0.1

gsn_panel(wks, plot, (/4,3/), pres)


end 


