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

    plot = new(dim_model,graphic)
    iplot= plot
    plot1 = plot

    yearstart = 195001
    yearend = 200412

    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)

    
    startyear_r = str_split_by_length(yearstart, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(yearend, 4 )

    yearnum = stringtoint(endyear_r(0)) - stringtoint(startyear_r(0)) + 1
    print(yearnum)
;;;


    wks = gsn_open_wks("pdf", "30个CMIP5模式nino-pattern-EOF")

    lat1 = -15
    lat2 = 15
    lon1 = 125
    lon2 = 280


do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
 
    ; fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/EOF/CP/CP-EOF-ssta-djf-"+model+".nc")
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/EOF/EOF-ssta-djf-"+model+".nc")

    f2 := addfile(fils2, "r")

    pattern_model = f2->pattern_model
    pcvar = f2->pcvar
    mean  = avg(pattern_model({-10:10},{180:})) ;;CP
    ; mean  = avg(pattern_model({-10:10},{160:210})) ;;CP
    ; mean  = avg(pattern_model({-10:10},{150:270})) ;;EP


    if mean.le.0 then 
     pattern_model = (/-1*pattern_model/)
    end if
    tos_enino = pattern_model

    enino_intensity(k) = max(tos_enino({-10:10},{120:280})) 

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn = False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025

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
    res@cnMinLevelValF       =   -1
    res@cnMaxLevelValF       =   1
    res@cnLevelSpacingF      =  0.2
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -15.0
    res@mpMaxLatF             = 15.0
    res@mpMinLonF             = 125.0
    res@mpMaxLonF             = 280.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    res@gsnRightString = pcvar + "%"
    res@gsnLeftString = model 

    res@gsnLeftStringFontHeightF = 0.03
       
    res@vpHeightF = 0.3
    res@vpWidthF = 0.7
      res@mpShapeMode = "FreeAspect"
    plot(k) = gsn_csm_contour_map(wks,tos_enino,res)



 print("********************************finish calculate"+" "+model+"************************************")
end do
  asciiwrite("Enino_intensity_EOF1_30model", enino_intensity)
   
    pres = True
    pres@gsnMaximize = True
    pres@gsnPanelLabelBar = True
    pres@pmLabelBarHeightF = 0.05
    pres@pmLabelBarWidthF  = 0.6
    pres@pmLabelBarOrthogonalPosF = -0.02
    ; pres@lbOrientation = "Vertical"
    pres@lbLabelFontHeightF = 0.01
    gsn_panel(wks, plot, (/6,5/), pres)
end 


