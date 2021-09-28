;;********************************************************
;;画散点图,得到关系
;;
;;lsh 2020-11-04
;**********************************************************
begin
;读入数据

    read_name = asciiread("/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/info.txt", -1,"string")
    print(read_name)
    model_name = read_name(0)

    dim_model = dimsizes(model_name)

    ; plot = new(dim_model,graphic)
    ; iplot= plot
    ; plot1 = plot

    ;;;read data
    lat1 = -5
    lat2 =  30
    lon1 =  100
    lon2 = 170

    ;;;
    yearstart = 185001
    yearend = 201412
    dbitcz = new(dim_model,float)
    enino_intensity = new(dim_model,float)
;;;

   plot = new(2, graphic)
   

do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
   
    fils1 =  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/pr_chazhi/pr_Amon_"+model+"_historical_r1i1p1_"+"*")

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    ind1 = ind(date.le.190012.and.date.ge.185001)
    ind2 = ind(date.le.200412.and.date.ge.195001)
    delete(date)     
    pr_r = lonFlip(f1->pr(timeind,:,:))
    pr = pr_r(:,{lat1:lat2},{lon1:lon2})
    delete(pr_r)
   
    pr = (/pr*86400/) ;;conver to mm/day
    printVarSummary(pr)
    printMinMax(pr, 1)

    pr_clm1 = dim_avg_n_Wrap(month_to_season(pr(ind1,:,:),"JJA"), 0)
    pr_clm2 = dim_avg_n_Wrap(month_to_season(pr(ind2,:,:),"JJA"), 0)

    pr_diff = pr_clm1

    pr_diff = (/pr_clm2 - pr_clm1/)


    print("*******************start calculate"+" "+model+"**********************************************")
 
    fils2 :=  systemfunc("ls /home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/"+model+"/tos_Omon_"+model+"_historical_r1i1p1_"+"*")

    f2 := addfiles(fils2, "r")

    time:=f2[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    ; timeind:=ind(date.le.200412.and.date.ge.195001)
       ind1 = ind(date.le.190012.and.date.ge.185001)
        ind2 = ind(date.le.200412.and.date.ge.195001)
    delete(date)     
    tos_r = lonFlip(f2[:]->tos(timeind,:,:))
    tos = tos_r(:,{lat1:lat2},{lon1:lon2})
    printVarSummary(tos)

    tos = (/tos - 273.6/)

    

    tos_clm1 = dim_avg_n_Wrap( month_to_season(tos(ind1,:,:),"JJA"),0)
    tos_clm2 = dim_avg_n_Wrap( month_to_season(tos(ind2,:,:),"JJA"),0)

    tos_diff = tos_clm1

    tos_diff = (/tos_clm2 - tos_clm1/)



   wks=gsn_open_wks("pdf", "SST-pre-diff-"+model)

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False


 
    res@lbLabelBarOn=True

    res@cnFillOn = True 
    res@cnLinesOn = False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF      = 0.025    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.025

    ; res@gsnContourNegLineDashPattern = 1
    ; res@gsnContourPosLineDashPattern = 0
    ; res@gsnContourZeroLineThicknessF = 2.0

    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "BlueDarkRed18"



    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = lat1
    res@mpMaxLatF             = lat2
    res@mpMinLonF             = lon1
    res@mpMaxLonF             = lon2
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    res@gsnRightString = ""
    res@gsnLeftString = model
    res@gsnLeftStringFontHeightF = 0.03
    plot(0) = gsn_csm_contour_map(wks, pr_diff,res)

    res@cnFillPalette = "BlueWhiteOrangeRed"
        res@cnMinLevelValF       =   -0.5
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.1

    plot(1) = gsn_csm_contour_map(wks, tos_diff,res)



    pres  = True
    gsn_panel(wks, plot, (/2,1/), pres)

end do 
end 



