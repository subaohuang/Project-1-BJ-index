load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =25101
    endmon = 28012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -90
   lat2 = 90
   lon1 = 0
   lon2 = 360
;;;;

;;;;读入控制实验结果
    ; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    ; fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    ; sst_r = f1->TEMP(timeind,0,:,:)  
    sst = f1->TEMP(timeind,0,:,:)  
   ;;;差值TEMP 
    ; sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    ; printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})

    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)

;;;读入B22000 WNPCLM海温数据

;;;;model name 
 ; modelname = "B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer"
 ; modelname = "lsh_B2000_WNP_heating_Idl"
; modelname = "lsh_B2000_WNP_heating_IDLthree"
   ; modelname = "lsh_B2000_WNP_heating_double"
     ; modelname = "lsh_B2000_WNP_heating_db"
 modelname = "HEAT_B_BP"
; ;;;;
; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
; fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")

; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"
; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"+ modelname+"/"+modelname+".TEMP.fvgrid.0251-0280.nc"
    ; fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils1 ="/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.TEMP.fvgrid.0251-0280.nc"

    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    print(date(0))

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    ; sst_r := f1->TEMP(timeind,:,:)  
    sst := f1->TEMP(timeind,0,:,:)  
   ; ;;;;差值TEMP 
   ;  sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
   ;  printVarSummary(sst)

    sst_wnp = sst(:,{lat1:lat2},{lon1:lon2})
    sst_wnp = (/rmMonAnnCycTLL(sst_wnp)/)


;;;;
   nino3_ct =  dim_avg_n_Wrap(sst_ct(:,{-5:5},{210:270}), (/1,2/))   
   ; nino3_sen =  dim_avg_n_Wrap(sst_sen(:,{-5:5},{210:270}), (/1,2/))   
   nino3_wnp =  dim_avg_n_Wrap(sst_wnp(:,{-5:5},{210:270}), (/1,2/))   

   std1_nino3 = dim_stddev_n(nino3_ct, 0)
   std2_nino3 = dim_stddev_n(nino3_wnp, 0)

   ; std1  = dim_stddev_n_Wrap(sst_ct, 0)
   ; std3  = dim_stddev_n_Wrap(sst_wnp, 0)

   ; var1 = dim_variance_n_Wrap(sst_ct, 0)
   ; var2 = dim_variance_n_Wrap(sst_wnp, 0)

    ; sst_ct_jja = month_to_season(sst_ct, "DJF")
    ; sst_wnp_jja = month_to_season(sst_wnp, "DJF")

    ; std1  = dim_stddev_n_Wrap(sst_ct_jja, 0)
    ; std3  = dim_stddev_n_Wrap(sst_wnp_jja, 0)

    sst_ct_jja = month_to_season(sst_ct, "DJF")
    sst_wnp_jja = month_to_season(sst_wnp, "DJF")

    var1 = dim_variance_n_Wrap(sst_ct_jja, 0)
    var2 = dim_variance_n_Wrap(sst_wnp_jja, 0)

    std1  = dim_stddev_n_Wrap(sst_ct_jja, 0)
    std3  = dim_stddev_n_Wrap(sst_wnp_jja, 0)


   s1 = var1
   s1 = dimsizes(sst_ct(:,1,1))
   s2 = s1
   prob = ftest(var1, s1, var2, s2, 0)
   
   
   std_diff = std1
   std_diff = (/std3 - std1/)

   std_diff = where(prob.le.0.1, std_diff, std_diff@_FillValue)
   

;;;

  wks = gsn_open_wks("pdf", "Nino3_std_diff_"+modelname)

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
    plot = gsn_csm_contour_map(wks,std_diff,res)
    ; res@gsnRightString = ev_sst@pcvar(1) +""
    ; plot(1) = gsn_csm_contour_map(wks,tos_enino2,res)

    draw(plot)
    frame(wks)
    ; pres = True
    ; gsn_panel(wks, plot, (/2,1/), pres)


end





