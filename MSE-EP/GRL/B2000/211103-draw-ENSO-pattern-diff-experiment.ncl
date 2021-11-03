load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =25601
    endmon = 034012

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
    ; fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h1.0251-0280.nc"
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
    ; sst = f1->TEMP(timeind,0,:,:)  
   ;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})

    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)

;;;读入B22000 WNPCLM海温数据
    modelname = "lsh_B2000_alt_north_db_year_80_20"
    ; ;;;;
    ; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    ; fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"+ modelname+"/"+modelname+".TEMP.h1.0251-0280.nc"
    ; fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    ; fils1 ="/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    print(date(0))

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst_r := f1->TEMP(timeind,0,:,:)  
    ; sst := f1->TEMP(timeind,0,:,:)  
    ; ;;;;差值TEMP 
    sst := PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    sst_wnp = sst(:,{lat1:lat2},{lon1:lon2})
    sst_wnp = (/rmMonAnnCycTLL(sst_wnp)/)
;;;;


;;;;;
    sst_ct_djf = month_to_season(sst_ct, "DJF")
    sst_ct_djf(0:dimsizes(sst_ct_djf&time)-2,:,:) = sst_ct_djf(1:dimsizes(sst_ct_djf&time)-1,:,:)
    sst_ct_djf(dimsizes(sst_ct_djf&time)-1,:,:) = sst_ct(dimsizes(sst_ct&time)-1,:,:)

    
    sst_wnp_djf = month_to_season(sst_wnp, "DJF")
    sst_wnp_djf(0:dimsizes(sst_wnp_djf&time)-2,:,:) = sst_wnp_djf(1:dimsizes(sst_wnp_djf&time)-1,:,:)
    sst_wnp_djf(dimsizes(sst_wnp_djf&time)-1,:,:) = sst_wnp(dimsizes(sst_wnp&time)-1,:,:)


    nino3_4_ct =  dim_avg_n_Wrap(sst_ct_djf(:,{-5:5},{190:240}), (/1,2/))   
    nino3_4_wnp =  dim_avg_n_Wrap(sst_wnp_djf(:,{-5:5},{190:240}), (/1,2/))   

    nino4_ct =  dim_avg_n_Wrap(sst_ct_djf(:,{-5:5},{150:200}), (/1,2/))   
    nino4_wnp =  dim_avg_n_Wrap(sst_wnp_djf(:,{-5:5},{150:200}), (/1,2/))  

    nino3_ct =  dim_avg_n_Wrap(sst_ct_djf(:,{-5:5},{210:270}), (/1,2/))   
    nino3_wnp =  dim_avg_n_Wrap(sst_wnp_djf(:,{-5:5},{210:270}), (/1,2/))  
    
    
    std_ct_nino3 = dim_stddev_n(nino3_ct, 0)
    std_wnp_nino3 = dim_stddev_n(nino3_wnp, 0)
    
    std_ct_nino3_4 = dim_stddev_n(nino3_4_ct, 0)
    std_wnp_nino3_4 = dim_stddev_n(nino3_4_wnp, 0)

    std_ct_nino4 = dim_stddev_n(nino4_ct, 0)
    std_wnp_nino4 = dim_stddev_n(nino4_wnp, 0)

    cp_enino_year_ct = ind(nino4_ct.ge.std_ct_nino4)
    ep_enino_year_ct = ind(nino3_ct.ge.std_ct_nino3)

    cp_enino_year_wnp = ind(nino4_wnp.ge.std_wnp_nino4)
    ep_enino_year_wnp = ind(nino3_wnp.ge.std_wnp_nino3)

    enino_year_ct = ind(nino3_4_ct.ge.std_ct_nino3_4)
    enino_year_wnp = ind(nino3_4_wnp.ge.std_wnp_nino3_4)

    ep_enino = ind(nino3_4_ct.ge.std_ct_nino3_4)

    ep_year_ct = ind(nino3_4_ct.ge.std_ct_nino3_4.and.(nino3_ct.ge.nino4_ct))

    print(ep_year_ct)
    print(ep_enino)
    exit()

    cp_year_ct = ind(nino3_4_ct.ge.std_ct_nino3_4.and.(nino3_ct.lt.nino4_ct))


    ep_year_wnp = ind(nino3_4_wnp.ge.std_wnp_nino3_4.and.(nino3_wnp.ge.nino4_wnp))

    cp_year_wnp = ind(nino3_4_wnp.ge.std_wnp_nino3_4.and.(nino3_wnp.lt.nino4_wnp))

    print(ep_year_ct)
    print(cp_year_ct)


    print(ep_year_wnp)
    print(cp_year_wnp)
    exit()
    

  
  
    sst_ct_djf_ep_enino = dim_avg_n_Wrap(sst_ct_djf(ep_enino_year_ct,:,:), 0)
    sst_wnp_djf_ep_enino = dim_avg_n_Wrap(sst_wnp_djf(ep_enino_year_wnp,:,:), 0)

    
    sst_ct_djf_cp_enino = dim_avg_n_Wrap(sst_ct_djf(cp_enino_year_ct,:,:), 0)
    sst_wnp_djf_cp_enino = dim_avg_n_Wrap(sst_wnp_djf(cp_enino_year_wnp,:,:), 0)
    
    
    data = new((/3,3/),float)
    data(0,:) = (/1,1.52,1.46/)
    data(1,:) = (/0.6579,1.00,1.55/)
    data(2,:) = (/0.6849,0.6452,1.00/)

    perc1 = sum(data(1,:))/sum(data)
    print(perc1)
    exit()

    ; sst_diff = sst_ct_djf_enino
    ; sst_diff = (/sst_wnp_djf_enino - sst_ct_djf_enino/)

    exit()
;;;




;;;;
    plot = new(3,graphic)

    wks = gsn_open_wks("pdf", "NINO4-Elnino_spatial_pattern_"+modelname)

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False

    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@pmLabelBarOrthogonalPosF = 0.3

    res@cnFillOn = True 
    res@cnLinesOn = False

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0
    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -2.0
    res@cnMaxLevelValF       =   2.0
    res@cnLevelSpacingF      =  0.4
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
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

    res@gsnLeftString = "El Nino SSTA in control run"
    res@gsnRightString = ""
    plot(0) = gsn_csm_contour_map(wks,sst_ct_djf_enino,res)
    res@gsnLeftString = "El Nino SSTA in ATL warming run"
    plot(1) = gsn_csm_contour_map(wks,sst_wnp_djf_enino,res)
    res@gsnLeftString = "ATL warming - CTRL"
    plot(2) = gsn_csm_contour_map(wks,sst_diff,res)
    ; res@gsnRightString = ev_sst@pcvar(1) +""
    ; plot(1) = gsn_csm_contour_map(wks,tos_enino2,res)

    pres = True
    gsn_panel(wks, plot, (/3,1/), pres)


end





