load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =25601
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
 
; modelname = "lsh_B2000_WNP_heating_Idl"
modelname = "lsh_B2000_WNP_heating_IDLthree"

;;;;读入控制实验结果
    fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CP_CTRL/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)
    timeind=ind(date.le.endmon.and.date.ge.startmon)
    sst = f1->TEMP(timeind,0,:,:)  
   ;;;差值TEMP 
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})
    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)
    
;;;读入B22000 WNPCLM海温数据
  
    ; fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.TEMP.fvgrid.0251-0280.nc"
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"+modelname+"/"+modelname+".TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)


    timeind:=ind(date.le.endmon.and.date.ge.startmon)
    sst := f1->TEMP(timeind,0,:,:)  

    sst_wnp = sst(:,{lat1:lat2},{lon1:lon2})
    sst_wnp = (/rmMonAnnCycTLL(sst_wnp)/)
;;;;


;;;读入OLR数据
  fils2 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CP_CTRL/B2000_f19g16_CP_CTRL.cam.h1.0251-0280_OLR_heatflux.nc"
   
   f2 = addfile(fils2, "r")
   time := f2->time 
   time = (/time-1/)
   date:=cd_calendar(time, 1)
   
   timeind:=ind(date.le.endmon.and.date.ge.startmon)
   olr_ctrl = f2->FLUT(timeind,:,:)  
   olr_ctrl = (/rmMonAnnCycTLL(olr_ctrl)/)


   ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.cam.h1.0251-0280_OLR_heatflux.nc"
   fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/"+modelname+"/"+modelname+".cam.h1.0251-0280_OLR_heatflux.nc"
   f2 = addfile(fils2, "r")
   time := f2->time 
   time = (/time-1/)
   date:=cd_calendar(time, 1)

   timeind:=ind(date.le.endmon.and.date.ge.startmon)
   olr_wnp = f2->FLUT(timeind,:,:)  
   olr_wnp = (/rmMonAnnCycTLL(olr_wnp)/)
;;;


;;;nino3 
   
   nino3_ctrl = dim_avg_n_Wrap(sst_ct(:,{-5:5},{210:270}), (/1,2/))
   nino3_wnp = dim_avg_n_Wrap(sst_wnp(:,{-5:5},{210:270}), (/1,2/))
;;;

;;;regression 单点回归
   
   ; rc_ctrl = regCoef_n(nino3_ctrl, olr_ctrl, 0, 0) 
   ; rc_wnp = regCoef_n(nino3_wnp, olr_wnp, 0, 0) 
   

    rc_ctrl = escorc_n(sst_ct, olr_ctrl, 0, 0) 
    rc_wnp = escorc_n(sst_wnp, olr_wnp, 0, 0) 

    var1 = dim_stddev_n(sst_ct, 0)
    var2 = dim_stddev_n(olr_ctrl, 0)

    rc_ctrl = (/-1*rc_ctrl*var2/var1/)


    var3 = dim_stddev_n(sst_wnp, 0)
    var4 = dim_stddev_n(olr_wnp, 0)

    rc_wnp = (/-1*rc_wnp*var4/var3/)

    copy_VarCoords(olr_ctrl(0,:,:), rc_ctrl)
    copy_VarCoords(olr_ctrl(0,:,:), rc_wnp)
   

   rc_diff = rc_ctrl
   rc_diff = (/rc_wnp - rc_ctrl/)
;;;

  wks = gsn_open_wks("pdf", "OLR-regression-onto-nino3-CESM-"+modelname)

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

    res@tmXBLabelFontHeightF      = 0.01    ;font height of tick labels
    res@tmYLLabelFontHeightF      = 0.01

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0

    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    ; res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "GMT_hot"
    res@cnFillPalette = "cmocean_matter"
    ; res@cnFillPalette = "BlueDarkRed18"
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   0.0
    res@cnMaxLevelValF       =   30.0
    res@cnLevelSpacingF      =  4.0

    ; res@cnMinLevelValF       =   -1.0
    ; res@cnMaxLevelValF       =   1.0
    ; res@cnLevelSpacingF      =  0.1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = -30.0
    res@mpMaxLatF             =30.0
    res@mpMinLonF             = 110.0
    res@mpMaxLonF             = 280.0

    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2

    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    res@gsnLeftStringFontHeightF = 0.03
    ; res@vpHeightF = 0.3
    ; res@vpWidthF = 0.7
    ; res@mpShapeMode = "FreeAspect"
    plot = new(3,graphic)
    res@gsnLeftString = "WNP heating"
    plot(0) = gsn_csm_contour_map(wks,rc_wnp,res)
    res@gsnLeftString = "CTRL"
    plot(1) = gsn_csm_contour_map(wks,rc_ctrl,res)
     res@cnMinLevelValF       =   -10.0
    res@cnMaxLevelValF       =   10.0
    res@cnLevelSpacingF      =  2.0
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnFillPalette = "BlueRed"
    res@gsnLeftString = "WNP heating - CTRL"
    plot(2) = gsn_csm_contour_map(wks,rc_diff,res)


    ; res@gsnRightString = ev_sst@pcvar(1) +""
    ; plot(1) = gsn_csm_contour_map(wks,tos_enino2,res)

    pres = True
    pres@gsnPanelLabelBar = True
    gsn_panel(wks, plot, (/3,1/), pres)


end





