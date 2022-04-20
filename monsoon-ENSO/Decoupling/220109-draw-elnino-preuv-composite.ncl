;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
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

function lead_corr(fils1,fils2,Leftstring, wks)
begin

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

   
    time1 := f1->time
    time1 = (/time1 -1/)

    date1 := cd_calendar(time1, 1)
    timeind1 := ind(date1.ge.25601.and.date1.le.30012)

    PRECL_CT  =  f1->PRECL(timeind1,:,:)
    PRECC_CT  =  f1->PRECC(timeind1,:,:)

    PRE_CT = PRECC_CT
    PRE_CT = (/PRECC_CT+PRECL_CT/)
    PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day
     
    ; slp = f1->PSL(timeind1,:,:)
    
    slp = PRE_CT
   
    ;;;;U,V
    v_CT_850  =  f1->V(timeind1,{850},:,:)
    u_CT_850  =  f1->U(timeind1,{850},:,:)

    u_CT_850 = (/rmMonAnnCycTLL(u_CT_850)/)
    v_CT_850 = (/rmMonAnnCycTLL(v_CT_850)/)


    slp  = (/rmMonAnnCycTLL(slp)/)    

    


;;;;
    time2 = f2->time
    time2 = (/time2 -1/)
    date2=cd_calendar(time2, 1)
    timeind2 = ind(date2.ge.25601.and.date2.le.30012)


    if fils2.eq."/home/yangsong3/CMIP6/linshh/CESM-data/B/Decouping/B2000_decoupling_indian_ocean/POP/B2000_decoupling_indian_ocean.sst.0251-0300_fvgrid.nc"  then
    sst_ctrl = f2->TEMP(timeind2,:,:)
    else   
    sst_ctrl = f2->TEMP(timeind2,0,:,:)
    end if
;;;

   sst_ctrl = (/rmMonAnnCycTLL(sst_ctrl)/)
;;;;Nino3.4
    Nino3_ctrl = dim_avg_n_Wrap(sst_ctrl(:,{-5:5},{190:240}), (/1,2/))

    dim = dimsizes(Nino3_ctrl)
    print(dim)

    Nino3_ctrl_DJF = month_to_season(Nino3_ctrl, "DJF") ;; 25601-25602, 25612-25701-02

    std = stddev(Nino3_ctrl_DJF)

    elnino_y = ind(Nino3_ctrl_DJF.ge.0.5*std)
    printVarSummary(Nino3_ctrl_DJF)
    elnino = elnino_y(0:dimsizes(elnino_y)-2)
;;;;



    ;;;season
    u_jja = month_to_season(u_CT_850,"JJA")
    v_jja = month_to_season(v_CT_850,"JJA")
    sst_jja = month_to_season(sst_ctrl,"JJA")
    slp_jja = month_to_season(slp,"JJA")


;;;;Nino year 
    
    u_elnino = dim_avg_n_Wrap(u_jja(elnino+1,:,:), 0)
    v_elnino = dim_avg_n_Wrap(v_jja(elnino+1,:,:), 0)

    sst_elnino = dim_avg_n_Wrap(sst_jja(elnino+1,:,:), 0)

    slp_elnino = dim_avg_n_Wrap(slp_jja(elnino+1,:,:), 0)
;;;;


;;;显著性检验
  prob_u = one_test(u_jja,elnino)
  prob_v = one_test(v_jja,elnino)
  prob_sst = one_test(sst_jja,elnino)
  prob_slp = one_test(slp_jja,elnino)
;;;

    plot = new(2,graphic)
    iplot = new(2,graphic)


printMinMax(slp_elnino, 1)

;;;
;;**************************************************************
;画图
;;***************************************************************
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    ; res@cnFillPalette="MPL_BrBG"
    res@cnFillOn=True
    res@cnLinesOn=False
    res@gsnAddCyclic = False
    res@cnLineLabelsOn=False


    res@lbLabelBarOn= True

    res@pmLabelBarOrthogonalPosF = 0.2

    res@lbLabelStride=1
    res@lbOrientation="horizontal"
    res@mpLimitMode           = "LatLon"


    ; res@mpMinLatF             = -40
    ; res@mpMaxLatF             = 40
    ; res@mpMinLonF             = 0
    ; res@mpMaxLonF             = 360

    res@mpMinLatF             = -50
    res@mpMaxLatF             = 50
    res@mpMinLonF             = 0
    res@mpMaxLonF             = 360


    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmBorderThicknessF=1.0

    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF = 1.5



  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 
  
  ; res@tmYLValues = (/-30,-20,-10,0,10,20,30/) 
  ; res@tmYLLabels = (/"30S","20S","10S","0","10N","20N","30N"/) 

  ; res@tmYLMinorValues = ispan(-30,30,5)
  res@gsnScalarContour        = True
  res@mpFillOn = False

  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "CBR_drywet"
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   2.0
  res@cnLevelSpacingF      =  0.4

  ; res@cnFillPalette = "precip_diff_12lev"
  ; ; res@cnFillPalette = "CBR_drywet"
  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ;   res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ;   res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5/)

  res@lbLabelStride=1


  res@cnExplicitLabelBarLabelsOn = True
  ; res@lbLabelStrings =  res@cnLevels
  res@vcExplicitLabelBarLabelsOn = True
  res@stExplicitLabelBarLabelsOn = True
  res@lbLabelFontHeightF = 0.016
  res@lbAutoManage = False

  res@vcGlyphStyle="LineArrow"
  res@vcMonoLineArrowColor = True
  res@vcLineArrowColor = "black"
  res@vcMinDistanceF           = 0.015          ; thin vectors

  res@vcRefLengthF             = 0.05    ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcVectorDrawOrder="PostDraw"
  
  ; resv@vcMonoLineArrowColor = True
  ; resv@vcLineArrowColor = "blue"
  res@vcLineArrowThicknessF=1.2
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

  res@vcRefMagnitudeF          = 3.0
  res@vcMinMagnitudeF          =.3
  ; res@gsnLeftStringFontHeightF = 0.0185
  res@gsnLeftStringFontHeightF = 0.02

  ; res@vcRefAnnoString1On = True
  ; res@vcRefAnnoString1 = "2 m/s"
  ; res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
  ; res@vcRefAnnoParallelPosF = 1.02
  ; res@vcRefAnnoFontHeightF = 0.013
  ; res@vcRefAnnoBackgroundColor = "Transparent"
  ; res@vcRefAnnoPerimOn = False

    res@gsnLeftString = Leftstring
    ; res@gsnLeftStringOrthogonalPosF  = -0.12

  
  plot(0) = gsn_csm_vector_scalar_map(wks, u_elnino ,v_elnino , slp_elnino , res)

  
  ; delete(res@cnLevelSelectionMode)
  ; delete(res@cnFillPalette)
  ; delete(res@cnLevels)
  ; delete(res@cnFillColors)
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "MPL_bwr"

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; res@lbLabelStride=1

  ; res@cnLevels := (/-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0/)*10

   ; res@cnLevels := (/-16,-12,-8,-4,-2,2,4,8,12,16/)

  res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res@cnMinLevelValF      = -0.6       ; set min contour level
  res@cnMaxLevelValF      = 0.6      ; set max contour level
  res@cnLevelSpacingF     = 0.1       ; set contour spacing

  plot(1) = gsn_csm_contour_map(wks, sst_elnino, res)


    ;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(0) = gsn_csm_contour(wks,prob_slp,res_pro)

    iplot(1) = gsn_csm_contour(wks,prob_sst,res_pro)
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    iplot(0) = gsn_contour_shade(iplot(0) , 0,0.1,opt)
    iplot(1) = gsn_contour_shade(iplot(1) , 0,0.1,opt)
        
    
    ; overlay(plot(0), iplot(0))
    ; overlay(plot(1), iplot(1))
 
   return plot
end 

begin
    ;;;cam4 输出结果
    fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/Cam/B2000_f19g16_CTRL.cam.h1.0250-0300_chazhi.nc"
    ;;;pop 输出结果
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CTRL.sst.0251-0300_fvgrid.nc"



    fils3 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/Decouping/B2000_decoupling_indian_ocean/B2000_decoupling_indian_ocean.cam.h1.0251-0300_chazhi.nc"
    fils4 ="/home/yangsong3/CMIP6/linshh/CESM-data/B/Decouping/B2000_decoupling_indian_ocean/POP/B2000_decoupling_indian_ocean.sst.0251-0300_fvgrid.nc" 
    

    ; fils5 =  "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CP_CTRL.cam.h1.0251-0340_chazhi.nc" 
    ; fils6 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_f19g16_CTRL/B2000_f19g16_CP_CTRL.TEMP.h1.0251-0340.nc"  
    

    ; fils7 =  "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer.cam.h0.0251-0280_chazhi.nc" 
    ; fils8 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer.TEMP_taux_tauy.h0.0251-0280.nc"  

    ; wks = gsn_open_wks("pdf", "elnino-decay-year-JJA-wind-slp-sst-ATrun")
    wks = gsn_open_wks("pdf", "PRE-elnino-decay-year-JJA-Indian-Decou")

    ; Leftstring = (/"ATFIX","ATFIX_warming","FREE","WNP&ATFIX"/)
    Leftstring = (/"All ocean Free","Indian Ocean decouping"/)
    corr_ATLyear_ctrl = lead_corr(fils1,fils2,Leftstring(0),wks)

    corr_ATLyear = lead_corr(fils3,fils4,Leftstring(1),wks)

    ; corr_ATLfree = lead_corr(fils5,fils6,Leftstring(2),wks)

    ; corr_ATLwnpclm = lead_corr(fils7,fils8,Leftstring(3),wks)


    plot_total = new(4,graphic)

    plot_total(0) = corr_ATLyear_ctrl(0)
    plot_total(1) = corr_ATLyear_ctrl(1)
    plot_total(2) = corr_ATLyear(0)
    plot_total(3) = corr_ATLyear(1)
    ; plot_total(4) = corr_ATLfree(0)
    ; plot_total(5) = corr_ATLfree(1)
    ; plot_total(6) = corr_ATLwnpclm(0)
    ; plot_total(7) = corr_ATLwnpclm(1)
  
  pres  = True
  pres@gsnPanelLabelBar = False

  gsn_panel(wks,(/plot_total(0),plot_total(2)/), (/2,1/), pres)

end  