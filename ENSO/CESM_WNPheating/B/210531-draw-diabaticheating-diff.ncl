load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
;;**************************************************************
;;;main program
;;***************************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      ; simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 
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


;;;读入控制实验组的海温数据
    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils1 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_ctrl = f1->PTENDT(timeind,{1000:100},:,:)
    sst_ctrl_avg = dim_avg_n_Wrap(sst_ctrl, 1)

    sst_ctrl_spring = month_to_season(sst_ctrl, "MAM")


    ; sst_ctrl_DTCOND = f1->DTCOND(timeind,{1000:100},:,:)  ;;;k/s
    ; sst_ctrl_DTCOND_avg = dim_avg_n_Wrap(sst_ctrl_DTCOND, 1)
    ; sst_ctrl_DTCOND_mam = month_to_season(sst_ctrl_DTCOND_avg, "MAM")
    ; sst_ctrl_DTCOND_clm = dim_avg_n_Wrap(sst_ctrl_DTCOND_mam, 0)
 
;;;读入B22000 WNPCLM海温数据
    
    ; fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/lsh_B2000_WNP_heating_Idl.PTENDT.cam.h1.0251-0280_chazhi.nc"
    fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/lsh_B2000_WNP_heating_double.PTENDT.cam.h1.0251-0280_chazhi.nc"
    ; fils2 = "/home/ys17-19/cjx/all_output_data/HEAT_B_BP_output/HEAT_B_BP.sst.250-289.nc"
    
    f2 = addfile(fils2, "r")
    time := f2->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)
    

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    sst_wnp = f2->PTENDT(timeind,{1000:100},:,:)  
    ; sst_wnp = f2->DTCOND(timeind,{1000:100},:,:)  
    sst_wnp_avg = dim_avg_n_Wrap(sst_wnp, 1)
      
;;;;
    sst_ctrl_mam = month_to_season(sst_ctrl_avg, "MAM")
    sst_wnp_mam = month_to_season(sst_wnp_avg, "MAM")

    sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl_mam, 0)
    sst_wnp_clm = dim_avg_n_Wrap(sst_wnp_mam, 0)

    ; sst_ctrl_clm = dim_avg_n_Wrap(sst_ctrl, 0)
    ; sst_wnp_clm = dim_avg_n_Wrap(sst_wnp, 0)
    
    SST_diff = sst_wnp_clm
    SST_diff = (/sst_wnp_clm - sst_ctrl_clm/)
    ; SST_diff = (/sst_ctrl_clm/)
    printVarSummary(SST_diff)

    SST_diff = (/86400*SST_diff/1004.0/)

;;; get variables: temperature, longitude, latitude


; heat_coff = SST_diff
; heat_coff = 0.0
; latS = 0.0
; latN = 20.0
; lonW = 125.0
; lonE = 160.0

; lon = SST_diff&lon
; lat = SST_diff&lat
; latC = (latS+latN)*0.5
; lonC = (lonW+lonE)*0.5
; ;modified heat

; do i=0,dimsizes(lon)-1,1
; do j=0,dimsizes(lat)-1,1
; if (lon(i).ge.lonW .and. lon(i).le.lonE .and. lat(j).ge.latS .and. lat(j).le.latN) then
; stat_lat=lat(j)
; stat_lon=lon(i)
; heat_coff(j,i) = tofloat((1 - abs(stat_lat-latC)/abs(latN-latC))*(1 - abs(stat_lon-lonC)/abs(lonE-lonC))) 
; end if
; end do
; end do


; heat_mar =(/1.00 ,1.00 \
;     ,1.00 ,1.15 ,1.28 ,1.29 ,1.26 ,1.24 ,1.23 ,1.21 ,1.20 ,1.20 ,1.17 ,1.08 \
;     /)

; heat_coff_3d = conform_dims(dimsizes(sst_ctrl_spring), heat_coff, (/2,3/))
; heat_mar_3d = conform_dims(dimsizes(sst_ctrl_spring), heat_mar, (/1/))


; heat_add =sst_ctrl_spring
; heat_add = sst_ctrl_spring + ((1004.0*2.0)/86400.0)*heat_coff_3d*heat_mar_3d

;   print(avg(heat_add(:,:,{0:20},{125:160})))
;   print(avg(sst_ctrl_spring(:,:,{0:20},{125:160})))
; exit()
    ; prob_p = t_test(sst_ctrl_mam ,sst_wnp_mam)
; SST_diff = where(prob_p.le.0.1, SST_diff, SST_diff@_FillValue)
;;;
   ; SST_diff = (/2*SST_diff*86400/)


   wks = gsn_open_wks("pdf", "diabaticheating-diff-WNPdbheating-minus-control")
     ; wks = gsn_open_wks("pdf", "DTCOND-diabaticheating-control")

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True
     res@trGridType = "TriangularMesh"

    ; res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnFillPalette = "BlRe"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -2.0
    res@cnMaxLevelValF       =   2.0
    res@cnLevelSpacingF      =  0.5
    ; ;  res@cnMinLevelValF       =   -0.02
    ; res@cnMaxLevelValF       =   0.02
    ; res@cnLevelSpacingF      =  0.005
    
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -50.0
    res@mpMaxLatF             = 50.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2.0
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3
    ; ;;
    res@gsnRightString = ""

    plot= gsn_csm_contour_map(wks,SST_diff,res)

    ; plot(i) = gsn_csm_contour_map(wks,PRE_diff,res)


    
    ;;;打点
    ; res_pro=True
    ; res_pro@gsnAddCyclic = False
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnFillOn = True

    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLevels = (/0,0.01/)
    ; res_pro@lbLabelBarOn = False
    ; res_pro@cnLinesOn=False
    ; ; res_pro@cnMonoLineColor =False
    ; ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    ; res_pro@cnMonoFillPattern = True
    ; res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.003
  
    ; res_pro@cnMonoFillColor = False
    ; res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False
   
    ; res_pro@gsnCenterString = ""
    ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; ; res_pro@cnLineThicknessF = 3.
    ; ; res_pro@cnFillOn        = True
    ; iplot(i)=gsn_csm_contour(wks,prob_p,res_pro)


     ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.01     ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    ; iplot= gsn_csm_contour(wks,prob_p,res_pro)

    ; opt                      = True
    ; opt@gsnShadeFillType     = "pattern"
    ; opt@gsnShadeMid     = 6         ; stipple
    ; ; opt@gsnShadeFillDotSizeF =  0.0025
    ; ; opt@gsnShadeHigh   = "green"
    ; opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    ; iplot = gsn_contour_shade(iplot , 0,0.01,opt)
    ; overlay(plot ,iplot)
    ; print("finish subset figure")
draw(plot)
frame(wks)
end



