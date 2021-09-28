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
   lat1 = -30
   lat2 = 30
   lon1 = 40
   lon2 = 120

;;;;
 

;;;;model name 
   modelname = "B2000_alt_north_year_CTRL"
;;;;


;;;;读入控制实验结果
    diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})

    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)
    
;;;读入B22000-SST海温数据

    ;;;;model name 
     modelname = "B2000_alt_north_year"
    ;;;;
    diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")
    f1=addfile(fils1, "r")
     time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    sst_sen = sst(:,{lat1:lat2},{lon1:lon2})
    sst_sen = (/rmMonAnnCycTLL(sst_sen)/)
;;;;





;;;读入B22000 WNPCLM海温数据

    ;;;;model name 
     ; modelname = "B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer"
    ;  modelname = "lsh_B2000_alt_north_year_WNP_nowindstress_new"
    ; ;;;;
    ; diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    ; fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")

    fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.h1.0210-0340.nc"

    f1=addfile(fils1, "r")
    time := f1->time 
    time = (/time-1/)
    date:=cd_calendar(time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    sst_wnp = sst(:,{lat1:lat2},{lon1:lon2})
    sst_wnp = (/rmMonAnnCycTLL(sst_wnp)/)
;;;;

;;;
  sst_wnp_jja = month_to_season(sst_wnp,"JJA")
  ; sst_sen_jja = month_to_season(sst_sen,"JJA")

 ;;;;EOF 
    optEOF      = True
    optEOF@jopt = 0
    neval  = 30                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
    ev_sst=eofunc_n_Wrap(sst_wnp_jja,neval,optEOF,0)    ;;;得到的标准化后的场，要使得不标准化
    printVarSummary(ev_sst)
    ev_sst&lon@units="degrees_east"
    ev_sst&lat@units="degrees_north"

    ev_timef=eofunc_ts_n_Wrap(sst_wnp_jja, ev_sst, False,0)
    printVarSummary(ev_timef)

    do ne=0,neval-1
       ev_sst(ne,:,:) = ev_sst(ne,:,:)*sqrt(ev_sst@eval(ne)) ;;;使得量级能够相同
    end do

    eval=ev_sst@eval
    tsmean=ev_timef@ts_mean

    ; ts_model1 = ev_timef(0,:)/eval(0)
      ; ts_model2 = ev_timef(1,:)/eval(1)
    ; ts_model = ev_timef(0,:)
    ; ts_model=ev_timef(1,:)/eval(1)
    tos_enino1 = ev_sst(0,:,:)
    ; ts_model_std = dim_standardize(ts_model, 0)
    ; tos_enino  =  regCoef_n(ts_model_std, sst_wnp_jja, 0, 0)

    tos_enino1 = (/tos_enino1*-1/)

      tos_enino2 = ev_sst(1,:,:)
    ; ts_model_std = dim_standardize(ts_model, 0)
    ; tos_enino  =  regCoef_n(ts_model_std, sst_wnp_jja, 0, 0)

    tos_enino2 = (/tos_enino2*-1/)

;;;;


; ;;;输出数据

;     diri = "/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/EOF/CP/"


;     system("rm "+diri+"CP-EOF-ssta-djf-"+model+".nc")
;     out = addfile(diri+"CP-EOF-ssta-djf-"+model+".nc", "c")
;     out->ts_model = ts_model
;     out->pattern_model = tos_enino
;     out->pcvar = ev_sst@pcvar(0)
; ;;;


   wks = gsn_open_wks("pdf", "Indian Ocean Pattern in CTRL run")

   ;  res               = True                   ; plot mods desired
   ;  res@gsnDraw = False
   ;  res@gsnFrame = False
   ;  res@tiMainString  = "Annual mean MSE anomalies (80W-150W,0-20N)"   ; title
   ;  res@tiMainFontHeightF = 0.02
   ;  res@tiYAxisFontHeightF = 0.02
   ;  res@tiXAxisFontHeightF = 0.02
   ;  res@tiYAxisString = "MSE anomaly"           ; y axis title
   ;  res@tiXAxisString = "Year"

   ;  ; res@xyLineColors      = (/"black","red","blue","green"/)  ; line colors
   ;   res@xyLineColors      = (/"black","red","blue"/)  ; line colors
   ;   res@xyMonoLineThickness = True
   ;  res@xyLineThicknessF = (/1.0/)        ; line thicknesses
   ;  res@xyDashPatterns    = (/0.0,0.0,0.0/)        ; line patterns
   ;  res@gsnYRefLine = 0.0

   ;  ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
   ;  ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
   ;  ; res@xyMarkers         = (/2,0/)          ; style

   ;  res@tmXBMode="Explicit"

   ;  res@trXMaxF = max(x)
   ;  ; res@trXMinF = 1950
   ;  res@trXMinF = 1
   ;  ; res@trYMinF = 0.2
   ;  ; res@trYMaxF = 1.8

   ;  res@tmXTOn = False
   ;  res@tmYROn  = False

   ;  plot = gsn_csm_xy(wks, x, (/nino3_ct,nino3_sen,nino3_wnp/), res)

;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False

    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
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
    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4
    res@cnLevelSpacingF      =  0.1
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 40.0
    res@mpMaxLonF             = 120.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2

    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    

    res@gsnRightString = ev_sst@pcvar(0) +""
    

    res@gsnLeftStringFontHeightF = 0.03

    plot= new(2,graphic)
    plot(0) = gsn_csm_contour_map(wks,tos_enino1,res)
    res@gsnRightString = ev_sst@pcvar(1) +""
    plot(1) = gsn_csm_contour_map(wks,tos_enino2,res)
 

 pres = True
 gsn_panel(wks, plot, (/2,1/), pres)
end





