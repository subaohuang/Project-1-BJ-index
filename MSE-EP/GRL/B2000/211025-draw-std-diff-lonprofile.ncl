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

;;;;读入控制实验结果
  modelname = "B2000_alt_north_year_CTRL"
    diri="/home/yangsong3/data-observation/linshh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".TEMP_taux_tauy.h1.0251-0280.nc")
    ; fils1 = "/home/yangsong3/data-observation/linshh/CESM-data/B/B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.TEMP_taux_tauy.h1.0251-0280.nc"
    ; fils1 = "/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/B2000_f19g16_CP_CTRL.TEMP.fvgrid.0251-0280.nc"
    f1=addfile(fils1, "r")
    time = f1->time 
    time = (/time-1/)
    date=cd_calendar(time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
    ; sst = f1->TEMP(timeind,0,:,:)  
   ; ;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    
    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})

    sst_ct = (/rmMonAnnCycTLL(sst_ct)/)

;;;读入B22000 WNPCLM海温数据
   ; modelname = "lsh_B2000_alt_north_db_year_80_20"
   modelname = "lsh_B2000_alt_north_year_80_20"
   ; modelname = "lsh_B2000_alt_north_db_year"
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
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)

    sst_wnp = sst(:,{lat1:lat2},{lon1:lon2})
    sst_wnp = (/rmMonAnnCycTLL(sst_wnp)/)


;;;;

; 
   ; nino3_ct =  dim_avg_n_Wrap(sst_ct(:,{-5:5},{210:270}), (/1,2/))   
   ; ; nino3_sen =  dim_avg_n_Wrap(sst_sen(:,{-5:5},{210:270}), (/1,2/))   
   ; nino3_wnp =  dim_avg_n_Wrap(sst_wnp(:,{-5:5},{210:270}), (/1,2/))   

   ; std1_nino3 = dim_stddev_n(nino3_ct, 0)
   ; std2_nino3 = dim_stddev_n(nino3_wnp, 0)

   ;  std1  = dim_stddev_n_Wrap(sst_ct_jja, 0)
   ;  std3  = dim_stddev_n_Wrap(sst_wnp_jja, 0)

    sst_ct_jja = month_to_annual(sst_ct, 1)
    sst_wnp_jja = month_to_annual(sst_wnp, 1)


   sst_tro_ct = dim_avg_n_Wrap(sst_ct_jja(:,{-5:5},:), (/1/))  
   sst_tro_wnp = dim_avg_n_Wrap(sst_wnp_jja(:,{-5:5},:), (/1/))  
   
   std1 = dim_stddev_n_Wrap(sst_tro_ct(:,{120:280}), 0)
   std2 = dim_stddev_n_Wrap(sst_tro_wnp(:,{120:280}), 0)

   ; s1 = var1
   ; s1 = dimsizes(sst_ct(:,1,1))
   ; s2 = s1
   ; prob = ftest(var1, s1, var2, s2, 0)
      
   std_diff = std1
   std_diff = (/std2 - std1/)

  x = std1&lon



   ; std_diff = where(prob.le.0.1, std_diff, std_diff@_FillValue)
   

  wks = gsn_open_wks("pdf", "1heating-STD-DIFF-tropical-pacific-longitude")
;;*************************************************
; create plot
;**************************************************
 

    res               = True                   ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    ; res@tiMainString  = model   ; title
    res@tiMainFontHeightF = 0.02
    res@tiYAxisString = "SST anomaly s.d."           ; y axis title
    res@tiXAxisString = ""
    

    res@tmXBLabelFontHeightF=0.01
    res@tmYLLabelFontHeightF=0.01
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0


    res@xyLineColors      = (/"black","red"/)  ; line colors
    res@xyLineThicknesses = (/4.0,4.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0/)        ; line patterns

    ; res@xyMarkLineModes   = (/"Markers","Lines"/) ; markers?
    ; res@xyMarkerColors    = (/"red",    "black"/)     ; color
    ; res@xyMarkers         = (/2,0/)          ; style

    res@tmXBMode="Explicit"

    res@trXMaxF = 280.0
    res@trXMinF = 120.0
    res@tmXTOn = False
    res@tmYROn  = False

    ; res@tmXBValues=(/0,90,180,270,357.5/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(0,360,10)
    res@vpWidthF = 0.6
    res@vpHeightF = 0.4

    res@tmXBValues=ispan(135,270,45)
    res@tmXBLabels=(/"135~S~o~N~E","180~S~o","135~S~o~N~W","90~S~o~N~W"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(120,280,15)
    res@tmYRMode="Automatic"

     plot = gsn_xy(wks,x,(/std1,std2/),res)          ; Draw an XY plot with 1 curve.
    
    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"Black","red"/) ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"CTRL","ATL warming"/)
    textres@lgLabelFontHeightF = 0.01
    ; textres@lgLabelOpacity = (/0/) 

    mid = simple_legend(wks,plot,gres,lineres,textres)   
   

     draw(plot)
     frame(wks)
end




