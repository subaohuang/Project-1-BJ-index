
;;********************************************************
;;画模式的海温趋势，并且其他海温不变
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;linshuheng 2019.10.12 
;**********************************************************
begin
;;;============读入数据===============================
    fils1 = "~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_pacifictrend.nc"
    ; fils1 = "~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc"
    
    f1 = addfile(fils1, "r")
    time = f1->time
    lon = f1->lon
    lat = f1->lat

    SST_cpl = f1->SST_cpl(:,:,:)
    SST_cpl_prediddle = f1->SST_cpl_prediddle(:,:,:)
    
    time = (/time-1/)
    date = cd_calendar(time, 1)
    timeind = ind(date.le.201212.and.date.ge.197801)
    
;;;===============计算趋势============================
    SST_ano = rmMonAnnCycTLL(SST_cpl)
    printVarSummary(SST_ano)
    
    sst_trend = regCoef_n(ispan(1,dimsizes(SST_ano&time),1),SST_ano,0,0)
    sst_trend = (sst_trend*120)
    
    copy_VarCoords(SST_ano(1,:,:), sst_trend)
  
;;;=================画图=============================
   
    ; wks = gsn_open_wks("pdf","trend_pattern_remove_atlantictrend")
    wks = gsn_open_wks("pdf","trend_pattern")
    
    res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False

    ; gsn_define_colormap(wks,"NCV_blue_red")
    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn = False

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.1
  
    res@lbOrientation = "Vertical"
    res@lbBoxEndCapStyle="TriangleBothEnds"
 
    res@tmYROn=False
    res@tmXTOn=False
    res@tmYLLabelFontHeightF=0.01
    res@tmXBLabelFontHeightF=0.01
    res@tmXBLabelFontThicknessF=0.02
    res@tmXBMinorOn=True
    ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W"/)
    ; res@tmYLLabels=(/"30~S~o~N~S","20~S~o~N~S","10~S~o~N~S","0~S~o~N","10~S~o~N~N","20~S~o~N~N","30~S~o~N~N"/)

    res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"white","transparent"/) ;填充位置相应的具体颜色

    res@mpOutlineOn = True; turn off outline 
    res@mpNationalLineThicknessF= 2.               ; double the thickness of national boundaries 
    res@mpOutlineSpecifiers = (/"Land"/)
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(sst_trend&lat)
    res@mpMaxLatF             = max(sst_trend&lat)
    res@mpMinLonF             = min(sst_trend&lon)
    res@mpMaxLonF             = max(sst_trend&lon)
    res@mpCenterLonF          =avg(sst_trend&lon)
    ; res@mpCenterLonF          = 0
   
    plot = gsn_csm_contour_map(wks, sst_trend, res)
    draw(plot)
    frame(wks)

end 






