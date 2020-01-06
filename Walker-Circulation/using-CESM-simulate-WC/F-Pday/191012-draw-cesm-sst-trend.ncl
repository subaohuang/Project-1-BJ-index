
;;********************************************************
;;画模式的海温趋势，并且其他海温不变
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;linshuheng 2019.10.12 
;**********************************************************
function draw_sst_trend(fils1,wks,gsnLeftString)
begin
;;;============读入数据===============================
    f1 = addfile(fils1, "r")
    time = f1->time
    lon = f1->lon
    lat = f1->lat

    SST_cpl = f1->SST_cpl(:,:,:)
    SST_cpl_prediddle = f1->SST_cpl_prediddle(:,:,:)
    
    time = (/time-1/)
    date = cd_calendar(time, 1)
    timeind = ind(date.le.201212.and.date.ge.197901)
    
;;;===============计算趋势============================
    SST_ano = rmMonAnnCycTLL(SST_cpl)
    printVarSummary(SST_ano)
    
    sst_trend = regCoef_n(ispan(1,dimsizes(SST_ano&time),1),SST_ano,0,0)
    sst_trend = (sst_trend*120)
    
    copy_VarCoords(SST_ano(1,:,:), sst_trend)
  
;;;=================画图=============================
   
    ; wks = gsn_open_wks("pdf","trend_pattern_remove_atlantictrend")
    

    ; colormap=read_colormap_file("MPL_RdBu")
    ; printVarSummary(colormap)
    ; colormap=(/colormap(::-1,:)/)
    ; gsn_define_colormap(wks, "cmp_b2r")
    ; gsn_define_colormap(wks, "BlueWhiteOrangeRed")
    
    res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False

    ; gsn_define_colormap(wks,"NCV_blue_red")
    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn = False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -0.5
    res@cnMaxLevelValF       =   0.5
    res@cnLevelSpacingF      =  0.05
    ; res@cnLevelSpacingF      =  0.1
    res@gsnLeftString = gsnLeftString
    res@gsnLeftStringFontHeightF = 20.0

  
    res@lbLabelBarOn = False
    ; res@lbOrientation = "Vertical"
    ; res@lbBoxEndCapStyle="TriangleBothEnds"
 
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
    ; res@mpFillDrawOrder = "PostDraw"


   
    plot = gsn_csm_contour_map(wks, sst_trend, res)
 
    ; getvalues plot@contour
    ; "cnFillColors": cnFillColors
    ; end getvalues

    ;  print(cnFillColors)
    ; dim = dimsizes(cnFillColors)
    ; ;;将中间位置的颜色设置成白色，白色的rgb值是1,1,1
    ;     NhlSetColor(wks, cnFillColors(dim/2-2),1, 1, 1)
    ;     NhlSetColor(wks, cnFillColors(dim/2-3),1, 1, 1)   

return plot
end 

begin 

wks = gsn_open_wks("pdf","rmindian_ssttrend_pattern")
   plot = new(4,graphic)

   fils = new(4,string)
   gsnLeftString = new(4,string)

   gsnLeftString(0) = "Control Run"
   gsnLeftString(1) = "remove Indian ocean warming "
   gsnLeftString(2) = "remove Atlantic warming"
   gsnLeftString(3) = "remove Pacific Lanina-like SST trend "


   fils(0) = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc"
   fils(1) = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_indiantrend.nc"
   fils(2) = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_atlantictrend.nc"
   ; fils(3) = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_pacifictrend.nc"
   fils(3) = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_1978_2012_savepacific.nc"
              
   do i = 0,3
    plot(i) = draw_sst_trend(fils((i)),wks,gsnLeftString(i))
   end do

   pres = True
   pres@gsnPanelLabelBar  = True

   gsn_panel(wks, plot, (/2,2/), pres)
   
end 






