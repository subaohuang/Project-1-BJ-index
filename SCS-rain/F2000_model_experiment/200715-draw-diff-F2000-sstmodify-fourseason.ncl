;;********************************************************
;;将模式的海温数据的印度洋地区的海温增加，增加的幅度是通过计算sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc这个文件的趋势
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;用来跑F2000实验
;;lnshuheng 2019.12.07 
;**********************************************************
begin
season = (/"MAM","JJA","SON","DJF"/)
  fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/F/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addscs_"+season+"_notest.nc" 
  fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/F/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addindianwarming_"+season+"_notest_final.nc" 
 
  plot = new(4,graphic)
  str3 = ispan(1,4,1) + ""
   wks=gsn_open_wks("pdf", "diff-cesmsstdata2modify-add-fourseason")

  se = (/"(a)","(b)","(c)","(d)"/)

do i = 0,3

    f1 = addfile(fils1, "r")
    SST_cpl = f1->SST_cpl

    f2 = addfile(fils2(i), "r")
    SST_cpl_m = f2->SST_cpl
;;;---------------
     
;;diff 
  sst_cpl_diff = month_to_season(SST_cpl_m, season(i)) - month_to_season(SST_cpl, season(i))
;;;;
  sst_cpl_diff_re = (/sst_cpl_diff(0,:,:)/)
;;;;
  copy_VarCoords(SST_cpl_m(1,:,:), sst_cpl_diff_re)


;;; draw 
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True

    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "BlueDarkRed18"
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1

    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0


    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -40.0
    res@mpMaxLatF             = 40.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 170.0
    ; res@mpCenterLonF          = 180.0

    res@gsnLeftString  = se(i) + " "+season(i)+ " SST warming add to CAM5"
  
    plot(i) =gsn_csm_contour_map(wks, sst_cpl_diff_re, res)
end do
    pres = True
    pres@gsnPanelLabelBar = True
    pres@gsnMaximize = True
    ; pres@lbOrientation = "Vertical"
   pres@pmLabelBarHeightF = 0.15
   ;  pres@pmLabelBarWidthF = 0.5
    pres@pmLabelBarHeightF = 0.1
    gsn_panel(wks, plot, (/2,2/), pres)

end 

     






